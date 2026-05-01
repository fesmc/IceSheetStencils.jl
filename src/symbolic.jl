# Symbolic derivation of the SSA momentum-balance stencil on an Arakawa C-grid.
#
# C-grid layout (cell index i, j):
#   H, η, β, s  -> cell centres (i, j)
#   u[i,j]      -> east face of cell (i,j),  i.e.  u_{i+1/2, j}
#   v[i,j]      -> north face of cell (i,j), i.e.  v_{i, j+1/2}
#
# SSA strong form (linear in u,v with η,H,β treated as known fields):
#
#   ∂/∂x [ 2 η H (2 u_x + v_y) ] + ∂/∂y [ η H (u_y + v_x) ] - β u = ρ g H s_x
#   ∂/∂y [ 2 η H (2 v_y + u_x) ] + ∂/∂x [ η H (u_y + v_x) ] - β v = ρ g H s_y

"""
    StencilSpec

Symbolic specification of one momentum-balance equation discretised on a
C-grid stencil.

- `residual`     symbolic expression that should equal zero
- `unknowns`     velocity unknowns appearing in the stencil
- `unknown_offsets` list of `(kind, di, dj)` tuples telling assembly where each
                  `unknowns[k]` lives globally (`kind ∈ (:u,:v)`).
- `H_vars / H_offsets`, etc. — cell-centred parameter arrays and their offsets
  relative to the equation's own grid index.
- `s_vars`, `s_offsets`     surface-elevation parameters
- `dx, dy, ρ, g`            scalar parameters
"""
struct StencilSpec
    residual::Num
    unknowns::Vector{Num}
    unknown_offsets::Vector{Tuple{Symbol,Int,Int}}
    H_vars::Vector{Num}
    H_offsets::Vector{Tuple{Int,Int}}
    η_vars::Vector{Num}
    η_offsets::Vector{Tuple{Int,Int}}
    β_vars::Vector{Num}
    β_offsets::Vector{Tuple{Int,Int}}
    s_vars::Vector{Num}
    s_offsets::Vector{Tuple{Int,Int}}
    dx::Num
    dy::Num
    ρ::Num
    g::Num
end

"""
    derive_x_residual() -> StencilSpec

Discrete x-momentum residual at u-point `(i, j)` (= `u_{i+1/2, j}`).
Uses centred differences and 4-point averaging of `η H` to corners.
"""
function derive_x_residual()
    @variables dx dy ρ g

    # u-stencil (5-point cross around centre)
    @variables u_m1_0 u_0_0 u_p1_0 u_0_m1 u_0_p1
    # v at the four corners of the u-cell (south/north × west/east)
    @variables v_w_s v_e_s v_w_n v_e_n

    # Cell-centred fields on a 2 (x) × 3 (y) patch around the u-point:
    # west cells at i, east cells at i+1; rows j-1, j, j+1.
    @variables H_w_s H_e_s H_w_c H_e_c H_w_n H_e_n
    @variables η_w_s η_e_s η_w_c η_e_c η_w_n η_e_n
    @variables β_w_c β_e_c
    @variables s_w_c s_e_c

    # Strain-rate components at the two flanking cell centres.
    dudx_W = (u_0_0  - u_m1_0) / dx
    dudx_E = (u_p1_0 - u_0_0 ) / dx
    dvdy_W = (v_w_n - v_w_s) / dy
    dvdy_E = (v_e_n - v_e_s) / dy

    # Membrane (longitudinal) stress at cell centres → x-derivative at u-point.
    σ_W = 2 * η_w_c * H_w_c * (2*dudx_W + dvdy_W)
    σ_E = 2 * η_e_c * H_e_c * (2*dudx_E + dvdy_E)
    dσdx = (σ_E - σ_W) / dx

    # Cross-shear stress at the two y-flanking corners (i+1/2, j±1/2).
    # η H is averaged from the four surrounding cell centres.
    ηH_NE = (η_w_c*H_w_c + η_e_c*H_e_c + η_w_n*H_w_n + η_e_n*H_e_n) / 4
    ηH_SE = (η_w_c*H_w_c + η_e_c*H_e_c + η_w_s*H_w_s + η_e_s*H_e_s) / 4
    dudy_N = (u_0_p1 - u_0_0 ) / dy
    dudy_S = (u_0_0  - u_0_m1) / dy
    dvdx_N = (v_e_n - v_w_n) / dx
    dvdx_S = (v_e_s - v_w_s) / dx
    τ_N = ηH_NE * (dudy_N + dvdx_N)
    τ_S = ηH_SE * (dudy_S + dvdx_S)
    dτdy = (τ_N - τ_S) / dy

    # Basal drag and driving stress, both interpolated to the u-point.
    β_u = (β_w_c + β_e_c) / 2
    drag = β_u * u_0_0
    H_u = (H_w_c + H_e_c) / 2
    dsdx = (s_e_c - s_w_c) / dx
    drive = ρ * g * H_u * dsdx

    residual = dσdx + dτdy - drag - drive

    unknowns = [u_m1_0, u_0_0, u_p1_0, u_0_m1, u_0_p1,
                v_w_s, v_e_s, v_w_n, v_e_n]
    # u[i,j] = u_{i+1/2, j}  ⇒  u_m1_0 corresponds to u[i-1, j], etc.
    # v[i,j] = v_{i, j+1/2}  ⇒  v_w_s   corresponds to v[i,   j-1] (offset (0,-1)).
    unknown_offsets = [
        (:u, -1, 0), (:u, 0, 0), (:u, 1, 0), (:u, 0, -1), (:u, 0, 1),
        (:v, 0, -1), (:v, 1, -1), (:v, 0, 0), (:v, 1, 0),
    ]

    H_vars = [H_w_s, H_e_s, H_w_c, H_e_c, H_w_n, H_e_n]
    H_offsets = [(0,-1), (1,-1), (0,0), (1,0), (0,1), (1,1)]
    η_vars = [η_w_s, η_e_s, η_w_c, η_e_c, η_w_n, η_e_n]
    η_offsets = H_offsets
    β_vars = [β_w_c, β_e_c]
    β_offsets = [(0,0), (1,0)]
    s_vars = [s_w_c, s_e_c]
    s_offsets = [(0,0), (1,0)]

    return StencilSpec(residual,
                       unknowns, unknown_offsets,
                       H_vars, H_offsets,
                       η_vars, η_offsets,
                       β_vars, β_offsets,
                       s_vars, s_offsets,
                       dx, dy, ρ, g)
end

"""
    derive_y_residual() -> StencilSpec

Discrete y-momentum residual at v-point `(i, j)` (= `v_{i, j+1/2}`).
Mirror image of [`derive_x_residual`](@ref).
"""
function derive_y_residual()
    @variables dx dy ρ g

    # v-stencil (5-point cross around centre)
    @variables v_0_m1 v_0_0 v_0_p1 v_m1_0 v_p1_0
    # u at the four corners of the v-cell (west/east × south/north)
    @variables u_w_s u_e_s u_w_n u_e_n

    # Cell-centred fields on a 3 (x) × 2 (y) patch around the v-point:
    # south cells at j, north cells at j+1; columns i-1, i, i+1.
    @variables H_w_s H_c_s H_e_s H_w_n H_c_n H_e_n
    @variables η_w_s η_c_s η_e_s η_w_n η_c_n η_e_n
    @variables β_c_s β_c_n
    @variables s_c_s s_c_n

    # Strain rates at the two flanking cell centres (south, north of v-point).
    dvdy_S = (v_0_0  - v_0_m1) / dy
    dvdy_N = (v_0_p1 - v_0_0 ) / dy
    dudx_S = (u_e_s - u_w_s) / dx
    dudx_N = (u_e_n - u_w_n) / dx

    # Longitudinal stress at cell centres → y-derivative at v-point.
    σ_S = 2 * η_c_s * H_c_s * (2*dvdy_S + dudx_S)
    σ_N = 2 * η_c_n * H_c_n * (2*dvdy_N + dudx_N)
    dσdy = (σ_N - σ_S) / dy

    # Cross-shear stress at the two x-flanking corners (i±1/2, j+1/2).
    ηH_NE = (η_c_s*H_c_s + η_e_s*H_e_s + η_c_n*H_c_n + η_e_n*H_e_n) / 4
    ηH_NW = (η_w_s*H_w_s + η_c_s*H_c_s + η_w_n*H_w_n + η_c_n*H_c_n) / 4
    dudy_E = (u_e_n - u_e_s) / dy
    dudy_W = (u_w_n - u_w_s) / dy
    dvdx_E = (v_p1_0 - v_0_0 ) / dx
    dvdx_W = (v_0_0  - v_m1_0) / dx
    τ_E = ηH_NE * (dudy_E + dvdx_E)
    τ_W = ηH_NW * (dudy_W + dvdx_W)
    dτdx = (τ_E - τ_W) / dx

    β_v = (β_c_s + β_c_n) / 2
    drag = β_v * v_0_0
    H_v = (H_c_s + H_c_n) / 2
    dsdy = (s_c_n - s_c_s) / dy
    drive = ρ * g * H_v * dsdy

    residual = dσdy + dτdx - drag - drive

    unknowns = [v_0_m1, v_0_0, v_0_p1, v_m1_0, v_p1_0,
                u_w_s, u_e_s, u_w_n, u_e_n]
    # v[i,j] = v_{i, j+1/2} ⇒ v_0_m1 corresponds to v[i, j-1], etc.
    # u[i,j] = u_{i+1/2, j} ⇒ u_w_s  corresponds to u[i-1, j] (offset (-1, 0)).
    unknown_offsets = [
        (:v, 0, -1), (:v, 0, 0), (:v, 0, 1), (:v, -1, 0), (:v, 1, 0),
        (:u, -1, 0), (:u, 0, 0), (:u, -1, 1), (:u, 0, 1),
    ]

    H_vars = [H_w_s, H_c_s, H_e_s, H_w_n, H_c_n, H_e_n]
    H_offsets = [(-1,0), (0,0), (1,0), (-1,1), (0,1), (1,1)]
    η_vars = [η_w_s, η_c_s, η_e_s, η_w_n, η_c_n, η_e_n]
    η_offsets = H_offsets
    β_vars = [β_c_s, β_c_n]
    β_offsets = [(0,0), (0,1)]
    s_vars = [s_c_s, s_c_n]
    s_offsets = [(0,0), (0,1)]

    return StencilSpec(residual,
                       unknowns, unknown_offsets,
                       H_vars, H_offsets,
                       η_vars, η_offsets,
                       β_vars, β_offsets,
                       s_vars, s_offsets,
                       dx, dy, ρ, g)
end
