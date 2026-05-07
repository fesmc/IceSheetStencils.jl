# Symbolic derivation of mass-conservation stencils on an Arakawa C-grid.
#
# Two equations are discretised:
#
#   ∂H/∂t + ∇·(ū H) = Ṁs - Ṁb              (thickness)
#   ∂φ/∂t + ū·∇φ    = c |∇φ|                (margin level set, c ≥ 0 inward)
#
# C-grid layout (cell index i, j):
#   H, φ, Ṁs, Ṁb, c -> cell centres (i, j)
#   u[i,j]           -> east face of cell (i,j),  i.e.  ū_{i+1/2, j}
#   v[i,j]           -> north face of cell (i,j), i.e.  v̄_{i, j+1/2}
#
# The derivations return a `MassStencilSpec`: a symbolic time-derivative
# expression `dq_dt = -rhs_expr` together with the variables and their
# grid offsets relative to the equation's own cell index. `compile_mass_stencil`
# turns that into a callable closure suitable for explicit time integration.

"""
    MassStencilSpec

Symbolic specification of an explicit, scalar mass-conservation stencil at one
cell centre `(i, j)`.

The advanced field `q` is either ice thickness `H` or margin level set `φ`.
The symbolic expression `dq_dt` evaluates to `∂q/∂t` at cell `(i, j)` given
local field values.

Fields
- `dq_dt`           symbolic expression for the time derivative at `(i, j)`
- `q_vars`          neighbour values of the advected scalar appearing in the
                    stencil
- `q_offsets`       grid offsets `(di, dj)` of `q_vars` relative to `(i, j)`
- `u_vars`, `u_offsets`  C-grid east-face velocities and their offsets
                    relative to `(i, j)` (offset `(0, 0)` means `u[i, j]`,
                    which is the east face of cell `(i, j)`)
- `v_vars`, `v_offsets`  C-grid north-face velocities and their offsets
- `source_vars`,
  `source_offsets`,
  `source_kinds`    additional cell-centred fields (e.g. `Ṁs`, `Ṁb`, `c`),
                    their offsets, and a label per entry for downstream
                    bookkeeping (`:Ṁs`, `:Ṁb`, `:c`, …)
- `dx`, `dy`        scalar grid-spacing parameters
"""
struct MassStencilSpec
    dq_dt::Num
    q_vars::Vector{Num}
    q_offsets::Vector{Tuple{Int,Int}}
    u_vars::Vector{Num}
    u_offsets::Vector{Tuple{Int,Int}}
    v_vars::Vector{Num}
    v_offsets::Vector{Tuple{Int,Int}}
    source_vars::Vector{Num}
    source_offsets::Vector{Tuple{Int,Int}}
    source_kinds::Vector{Symbol}
    dx::Num
    dy::Num
end

"""
    derive_thickness_residual() -> MassStencilSpec

First-order upwind discretisation of the depth-integrated thickness equation
`∂H/∂t + ∇·(ū H) = Ṁs - Ṁb` at cell `(i, j)`. Uses a 5-point cross stencil on
`H`, the four surrounding C-grid face velocities, and the two source fields
at the cell centre.
"""
function derive_thickness_residual()
    @variables dx dy
    @variables H_m1_0 H_0_0 H_p1_0 H_0_m1 H_0_p1
    @variables u_w u_e v_s v_n            # C-grid face velocities
    @variables Mdot_s Mdot_b              # surface and basal mass-balance terms

    # Upwind face fluxes — pick the cell on the upwind side of each face.
    Fx_e = ifelse(u_e >= 0, u_e * H_0_0,  u_e * H_p1_0)
    Fx_w = ifelse(u_w >= 0, u_w * H_m1_0, u_w * H_0_0 )
    Fy_n = ifelse(v_n >= 0, v_n * H_0_0,  v_n * H_0_p1)
    Fy_s = ifelse(v_s >= 0, v_s * H_0_m1, v_s * H_0_0 )

    div_flux = (Fx_e - Fx_w)/dx + (Fy_n - Fy_s)/dy
    dq_dt = -div_flux + (Mdot_s - Mdot_b)

    q_vars     = [H_m1_0, H_0_0, H_p1_0, H_0_m1, H_0_p1]
    q_offsets  = [(-1,0), (0,0), (1,0), (0,-1), (0,1)]

    # u[i,j] is the east face of cell (i,j). The west face of (i,j) is u[i-1,j].
    u_vars     = [u_w, u_e]
    u_offsets  = [(-1, 0), (0, 0)]
    # v[i,j] is the north face of cell (i,j). The south face of (i,j) is v[i,j-1].
    v_vars     = [v_s, v_n]
    v_offsets  = [(0, -1), (0, 0)]

    source_vars    = [Mdot_s, Mdot_b]
    source_offsets = [(0, 0), (0, 0)]
    source_kinds   = [:Mdot_s, :Mdot_b]

    return MassStencilSpec(dq_dt,
                           q_vars, q_offsets,
                           u_vars, u_offsets,
                           v_vars, v_offsets,
                           source_vars, source_offsets, source_kinds,
                           dx, dy)
end

# WENO5 reconstruction of the value of `f` at a face, taking the upwind side
# of the face. `fm2, fm1, f0, fp1, fp2` are five consecutive cell-centred
# values on the upwind half-stencil ordered furthest-upwind → furthest-downwind.
#
# Standard Jiang–Shu (1996) coefficients.
function _weno5_reconstruct(fm2, fm1, f0, fp1, fp2)
    ε = 1e-6

    # Three candidate third-order stencils.
    p0 = ( 2*fm2 - 7*fm1 + 11*f0 ) / 6
    p1 = (  -fm1 +  5*f0  +  2*fp1) / 6
    p2 = ( 2*f0  +  5*fp1 -    fp2) / 6

    # Smoothness indicators.
    β0 = (13//12)*(fm2 - 2*fm1 + f0 )^2 + (1//4)*(fm2 - 4*fm1 + 3*f0)^2
    β1 = (13//12)*(fm1 - 2*f0  + fp1)^2 + (1//4)*(fm1 - fp1)^2
    β2 = (13//12)*(f0  - 2*fp1 + fp2)^2 + (1//4)*(3*f0 - 4*fp1 + fp2)^2

    # Linear weights (optimal third-stencil combination → 5th-order accuracy).
    d0, d1, d2 = 1//10, 6//10, 3//10
    α0 = d0 / (ε + β0)^2
    α1 = d1 / (ε + β1)^2
    α2 = d2 / (ε + β2)^2
    α = α0 + α1 + α2
    return (α0/α)*p0 + (α1/α)*p1 + (α2/α)*p2
end

"""
    derive_thickness_residual_weno5() -> MassStencilSpec

Fifth-order WENO discretisation of the thickness equation. The face flux at
`i+1/2, j` reconstructs `H` from the upwind side using a 5-cell WENO5 stencil;
the side is selected by the sign of the face velocity. The y-flux is built
analogously, giving a 13-point cross stencil on `H` (5 in x, 5 in y, sharing
the centre cell).
"""
function derive_thickness_residual_weno5()
    @variables dx dy
    @variables H_m2_0 H_m1_0 H_0_0 H_p1_0 H_p2_0 H_p3_0
    @variables H_0_m2 H_0_m1 H_0_p1 H_0_p2 H_0_p3
    @variables u_w u_e v_s v_n
    @variables Mdot_s Mdot_b

    # East face i+1/2: positive ū uses cells [i-2..i+2]; negative ū uses
    # cells [i-1..i+3] reversed so the WENO half-stencil still runs
    # furthest-upwind → furthest-downwind.
    H_e_pos = _weno5_reconstruct(H_m2_0, H_m1_0, H_0_0,  H_p1_0, H_p2_0)
    H_e_neg = _weno5_reconstruct(H_p3_0, H_p2_0, H_p1_0, H_0_0,  H_m1_0)
    H_e = ifelse(u_e >= 0, H_e_pos, H_e_neg)
    Fx_e = u_e * H_e

    # West face i-1/2: positive ū uses cells [i-3..i+1]; negative ū uses
    # cells [i-2..i+2] reversed.
    @variables H_m3_0
    H_w_pos = _weno5_reconstruct(H_m3_0, H_m2_0, H_m1_0, H_0_0,  H_p1_0)
    H_w_neg = _weno5_reconstruct(H_p2_0, H_p1_0, H_0_0,  H_m1_0, H_m2_0)
    H_w = ifelse(u_w >= 0, H_w_pos, H_w_neg)
    Fx_w = u_w * H_w

    # North / south faces — y direction.
    H_n_pos = _weno5_reconstruct(H_0_m2, H_0_m1, H_0_0,  H_0_p1, H_0_p2)
    H_n_neg = _weno5_reconstruct(H_0_p3, H_0_p2, H_0_p1, H_0_0,  H_0_m1)
    H_n = ifelse(v_n >= 0, H_n_pos, H_n_neg)
    Fy_n = v_n * H_n

    @variables H_0_m3
    H_s_pos = _weno5_reconstruct(H_0_m3, H_0_m2, H_0_m1, H_0_0,  H_0_p1)
    H_s_neg = _weno5_reconstruct(H_0_p2, H_0_p1, H_0_0,  H_0_m1, H_0_m2)
    H_s = ifelse(v_s >= 0, H_s_pos, H_s_neg)
    Fy_s = v_s * H_s

    div_flux = (Fx_e - Fx_w)/dx + (Fy_n - Fy_s)/dy
    dq_dt = -div_flux + (Mdot_s - Mdot_b)

    q_vars = [H_m3_0, H_m2_0, H_m1_0, H_0_0, H_p1_0, H_p2_0, H_p3_0,
              H_0_m3, H_0_m2, H_0_m1,         H_0_p1, H_0_p2, H_0_p3]
    q_offsets = [(-3,0),(-2,0),(-1,0),(0,0),(1,0),(2,0),(3,0),
                 (0,-3),(0,-2),(0,-1),       (0,1),(0,2),(0,3)]

    u_vars     = [u_w, u_e]
    u_offsets  = [(-1, 0), (0, 0)]
    v_vars     = [v_s, v_n]
    v_offsets  = [(0, -1), (0, 0)]

    source_vars    = [Mdot_s, Mdot_b]
    source_offsets = [(0, 0), (0, 0)]
    source_kinds   = [:Mdot_s, :Mdot_b]

    return MassStencilSpec(dq_dt,
                           q_vars, q_offsets,
                           u_vars, u_offsets,
                           v_vars, v_offsets,
                           source_vars, source_offsets, source_kinds,
                           dx, dy)
end

"""
    derive_levelset_residual() -> MassStencilSpec

Discretisation of the 2-D margin level-set equation
`∂φ/∂t + ū·∇φ = c |∇φ|` at cell centre `(i, j)`. Uses scalar upwinding for
the advective term and the Osher–Sethian Godunov formula for the calving
term (an inward-propagating Hamilton–Jacobi front).

Convention: `φ < 0` inside the ice; `c ≥ 0` is the inward normal speed.
"""
function derive_levelset_residual()
    @variables dx dy
    @variables φ_m1_0 φ_0_0 φ_p1_0 φ_0_m1 φ_0_p1
    @variables u_w u_e v_s v_n
    @variables c_0_0

    # Cell-centred velocities, averaged from the two flanking C-grid faces.
    ū = (u_w + u_e) / 2
    v̄ = (v_s + v_n) / 2

    # One-sided gradients.
    Dxm = (φ_0_0  - φ_m1_0) / dx
    Dxp = (φ_p1_0 - φ_0_0 ) / dx
    Dym = (φ_0_0  - φ_0_m1) / dy
    Dyp = (φ_0_p1 - φ_0_0 ) / dy

    # Scalar upwind advection: ū · ∇φ.
    advx = max(ū, 0)*Dxm + min(ū, 0)*Dxp
    advy = max(v̄, 0)*Dym + min(v̄, 0)*Dyp

    # Osher–Sethian Godunov |∇φ| for an inward-propagating front (effective
    # normal speed F = -c ≤ 0): pick min(D⁻,0)² + max(D⁺,0)² in each direction.
    gx2 = min(Dxm, 0)^2 + max(Dxp, 0)^2
    gy2 = min(Dym, 0)^2 + max(Dyp, 0)^2
    grad_god = sqrt(gx2 + gy2)

    dq_dt = -(advx + advy) + c_0_0 * grad_god

    q_vars    = [φ_m1_0, φ_0_0, φ_p1_0, φ_0_m1, φ_0_p1]
    q_offsets = [(-1,0), (0,0), (1,0), (0,-1), (0,1)]

    u_vars     = [u_w, u_e]
    u_offsets  = [(-1, 0), (0, 0)]
    v_vars     = [v_s, v_n]
    v_offsets  = [(0, -1), (0, 0)]

    source_vars    = [c_0_0]
    source_offsets = [(0, 0)]
    source_kinds   = [:c]

    return MassStencilSpec(dq_dt,
                           q_vars, q_offsets,
                           u_vars, u_offsets,
                           v_vars, v_offsets,
                           source_vars, source_offsets, source_kinds,
                           dx, dy)
end

"""
    CompiledMassStencil

In-place RHS closure plus offset metadata produced by
[`compile_mass_stencil`](@ref). The closure has the signature
`rhs!(out, params)` and writes a single scalar — the time derivative at one
cell centre — into `out[1]`. `params` is the concatenation
`[q_vars; u_vars; v_vars; source_vars; dx; dy]` in the order recorded by the
`MassStencilSpec`.
"""
struct CompiledMassStencil
    rhs!::Function
    nq::Int
    nu::Int
    nv::Int
    nsrc::Int
    q_offsets::Vector{Tuple{Int,Int}}
    u_offsets::Vector{Tuple{Int,Int}}
    v_offsets::Vector{Tuple{Int,Int}}
    source_offsets::Vector{Tuple{Int,Int}}
    source_kinds::Vector{Symbol}
    p_q::UnitRange{Int}
    p_u::UnitRange{Int}
    p_v::UnitRange{Int}
    p_src::UnitRange{Int}
    p_scalar::UnitRange{Int}   # dx, dy
end

"""
    compile_mass_stencil(spec::MassStencilSpec) -> CompiledMassStencil

`build_function`-compile the symbolic RHS of `spec` into a fast in-place
Julia closure. The compiled stencil evaluates the time derivative at one
cell centre given the local field values gathered from the grid.
"""
function compile_mass_stencil(spec::MassStencilSpec)
    nq = length(spec.q_vars)
    nu = length(spec.u_vars)
    nv = length(spec.v_vars)
    nsrc = length(spec.source_vars)
    p_q   = 1:nq
    p_u   = (nq+1):(nq+nu)
    p_v   = (nq+nu+1):(nq+nu+nv)
    p_src = (nq+nu+nv+1):(nq+nu+nv+nsrc)
    p_scalar = (nq+nu+nv+nsrc+1):(nq+nu+nv+nsrc+2)

    params = vcat(spec.q_vars, spec.u_vars, spec.v_vars, spec.source_vars,
                  Num[spec.dx, spec.dy])

    _, rhs_ip = build_function([spec.dq_dt], params; expression=Val{false})

    return CompiledMassStencil(rhs_ip, nq, nu, nv, nsrc,
                               spec.q_offsets,
                               spec.u_offsets, spec.v_offsets,
                               spec.source_offsets, spec.source_kinds,
                               p_q, p_u, p_v, p_src, p_scalar)
end
