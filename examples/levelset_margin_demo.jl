# 2-D level-set margin retreat under prescribed inward calving speed.
#
# Setup:
#   - Square Nx×Ny grid, periodic BCs (just for simplicity of indexing).
#   - Initial circular ice cap of radius R₀ at the domain centre.
#   - Zero ice flow (ū = v̄ = 0): the margin moves under calving alone.
#   - Spatially uniform inward calving speed c = 200 m/yr.
#
# Expected result: the zero contour of φ shrinks at speed c, so the radius
# decreases linearly with time, R(t) = R₀ − c·t. The script reports the
# observed radius (estimated from the area enclosed by φ ≤ 0) at three
# checkpoints and compares with this prediction.
#
# Run from the repo root:
#   julia --project examples/levelset_margin_demo.jl

using IceSheetStencils
using Symbolics
using Printf

println("=== 2-D level-set margin retreat ===\n")

spec = derive_levelset_residual()
println("Symbolic ∂φ/∂t at cell (i,j):")
println("  ", spec.dq_dt)
println()
println("Stencil offsets:")
println("  φ neighbours : ", spec.q_offsets)
println("  u (east face): ", spec.u_offsets)
println("  v (north face): ", spec.v_offsets)
println("  sources      : ", collect(zip(spec.source_kinds, spec.source_offsets)))
println()

compiled = compile_mass_stencil(spec)

# ─────────────────────────────────────────────────────────────────────────────
# Grid + initial condition
# ─────────────────────────────────────────────────────────────────────────────

Nx, Ny = 96, 96
dx, dy = 1.0e3, 1.0e3                        # 1 km
Lx = Nx * dx; Ly = Ny * dy
R0 = 30.0e3                                  # 30 km initial ice cap radius
xc, yc = Lx/2, Ly/2

xs = (collect(1:Nx) .- 0.5) .* dx
ys = (collect(1:Ny) .- 0.5) .* dy

φ = zeros(Nx, Ny)
for j in 1:Ny, i in 1:Nx
    r = sqrt((xs[i]-xc)^2 + (ys[j]-yc)^2)
    φ[i, j] = r - R0                         # signed distance to circle
end

u = zeros(Nx, Ny)        # no flow
v = zeros(Nx, Ny)
c_field = fill(200.0/(365.25*86400), Nx, Ny) # 200 m/yr inward retreat speed

mod1x(i) = mod(i-1, Nx) + 1
mod1y(j) = mod(j-1, Ny) + 1

# ─────────────────────────────────────────────────────────────────────────────
# Time stepping
# ─────────────────────────────────────────────────────────────────────────────

params = zeros(compiled.nq + compiled.nu + compiled.nv + compiled.nsrc + 2)
out = zeros(1)

function levelset_rhs!(dφ, φ, u, v, c_field)
    for j in 1:Ny, i in 1:Nx
        for k in 1:compiled.nq
            di, dj = compiled.q_offsets[k]
            params[compiled.p_q[k]] = φ[mod1x(i+di), mod1y(j+dj)]
        end
        for k in 1:compiled.nu
            di, dj = compiled.u_offsets[k]
            params[compiled.p_u[k]] = u[mod1x(i+di), mod1y(j+dj)]
        end
        for k in 1:compiled.nv
            di, dj = compiled.v_offsets[k]
            params[compiled.p_v[k]] = v[mod1x(i+di), mod1y(j+dj)]
        end
        for k in 1:compiled.nsrc                 # the only source is c
            di, dj = compiled.source_offsets[k]
            params[compiled.p_src[k]] = c_field[mod1x(i+di), mod1y(j+dj)]
        end
        params[compiled.p_scalar[1]] = dx
        params[compiled.p_scalar[2]] = dy

        compiled.rhs!(out, params)
        dφ[i, j] = out[1]
    end
end

# CFL for the calving Hamilton-Jacobi term: c·dt/dx ≤ 1.
c_max = maximum(c_field)
dt = 0.4 * dx / c_max
total_time = 100.0 * 365.25 * 86400          # 100 yr
n_steps = Int(round(total_time / dt))
@printf("dt = %.3e s (%.3f yr), %d steps over %.1f yr\n",
        dt, dt/(365.25*86400), n_steps, total_time/(365.25*86400))

dφ = zeros(Nx, Ny)

# Estimate radius from area of {φ ≤ 0}.
function radius_from_area(φ, dx, dy)
    A = sum(φ .<= 0) * dx * dy
    return sqrt(A / π)
end

# Periodically reinitialise φ to a signed distance using the same Godunov
# |∇φ| as the level-set spec — we re-evaluate the spec.dq_dt expression with
# c set to sign(φ₀) and the velocity terms zeroed by re-using the symbolic
# structure: the calving term `c·|∇φ|` of the spec already provides the
# Godunov |∇φ|, so we can compile a separate stencil with c standing in for
# `sign(φ₀)·1` and the original velocity / source paths unchanged.
#
# To keep the demo short, we use a simple, hand-rolled reinitialisation in
# pseudo-time directly on the array (5 sub-iterations every 10 real steps).
function reinit_step!(φ, φ0, dx, dy, dτ)
    Nx, Ny = size(φ)
    φnew = similar(φ)
    @inbounds for j in 1:Ny, i in 1:Nx
        Dxm = (φ[i,j] - φ[mod1x(i-1), j]) / dx
        Dxp = (φ[mod1x(i+1), j] - φ[i,j]) / dx
        Dym = (φ[i,j] - φ[i, mod1y(j-1)]) / dy
        Dyp = (φ[i, mod1y(j+1)] - φ[i,j]) / dy
        s = sign(φ0[i,j])
        if s >= 0
            gx2 = max(Dxm,0)^2 + min(Dxp,0)^2
            gy2 = max(Dym,0)^2 + min(Dyp,0)^2
        else
            gx2 = min(Dxm,0)^2 + max(Dxp,0)^2
            gy2 = min(Dym,0)^2 + max(Dyp,0)^2
        end
        g = sqrt(gx2 + gy2)
        φnew[i,j] = φ[i,j] - dτ * s * (g - 1)
    end
    φ .= φnew
end

R_init = radius_from_area(φ, dx, dy)
@printf("Initial radius (from area): %.2f km   (target %.2f km)\n", R_init/1e3, R0/1e3)

checkpoints = (n_steps ÷ 4, n_steps ÷ 2, n_steps)
let ic = 1
    for step in 1:n_steps
        levelset_rhs!(dφ, φ, u, v, c_field)
        @. φ += dt * dφ

        # Reinitialise every 10 steps.
        if step % 10 == 0
            φ0 = copy(φ)
            dτ = 0.4 * dx
            for _ in 1:5
                reinit_step!(φ, φ0, dx, dy, dτ)
            end
        end

        if step == checkpoints[ic]
            t = step * dt
            R_obs = radius_from_area(φ, dx, dy)
            R_pred = max(R0 - c_max * t, 0.0)
            @printf("  t = %5.1f yr   R_obs = %6.2f km   R_pred = %6.2f km   error = %+5.2f km\n",
                    t/(365.25*86400), R_obs/1e3, R_pred/1e3, (R_obs - R_pred)/1e3)
            ic = min(ic + 1, length(checkpoints))
        end
    end
end

println("\nDone.")
