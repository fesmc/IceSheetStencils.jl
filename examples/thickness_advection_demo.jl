# First-order upwind discretisation of the thickness equation.
#
# This script demonstrates:
#   1. The symbolic RHS produced by `derive_thickness_residual()` and the
#      offset metadata that locates each stencil reference on the grid.
#   2. A small time-stepping demo: a Gaussian thickness blob is advected
#      across a periodic 64×64 domain by a uniform flow, with a constant
#      surface accumulation. Total mass should grow linearly with time.
#
# Run from the repo root:
#   julia --project examples/thickness_advection_demo.jl

using IceSheetStencils
using Symbolics
using Printf

println("=== First-order upwind thickness equation ===\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Symbolic stencil
# ─────────────────────────────────────────────────────────────────────────────

spec = derive_thickness_residual()

println("Symbolic ∂H/∂t at cell (i,j):")
println("  ", spec.dq_dt)
println()
println("Stencil offsets:")
println("  H neighbours : ", spec.q_offsets)
println("  u (east face): ", spec.u_offsets)
println("  v (north face): ", spec.v_offsets)
println("  sources      : ", collect(zip(spec.source_kinds, spec.source_offsets)))
println()

# ─────────────────────────────────────────────────────────────────────────────
# 2. Compile the stencil and time-step it on a periodic grid
# ─────────────────────────────────────────────────────────────────────────────

compiled = compile_mass_stencil(spec)

Nx, Ny = 64, 64
dx, dy = 1.0e3, 1.0e3                        # 1 km
ux_const = 50.0 / (365.25 * 86400)           # 50 m/yr eastward, in m/s
vy_const = 0.0
Mdot_s_const = 0.3 / (365.25 * 86400)        # 0.3 m/yr accumulation
Mdot_b_const = 0.0

# Initial Gaussian thickness blob centred at (Nx/4, Ny/2).
xs = (1:Nx) .* dx
ys = (1:Ny) .* dy
H = zeros(Nx, Ny)
x0, y0, σ = Nx*dx/4, Ny*dy/2, 8*dx
for j in 1:Ny, i in 1:Nx
    H[i, j] = 800.0 * exp(-((xs[i]-x0)^2 + (ys[j]-y0)^2)/(2σ^2))
end

u = fill(ux_const, Nx, Ny)   # u[i,j] = ū at east face of cell (i,j)
v = fill(vy_const, Nx, Ny)   # v[i,j] = v̄ at north face of cell (i,j)
Mdot_s = fill(Mdot_s_const, Nx, Ny)
Mdot_b = fill(Mdot_b_const, Nx, Ny)

mod1x(i) = mod(i-1, Nx) + 1
mod1y(j) = mod(j-1, Ny) + 1

# Pre-allocated scratch buffers — `params` order matches `compile_mass_stencil`.
params = zeros(compiled.nq + compiled.nu + compiled.nv + compiled.nsrc + 2)
out    = zeros(1)

function rhs!(dH, H, u, v, Mdot_s, Mdot_b)
    for j in 1:Ny, i in 1:Nx
        # Gather H neighbours.
        for k in 1:compiled.nq
            di, dj = compiled.q_offsets[k]
            params[compiled.p_q[k]] = H[mod1x(i+di), mod1y(j+dj)]
        end
        # Velocities.
        for k in 1:compiled.nu
            di, dj = compiled.u_offsets[k]
            params[compiled.p_u[k]] = u[mod1x(i+di), mod1y(j+dj)]
        end
        for k in 1:compiled.nv
            di, dj = compiled.v_offsets[k]
            params[compiled.p_v[k]] = v[mod1x(i+di), mod1y(j+dj)]
        end
        # Sources, in the order recorded by source_kinds.
        for k in 1:compiled.nsrc
            kind = compiled.source_kinds[k]
            di, dj = compiled.source_offsets[k]
            field = kind === :Mdot_s ? Mdot_s : Mdot_b
            params[compiled.p_src[k]] = field[mod1x(i+di), mod1y(j+dj)]
        end
        params[compiled.p_scalar[1]] = dx
        params[compiled.p_scalar[2]] = dy

        compiled.rhs!(out, params)
        dH[i, j] = out[1]
    end
end

# CFL allows dt ≈ dx/u ≈ 10 yr for this flow; cap dt so we resolve the run
# at ~100 substeps and stay well below CFL.
n_years = 5.0
total_time = n_years * 365.25 * 86400
dt_cfl = 0.5 * dx / abs(ux_const + 1e-30)
dt = min(dt_cfl, total_time / 100)
n_steps = Int(round(total_time / dt))
println("Time stepping: dt = $(round(dt; digits=1)) s, $n_steps steps over $n_years yr.")

dH = zeros(Nx, Ny)
mass_initial = sum(H) * dx * dy
for step in 1:n_steps
    rhs!(dH, H, u, v, Mdot_s, Mdot_b)
    @. H += dt * dH
    if step == 1 || step == n_steps ÷ 2 || step == n_steps
        mass = sum(H) * dx * dy
        @printf("  step %5d  max H = %7.2f m   total mass = %.4e m³\n",
                step, maximum(H), mass)
    end
end

mass_final = sum(H) * dx * dy
domain_area = Nx*dx * Ny*dy
expected_growth = Mdot_s_const * (n_steps * dt) * domain_area
mass_growth = mass_final - mass_initial

println()
@printf("Mass-balance check (uniform Ṁs over periodic domain):\n")
@printf("  observed Δmass = %.4e m³\n", mass_growth)
@printf("  expected Δmass = %.4e m³\n", expected_growth)
@printf("  relative error = %.3e\n", abs(mass_growth - expected_growth)/expected_growth)

println("\nDone.")
