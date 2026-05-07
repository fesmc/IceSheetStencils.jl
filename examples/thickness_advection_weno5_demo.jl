# Fifth-order WENO discretisation of the thickness equation.
#
# This script demonstrates that Symbolics.jl can derive the (large) WENO5
# flux expression mechanically, and then compares the diffusion of an
# advected Gaussian blob between WENO5 and first-order upwind on the same
# grid and time step. With the velocity uniform and prescribed, the only
# difference between the two runs is the spatial reconstruction order.
#
# Run from the repo root:
#   julia --project examples/thickness_advection_weno5_demo.jl

using IceSheetStencils
using Symbolics
using Printf

println("=== WENO5 thickness equation ===\n")

# ─────────────────────────────────────────────────────────────────────────────
# 1. Derive both stencils symbolically
# ─────────────────────────────────────────────────────────────────────────────

spec_w = derive_thickness_residual_weno5()
spec_1 = derive_thickness_residual()

println("WENO5 stencil width: $(length(spec_w.q_offsets)) cells (vs $(length(spec_1.q_offsets)) for first-order).")
println("WENO5 H-offsets:")
println("  ", spec_w.q_offsets)

# Symbolic expression size as a rough proxy for "Symbolics is doing real work".
expr_size(ex) = length(string(ex))
println()
println("Symbolic expression length:")
@printf("  first-order upwind  : %6d chars\n", expr_size(spec_1.dq_dt))
@printf("  WENO5               : %6d chars\n", expr_size(spec_w.dq_dt))
println()

# ─────────────────────────────────────────────────────────────────────────────
# 2. Compile both and run the same advection problem with each
# ─────────────────────────────────────────────────────────────────────────────

compiled_w = compile_mass_stencil(spec_w)
compiled_1 = compile_mass_stencil(spec_1)

Nx, Ny = 64, 64
dx, dy = 1.0e3, 1.0e3
ux_const = 50.0 / (365.25 * 86400)           # 50 m/yr eastward
vy_const = 0.0

mod1x(i) = mod(i-1, Nx) + 1
mod1y(j) = mod(j-1, Ny) + 1

function make_initial()
    H = zeros(Nx, Ny)
    xs = (1:Nx) .* dx
    ys = (1:Ny) .* dy
    x0, y0, σ = Nx*dx/4, Ny*dy/2, 6*dx
    for j in 1:Ny, i in 1:Nx
        H[i, j] = 800.0 * exp(-((xs[i]-x0)^2 + (ys[j]-y0)^2)/(2σ^2))
    end
    return H
end

u = fill(ux_const, Nx, Ny)
v = fill(vy_const, Nx, Ny)
Mdot_s = zeros(Nx, Ny)
Mdot_b = zeros(Nx, Ny)

function step!(H, dH, compiled, dt)
    params = zeros(compiled.nq + compiled.nu + compiled.nv + compiled.nsrc + 2)
    out = zeros(1)
    for j in 1:Ny, i in 1:Nx
        for k in 1:compiled.nq
            di, dj = compiled.q_offsets[k]
            params[compiled.p_q[k]] = H[mod1x(i+di), mod1y(j+dj)]
        end
        for k in 1:compiled.nu
            di, dj = compiled.u_offsets[k]
            params[compiled.p_u[k]] = u[mod1x(i+di), mod1y(j+dj)]
        end
        for k in 1:compiled.nv
            di, dj = compiled.v_offsets[k]
            params[compiled.p_v[k]] = v[mod1x(i+di), mod1y(j+dj)]
        end
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
    @. H += dt * dH
end

# Advect for ~10 yr, taking 200 steps.
total_time = 10.0 * 365.25 * 86400
n_steps = 200
dt = total_time / n_steps

H1 = make_initial()
Hw = make_initial()
H0_max = maximum(H1)
dH = zeros(Nx, Ny)

@printf("Advecting for %.1f yr at dt = %.3e s (%d steps).\n",
        total_time/(365.25*86400), dt, n_steps)
println("Pre-compiling kernels...")
step!(H1, dH, compiled_1, 0.0)   # warm up build_function-d closure
step!(Hw, dH, compiled_w, 0.0)
println("Running...")

@time for _ in 1:n_steps
    step!(H1, dH, compiled_1, dt)
end
@time for _ in 1:n_steps
    step!(Hw, dH, compiled_w, dt)
end

@printf("\nPeak thickness (initial = %.2f m):\n", H0_max)
@printf("  first-order upwind  : %.2f m   (%.1f%% of initial)\n",
        maximum(H1), 100*maximum(H1)/H0_max)
@printf("  WENO5               : %.2f m   (%.1f%% of initial)\n",
        maximum(Hw), 100*maximum(Hw)/H0_max)
@printf("\nMass conservation (zero source, periodic):\n")
mass_init = sum(make_initial()) * dx * dy
@printf("  first-order upwind  Δmass / mass = %.2e\n",
        (sum(H1)*dx*dy - mass_init)/mass_init)
@printf("  WENO5               Δmass / mass = %.2e\n",
        (sum(Hw)*dx*dy - mass_init)/mass_init)

println("\nDone.")
