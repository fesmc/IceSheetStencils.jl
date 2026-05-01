using IceSheetStencils
using LinearAlgebra
using SparseArrays
using Symbolics
using Test

const ATOL = 1e-8

@testset "IceSheetStencils.jl" begin

    @testset "stencils compile" begin
        s = ssa_stencils()
        @test s.x.nunknowns == 9
        @test s.y.nunknowns == 9
        # Cached on second call
        @test ssa_stencils() === s
    end

    @testset "constant fields, flat surface" begin
        Nx, Ny = 8, 6
        H  = fill(1000.0, Nx, Ny)
        η  = fill(1e15,   Nx, Ny)
        β0 = 1e9
        β  = fill(β0,     Nx, Ny)
        s  = fill(2000.0, Nx, Ny)
        A, b = assemble_ssa(SSAFields(H, η, β, s); dx=1e3, dy=1e3)

        # No driving stress
        @test maximum(abs, b) < ATOL

        # For uniform u with v = 0 the only surviving term in the residual is
        # the drag −β·u, so A·[1; 0] should be −β at every u-row, ≈0 elsewhere.
        N = Nx * Ny
        x_u = [ones(N); zeros(N)]
        Ax = A * x_u
        @test all(≈(-β0; atol=1e-3), Ax[1:N])
        @test maximum(abs, Ax[N+1:end]) < 1e-3

        x_v = [zeros(N); ones(N)]
        Av = A * x_v
        @test maximum(abs, Av[1:N]) < 1e-3
        @test all(≈(-β0; atol=1e-3), Av[N+1:end])
    end

    @testset "operator is symmetric on periodic grid" begin
        Nx, Ny = 8, 6
        H = [1000.0 + 100*sin(i)*cos(j) for i in 1:Nx, j in 1:Ny]
        η = [1e15 * (1 + 0.1*sin(i+j))  for i in 1:Nx, j in 1:Ny]
        β = [1e9  * (1 + 0.2*cos(i))    for i in 1:Nx, j in 1:Ny]
        s = [2000.0 + 10*sin(i+2j)      for i in 1:Nx, j in 1:Ny]
        A, _ = assemble_ssa(SSAFields(H, η, β, s); dx=1e3, dy=1e3)
        # Mathematically symmetric; floating-point rounding from differing
        # Symbolics-generated operation orders breaks strict equality.
        @test norm(A - A') < 1e-12 * norm(A)
    end

    @testset "1-D balance: surface sin(x), uniform fields ⇒ v ≡ 0" begin
        Nx, Ny = 16, 12
        H = fill(1000.0, Nx, Ny)
        η = fill(1e15,   Nx, Ny)
        β = fill(1e9,    Nx, Ny)
        # surface depends only on i ⇒ y-translation symmetry ⇒ v exactly 0
        s = [10*sin(2π*i/Nx) for i in 1:Nx, j in 1:Ny]
        A, b = assemble_ssa(SSAFields(H, η, β, s); dx=1e3, dy=1e3)

        x = A \ b
        @test norm(A*x - b) < 1e-6 * norm(b)

        N = Nx * Ny
        u = reshape(x[1:N],     Nx, Ny)
        v = reshape(x[N+1:end], Nx, Ny)
        # v must vanish; u must not depend on j
        @test maximum(abs, v) < 1e-6 * maximum(abs, u)
        for j in 2:Ny
            @test maximum(abs, u[:, j] - u[:, 1]) < 1e-8 * maximum(abs, u[:, 1])
        end
    end

    @testset "MMS convergence (variable η, H, β)" begin
        # Method of manufactured solutions: pick analytic (u, v) and variable
        # (η, H, β); compute the strong-form SSA LHS analytically; solve the
        # discrete system with that LHS as forcing; verify error decreases at
        # ~2nd order under grid refinement.

        Lx = 1e6;  Ly = 1e6
        @variables xs ys
        Dx = Differential(xs);  Dy = Differential(ys)

        # Manufactured velocity (periodic on [0, Lx] × [0, Ly]).
        u_sym = sin(2π*xs/Lx) * cos(2π*ys/Ly)
        v_sym = -cos(2π*xs/Lx) * sin(2π*ys/Ly)

        # Strictly-positive periodic coefficient fields with non-trivial variation.
        η_sym = 1e14 * (1.5 + 0.5*sin(2π*xs/Lx)*sin(2π*ys/Ly))
        H_sym = 1000.0 * (1.3 + 0.3*cos(2π*xs/Lx)*cos(2π*ys/Ly))
        β_sym = 1e9 * (1.4 + 0.4*sin(2π*xs/Lx))

        # Strong-form SSA LHS evaluated for the manufactured velocity.
        ux = expand_derivatives(Dx(u_sym))
        uy = expand_derivatives(Dy(u_sym))
        vx = expand_derivatives(Dx(v_sym))
        vy = expand_derivatives(Dy(v_sym))
        σx = 2*η_sym*H_sym*(2*ux + vy)
        σy = 2*η_sym*H_sym*(2*vy + ux)
        τ  =   η_sym*H_sym*(uy + vx)
        fx_sym = expand_derivatives(Dx(σx) + Dy(τ)) - β_sym*u_sym
        fy_sym = expand_derivatives(Dy(σy) + Dx(τ)) - β_sym*v_sym

        fx_fn = build_function(fx_sym, xs, ys; expression=Val{false})
        fy_fn = build_function(fy_sym, xs, ys; expression=Val{false})
        u_fn  = build_function(u_sym,  xs, ys; expression=Val{false})
        v_fn  = build_function(v_sym,  xs, ys; expression=Val{false})
        η_fn  = build_function(η_sym,  xs, ys; expression=Val{false})
        H_fn  = build_function(H_sym,  xs, ys; expression=Val{false})
        β_fn  = build_function(β_sym,  xs, ys; expression=Val{false})

        function solve_mms(N)
            dx = Lx/N;  dy = Ly/N
            xc(i) = (i - 0.5) * dx;  yc(j) = (j - 0.5) * dy
            xu(i) = i * dx                   # u-point: east face of cell i
            yv(j) = j * dy                   # v-point: north face of cell j

            H = [H_fn(xc(i), yc(j)) for i in 1:N, j in 1:N]
            η = [η_fn(xc(i), yc(j)) for i in 1:N, j in 1:N]
            β = [β_fn(xc(i), yc(j)) for i in 1:N, j in 1:N]
            s = zeros(N, N)            # zero out the built-in driving stress

            A, _ = assemble_ssa(SSAFields(H, η, β, s); dx=dx, dy=dy)

            b = zeros(2*N*N)
            for j in 1:N, i in 1:N
                row_u = i + (j-1)*N
                row_v = N*N + row_u
                b[row_u] = fx_fn(xu(i), yc(j))
                b[row_v] = fy_fn(xc(i), yv(j))
            end

            x = A \ b
            u_h = reshape(x[1:N*N],     N, N)
            v_h = reshape(x[N*N+1:end], N, N)
            u_ref = [u_fn(xu(i), yc(j)) for i in 1:N, j in 1:N]
            v_ref = [v_fn(xc(i), yv(j)) for i in 1:N, j in 1:N]
            return max(maximum(abs, u_h - u_ref), maximum(abs, v_h - v_ref))
        end

        Ns    = [16, 32, 64]
        errs  = [solve_mms(N) for N in Ns]
        rates = [log2(errs[k]/errs[k+1]) for k in 1:length(Ns)-1]
        @info "MMS convergence" Ns errs rates

        @test all(rates .> 1.8)              # ~2nd-order
        @test errs[end] < errs[1] / 8        # absolute error drops ≥ 8× over 4× refinement
    end

end
