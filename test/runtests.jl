using IceSheetStencils
using LinearAlgebra
using SparseArrays
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

end
