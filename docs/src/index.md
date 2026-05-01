```@meta
CurrentModule = IceSheetStencils
```

# IceSheetStencils.jl

Symbolic derivation and sparse-matrix assembly of finite-difference stencils
for ice-sheet model approximations, built on
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

The first stencil implemented is the [Shallow Shelf Approximation
(SSA)](theory.md) momentum balance on an [Arakawa C-grid](discretization/c_grid.md)
with periodic boundary conditions. See the
[discretization overview](discretization/index.md) for the design pattern that
will be reused for additional approximations.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/fesmc/IceSheetStencils.jl")
```

## Quick start

```julia
using IceSheetStencils

Nx, Ny = 64, 48
H = fill(1000.0, Nx, Ny)             # ice thickness, m
η = fill(1e15,   Nx, Ny)             # effective viscosity, Pa·s
β = fill(1e9,    Nx, Ny)             # basal-friction coefficient, Pa·s·m⁻¹
s = [10*sin(2π*i/Nx) for i in 1:Nx, j in 1:Ny]   # surface elevation, m

A, b = assemble_ssa(SSAFields(H, η, β, s); dx = 1e3, dy = 1e3)
x = A \ b

u = reshape(x[1:Nx*Ny],     Nx, Ny)
v = reshape(x[Nx*Ny+1:end], Nx, Ny)
```

`assemble_ssa` returns the sparse linear operator `A` of size ``2N \times 2N``
with ``N = N_x N_y`` and the driving-stress right-hand side `b`. The state
vector is `[u_flat; v_flat]` with column-major flattening of each `Nx × Ny`
array.

## How it is built

1. [`derive_x_residual`](@ref) and [`derive_y_residual`](@ref) construct the
   per-stencil-point residual symbolically using `@variables` for every
   neighbour appearing in the centred-difference expression.
2. [`compile_stencil`](@ref) takes the `Symbolics.jacobian` w.r.t. the
   velocity unknowns to extract the 9 stencil coefficients, computes the
   constant (driving-stress) term as ``-\mathrm{residual}\bigl|_{u=v=0}``,
   and `build_function`s both into compiled Julia closures.
3. [`assemble_ssa`](@ref) loops over the grid, gathers local field values into
   a parameter buffer, evaluates the compiled stencil, and pushes triplets
   into a `SparseMatrixCSC`.

The compiled stencils are cached lazily on first call to [`ssa_stencils`](@ref).

## Validation

The discrete operator is validated by four test sets:

| Test | What it proves |
|---|---|
| Constant-field consistency | All viscous terms cancel for uniform velocity; drag interpolates correctly |
| Discrete adjointness | Operator is symmetric to machine precision under spatially variable ``(\eta, H, \beta)`` |
| 1-D sin(x) surface | Driving stress is wired correctly; `y`-translation symmetry preserved exactly |
| Method of manufactured solutions | Variable-coefficient discretisation converges at second order (rates 1.98, 2.00 at ``N = 16, 32, 64``) |

## License

MIT.
