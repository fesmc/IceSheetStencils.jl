# IceSheetStencils.jl

Symbolic derivation and sparse-matrix assembly of finite-difference stencils
for ice-sheet model approximations, built on
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## SSA momentum balance

The first stencil implemented is the [Shallow Shelf
Approximation](https://en.wikipedia.org/wiki/Ice-sheet_dynamics#Shallow_shelf_approximation)
discretised on an Arakawa C-grid:

```
∂ₓ[ 2 η H (2 uₓ + v_y) ] + ∂_y[ η H (u_y + vₓ) ] − β u = ρ g H sₓ
∂_y[ 2 η H (2 v_y + uₓ) ] + ∂ₓ[ η H (u_y + vₓ) ] − β v = ρ g H s_y
```

with periodic boundary conditions. The system is linear in `(u, v)` for fixed
`(η, H, β)`; the outer Picard / Newton iteration on `η = η(u, v)` lives in the
caller. Each row carries 9 nonzeros (5-point cross in the same component plus
4 corner couplings to the other component).

## Usage

```julia
using IceSheetStencils

Nx, Ny = 64, 48
H = fill(1000.0, Nx, Ny)
η = fill(1e15,   Nx, Ny)
β = fill(1e9,    Nx, Ny)
s = [10*sin(2π*i/Nx) for i in 1:Nx, j in 1:Ny]

A, b = assemble_ssa(SSAFields(H, η, β, s); dx=1e3, dy=1e3)
x = A \ b
u = reshape(x[1:Nx*Ny],     Nx, Ny)
v = reshape(x[Nx*Ny+1:end], Nx, Ny)
```

`assemble_ssa` returns the sparse operator `A` (`2 N × 2 N`, with `N = Nx*Ny`)
and the driving-stress RHS `b`. The state vector is `[u_flat; v_flat]` with
column-major flattening.

## How the stencil is derived

1. `derive_x_residual()` and `derive_y_residual()` build the symbolic residual
   at one stencil point, using `@variables` for every neighbour appearing in
   the centred-difference expression.
2. `compile_stencil(spec)` takes the `Symbolics.jacobian` w.r.t. the velocity
   unknowns to extract the 9 stencil coefficients, computes the constant
   (driving-stress) term as `−residual|_{u=v=0}`, and `build_function`s both
   into compiled Julia closures.
3. `assemble_ssa` loops over the grid, gathers local field values into a
   parameter buffer, evaluates the compiled stencil, and pushes triplets into
   a `SparseMatrixCSC`.

The compiled stencils are cached lazily on first call to `ssa_stencils()`.
