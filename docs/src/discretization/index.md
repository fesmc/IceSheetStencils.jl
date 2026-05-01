```@meta
CurrentModule = IceSheetStencils
```

# Discretization overview

`IceSheetStencils.jl` separates the *symbolic derivation* of a discrete
stencil from its *numerical assembly* into a sparse matrix, with
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) bridging the
two. The same machinery is intended to serve every approximation we add.

## The pipeline

For each momentum-balance equation a discretisation is a three-stage
pipeline:

```
  derive_*_residual()         compile_stencil()              assemble_ssa()
 ──────────────────────► ────────────────────────────► ────────────────────────►
  symbolic residual at        compiled coefficient + RHS         sparse A, b for
   one stencil point          closures (build_function)        the periodic grid
```

The intermediate object is [`StencilSpec`](@ref): a record of

- the symbolic **residual** at one stencil centre,
- the velocity **unknowns** appearing in it, plus their grid offsets,
- the cell-centred **field parameters** (``H``, ``\eta``, ``\beta``, ``s``)
  and their offsets, and
- the **scalar parameters** ``dx``, ``dy``, ``\rho``, ``g``.

`compile_stencil` turns this into a [`CompiledStencil`](@ref) by

1. taking `Symbolics.jacobian(residual, unknowns)` to read off the linear
   stencil coefficients (one per unknown, each a function of the local
   field values), and
2. computing the RHS as ``-\mathrm{residual}\bigl|_{u=v=0}``,

then `build_function`-ing both into in-place Julia closures that can be
invoked once per grid point.

`assemble_ssa` walks the grid, gathers the local field values from the
periodic neighbourhood at each row, evaluates the compiled stencils, and
pushes triplets into a `SparseMatrixCSC`.

## Why this split?

- **Correctness**: the discrete coefficients are *derived*, not
  transcribed by hand. A change to the residual is automatically reflected
  in every entry of `A`. The [C-grid SSA page](c_grid.md) renders the
  derived coefficients directly.
- **Reuse**: a new approximation (SIA, BP, DIVA, …) only needs to provide
  its `derive_*_residual` function. `compile_stencil` and `assemble_ssa`
  are independent of the physics.
- **Performance**: the assembly loop calls compiled closures, not the
  symbolic engine. Symbolics overhead is paid once, on first use of
  [`ssa_stencils`](@ref).

## Available discretisations

| Name | Status | Description |
|---|---|---|
| [Arakawa C-grid, SSA, periodic BCs](c_grid.md) | implemented | Centred FD on a staggered grid; the only stencil currently shipped. |
| C-grid, SSA, Dirichlet / Neumann BCs | planned | Will reuse `StencilSpec` with boundary-marker logic in `assemble_ssa`. |
| C-grid, SIA / BP / DIVA            | future | Each adds a new `derive_*_residual` function; rest of pipeline is unchanged. |

If you want to inspect the existing implementation in detail, jump to the
[C-grid SSA stencil](c_grid.md).
