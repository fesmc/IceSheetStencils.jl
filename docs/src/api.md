# API reference

```@meta
CurrentModule = IceSheetStencils
```

## Symbolic derivation

```@docs
StencilSpec
derive_x_residual
derive_y_residual
```

## Compilation

```@docs
CompiledStencil
compile_stencil
SSAStencils
ssa_stencils
```

## Assembly

```@docs
SSAFields
assemble_ssa
```

## Preconditioners and PCG solvers

```@docs
AbstractPreconditioner
apply!
IdentityPreconditioner
JacobiPreconditioner
AbstractPCGSolver
PCGSolver
ChronopoulosGearPCGSolver
solve!
```

## Mass conservation

```@docs
MassStencilSpec
derive_thickness_residual
derive_thickness_residual_weno5
derive_levelset_residual
CompiledMassStencil
compile_mass_stencil
```

## Index

```@index
```
