module IceSheetStencils

using Symbolics
using SparseArrays
using LinearAlgebra

include("symbolic.jl")
include("assembly.jl")
include("preconditioners.jl")
include("pcg.jl")
include("mass_conservation.jl")

export derive_x_residual, derive_y_residual
export StencilSpec, compile_stencil
export SSAStencils, ssa_stencils
export SSAFields, assemble_ssa
export AbstractPreconditioner, IdentityPreconditioner, JacobiPreconditioner, apply!
export AbstractPCGSolver, PCGSolver, ChronopoulosGearPCGSolver, solve!
export derive_thickness_residual, derive_thickness_residual_weno5
export derive_levelset_residual
export MassStencilSpec, CompiledMassStencil, compile_mass_stencil

end
