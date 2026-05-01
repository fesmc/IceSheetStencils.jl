module IceSheetStencils

using Symbolics
using SparseArrays
using LinearAlgebra

include("symbolic.jl")
include("assembly.jl")

export derive_x_residual, derive_y_residual
export StencilSpec, compile_stencil
export SSAStencils, ssa_stencils
export SSAFields, assemble_ssa

end
