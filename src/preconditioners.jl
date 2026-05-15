# Preconditioners for the SPD SSA system produced by the energy-functional
# (weak-form) discretisation.  The Krylov solvers in `pcg.jl` are written
# against the abstract interface defined here so that ILU/AMG implementations
# can be slotted in without touching the solver loops.

"""
    AbstractPreconditioner

Supertype for left preconditioners ``M^{-1}`` applied inside the PCG iteration.

A concrete subtype must implement

```julia
apply!(z::AbstractVector, P::AbstractPreconditioner, r::AbstractVector)
```

which writes ``z \\leftarrow M^{-1} r`` and returns `z`.
"""
abstract type AbstractPreconditioner end

"""
    apply!(z, P::AbstractPreconditioner, r)

Apply the preconditioner in place: ``z \\leftarrow M^{-1} r``.
"""
function apply! end

"""
    IdentityPreconditioner()

The no-op preconditioner, ``M^{-1} = I``.  Reduces PCG to plain CG.
"""
struct IdentityPreconditioner <: AbstractPreconditioner end

@inline function apply!(z::AbstractVector, ::IdentityPreconditioner,
                        r::AbstractVector)
    @. z = r
    return z
end

"""
    JacobiPreconditioner(A) <: AbstractPreconditioner
    JacobiPreconditioner(inv_diag::AbstractVector)

Diagonal (Jacobi) preconditioner, ``M = \\mathrm{diag}(A)``.  The constructor
that takes a matrix extracts and inverts the diagonal once; the constructor
that takes a vector treats it as ``1 / a_{ii}`` directly.

Matches the default `which_ho_precond = 'diagonal'` option in CISM
(Lipscomb et al. 2019); cheap to apply, requires no extra storage beyond a
single vector, and is sufficient to remove the worst conditioning from
spatially variable ``\\eta H`` and ``\\beta``.
"""
struct JacobiPreconditioner{T,V<:AbstractVector{T}} <: AbstractPreconditioner
    inv_diag::V
end

function JacobiPreconditioner(A::AbstractMatrix)
    d = diag(A)
    any(iszero, d) &&
        throw(ArgumentError("JacobiPreconditioner: matrix has zero on diagonal"))
    return JacobiPreconditioner(one.(d) ./ d)
end

@inline function apply!(z::AbstractVector, P::JacobiPreconditioner,
                        r::AbstractVector)
    @. z = P.inv_diag * r
    return z
end
