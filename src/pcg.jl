# Preconditioned Conjugate Gradient solvers for the SPD SSA system produced by
# the energy-functional (weak-form) discretisation.  Two variants are provided:
#
#   PCGSolver                     — textbook PCG (Hestenes–Stiefel)
#   ChronopoulosGearPCGSolver     — communication-reducing variant from
#                                   D'Azevedo & Romine (1992), used in CISM
#                                   (Lipscomb et al. 2019)
#
# Both expect a symmetric positive-definite operator.  Use
# `assemble_ssa(...; spd = true)` to obtain a system that satisfies this.

"""
    AbstractPCGSolver

Supertype for SSA-SPD Krylov solvers in this module.  Concrete subtypes carry
their own preallocated work vectors and a preconditioner, and dispatch through
[`solve!`](@ref).
"""
abstract type AbstractPCGSolver end

"""
    PCGSolver{T,P}(A; P=JacobiPreconditioner(A), tol=1e-8, maxiter=size(A,1))

Textbook preconditioned conjugate gradients (Hestenes–Stiefel) for ``Ax = b``
with `A` symmetric positive-definite.

Per iteration:

```
q  = A p
α  = (r, z) / (p, q)                ← global sync 1
x  = x + α p
r  = r − α q
z  = M⁻¹ r
β  = (r_new, z_new) / (r_old, z_old)
p  = z + β p
ρ  = ‖r‖                             ← global sync 2
```

Two global reductions per iteration.  Cheap, correct, and the reference against
which the communication-reduced variant in
[`ChronopoulosGearPCGSolver`](@ref) is validated.

The struct holds preallocated vectors so that an outer Picard / Newton loop can
re-use the solver across viscosity updates without re-allocating.
"""
struct PCGSolver{T,P<:AbstractPreconditioner} <: AbstractPCGSolver
    P::P
    tol::Float64
    maxiter::Int
    r::Vector{T}
    z::Vector{T}
    p::Vector{T}
    q::Vector{T}
end

function PCGSolver(A::AbstractMatrix;
                   P::AbstractPreconditioner = JacobiPreconditioner(A),
                   tol::Real = 1e-8,
                   maxiter::Integer = size(A, 1))
    n = size(A, 1)
    size(A, 2) == n || throw(ArgumentError("PCGSolver: matrix must be square"))
    T = float(eltype(A))
    return PCGSolver{T,typeof(P)}(P, Float64(tol), Int(maxiter),
                                  Vector{T}(undef, n), Vector{T}(undef, n),
                                  Vector{T}(undef, n), Vector{T}(undef, n))
end

"""
    ChronopoulosGearPCGSolver{T,P}(A; P=JacobiPreconditioner(A),
                                   tol=1e-8, maxiter=size(A,1))

Chronopoulos–Gear (one-sync) preconditioned conjugate gradients.

Compared with the textbook variant, the two inner products ``(r, z)`` and
``(p, Ap)`` of [`PCGSolver`](@ref) are replaced by ``\\gamma_k = (r_k, z_k)``
and ``\\delta_k = (Az_k, z_k)``, which can be merged into a **single** global
reduction per iteration.  The step lengths use the recurrence

```
β_k = γ_k / γ_{k-1}
α_k = γ_k / (δ_k − β_k γ_k / α_{k-1})
```

so that ``\\alpha_k`` is obtained from quantities already in flight rather than
from a separate ``(p_k, Ap_k)`` reduction.  This is exactly the variant
described in §3.3 of Lipscomb et al. (2019) for CISM's `cg_chronopoulos_gear`
option and originally due to D'Azevedo & Romine (1992).

The local arithmetic cost is slightly higher than textbook PCG (two extra
AXPYs to update auxiliary vectors `p, s = w + β s`), but on a distributed
machine the saving of one global reduction per iteration typically dominates.
"""
struct ChronopoulosGearPCGSolver{T,P<:AbstractPreconditioner} <: AbstractPCGSolver
    P::P
    tol::Float64
    maxiter::Int
    r::Vector{T}
    u::Vector{T}   # M⁻¹ r
    w::Vector{T}   # A u
    p::Vector{T}
    s::Vector{T}   # ≈ A p, maintained by recurrence
end

function ChronopoulosGearPCGSolver(A::AbstractMatrix;
                                   P::AbstractPreconditioner = JacobiPreconditioner(A),
                                   tol::Real = 1e-8,
                                   maxiter::Integer = size(A, 1))
    n = size(A, 1)
    size(A, 2) == n ||
        throw(ArgumentError("ChronopoulosGearPCGSolver: matrix must be square"))
    T = float(eltype(A))
    return ChronopoulosGearPCGSolver{T,typeof(P)}(
        P, Float64(tol), Int(maxiter),
        Vector{T}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n),
        Vector{T}(undef, n), Vector{T}(undef, n))
end

"""
    solve!(x, solver::AbstractPCGSolver, A, b)
        -> (; iterations, residual_norm, converged)

Solve ``Ax = b`` in place starting from the initial guess in `x` (pass
`fill!(x, 0)` for a cold start).  Returns a NamedTuple summarising the run.
"""
function solve! end

# ────────────────────────────────────────────────────────────────────────────
#  Standard PCG (Hestenes–Stiefel)
# ────────────────────────────────────────────────────────────────────────────

function solve!(x::AbstractVector, sol::PCGSolver, A::AbstractMatrix,
                b::AbstractVector)
    length(x) == length(b) == size(A, 1) ||
        throw(DimensionMismatch("solve!: x, b, A sizes inconsistent"))

    r, z, p, q = sol.r, sol.z, sol.p, sol.q

    # r = b - A x
    mul!(r, A, x)
    @. r = b - r

    bnorm = norm(b)
    target = sol.tol * (bnorm > 0 ? bnorm : one(bnorm))

    rnorm = norm(r)
    if rnorm ≤ target
        return (iterations = 0, residual_norm = rnorm, converged = true)
    end

    # z₀ = M⁻¹ r₀ ;  p₀ = z₀
    apply!(z, sol.P, r)
    copyto!(p, z)
    rz = dot(r, z)

    iter = 0
    converged = false
    for k in 1:sol.maxiter
        iter = k
        mul!(q, A, p)
        pq = dot(p, q)
        α = rz / pq

        @. x = x + α * p
        @. r = r - α * q

        rnorm = norm(r)
        if rnorm ≤ target
            converged = true
            break
        end

        apply!(z, sol.P, r)
        rz_new = dot(r, z)
        β = rz_new / rz
        @. p = z + β * p
        rz = rz_new
    end

    return (iterations = iter, residual_norm = rnorm, converged = converged)
end

# ────────────────────────────────────────────────────────────────────────────
#  Chronopoulos–Gear PCG (one-sync)
# ────────────────────────────────────────────────────────────────────────────

function solve!(x::AbstractVector, sol::ChronopoulosGearPCGSolver,
                A::AbstractMatrix, b::AbstractVector)
    length(x) == length(b) == size(A, 1) ||
        throw(DimensionMismatch("solve!: x, b, A sizes inconsistent"))

    r, u, w, p, s = sol.r, sol.u, sol.w, sol.p, sol.s

    # r₀ = b - A x₀
    mul!(r, A, x)
    @. r = b - r

    bnorm = norm(b)
    target = sol.tol * (bnorm > 0 ? bnorm : one(bnorm))

    rnorm = norm(r)
    if rnorm ≤ target
        return (iterations = 0, residual_norm = rnorm, converged = true)
    end

    # Initial step (k = 0).  γ₀ = (r₀, u₀), δ₀ = (w₀, u₀) → single reduction.
    apply!(u, sol.P, r)
    mul!(w, A, u)

    γ_prev = dot(r, u)
    δ      = dot(w, u)

    α = γ_prev / δ
    copyto!(p, u)
    copyto!(s, w)

    @. x = x + α * p
    @. r = r - α * s

    α_prev = α
    iter = 0
    converged = false
    for k in 1:sol.maxiter
        iter = k

        # Single reduction: rnorm², γ, δ
        rnorm = norm(r)
        if rnorm ≤ target
            converged = true
            break
        end

        apply!(u, sol.P, r)
        mul!(w, A, u)

        γ = dot(r, u)
        δ = dot(w, u)

        β = γ / γ_prev
        α = γ / (δ - β * γ / α_prev)

        @. p = u + β * p
        @. s = w + β * s
        @. x = x + α * p
        @. r = r - α * s

        γ_prev = γ
        α_prev = α
    end

    return (iterations = iter, residual_norm = rnorm, converged = converged)
end
