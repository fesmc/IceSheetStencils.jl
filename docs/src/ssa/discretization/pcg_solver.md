```@meta
CurrentModule = IceSheetStencils
```

# PCG and Chronopoulos–Gear PCG solvers

The [energy-functional formulation](energy_functional.md) gives the discrete
SSA system the structural property that makes a whole family of Krylov solvers
applicable: **symmetric positive-definiteness**. ``A`` is the Hessian of the
convex discrete energy ``E_{\mathrm{disc}}`` and is therefore SPD for strictly
positive ``\beta``.

This page describes the two preconditioned-conjugate-gradient solvers that
exploit this — the textbook (Hestenes–Stiefel) variant, and the one-sync
**Chronopoulos–Gear** variant from D'Azevedo & Romine (1992) that is used as
the higher-order solver in CISM (Lipscomb et al. 2019, §3.3).

## From the energy form to an SPD linear system

The strong-form residual ``R(u,v)`` is the gradient of ``-E_{\mathrm{disc}}``
with respect to the velocity DOFs (see the
[gradient identity](energy_functional.md#Gradient-identity)).  The Jacobian
``\partial R / \partial(u,v)`` is therefore the **negated** Hessian of
``E_{\mathrm{disc}}``, and is symmetric *negative*-definite for ``\beta > 0``.

[`assemble_ssa`](@ref) returns this strong-form operator by default, which is
fine for the existing direct (`A \ b`) path.  For CG-family solvers we want
the opposite sign — directly the Hessian of ``E_{\mathrm{disc}}`` — which is
recovered by passing `spd = true`:

```julia
A, b = assemble_ssa(SSAFields(H, η, β, s); dx, dy, spd = true)
# A is SPD: A = +∂²E_disc / ∂(u,v)²
```

This is the only difference at the assembly level.  Internally the routine
just negates both `A` and `b` once it is finished filling the triplets, so the
cost is one in-place sign flip on the nonzero array and on `b`.

## Standard PCG (Hestenes–Stiefel)

Given the SPD system ``A x = b`` and a preconditioner ``M``, textbook PCG
walks the conjugate directions

```math
\begin{aligned}
r_0 &= b - A x_0,\quad z_0 = M^{-1} r_0,\quad p_0 = z_0,\\
\text{for } k = 0, 1, \dots\\
  q_k &= A p_k,\\
  \alpha_k &= \frac{(r_k, z_k)}{(p_k, q_k)},\\
  x_{k+1} &= x_k + \alpha_k p_k,\\
  r_{k+1} &= r_k - \alpha_k q_k,\\
  z_{k+1} &= M^{-1} r_{k+1},\\
  \beta_k &= \frac{(r_{k+1}, z_{k+1})}{(r_k, z_k)},\\
  p_{k+1} &= z_{k+1} + \beta_k p_k.
\end{aligned}
```

There are **two global reductions** per iteration: one for ``(p_k, q_k)`` and
one for ``\|r_{k+1}\|`` (the convergence test).  Implemented by
[`PCGSolver`](@ref):

```julia
using IceSheetStencils

A, b = assemble_ssa(SSAFields(H, η, β, s); dx = 1e3, dy = 1e3, spd = true)
sol  = PCGSolver(A; tol = 1e-10, maxiter = 2 * size(A, 1))
x    = zeros(size(A, 1))
info = solve!(x, sol, A, b)
@show info.iterations info.residual_norm info.converged
```

The struct holds preallocated work vectors (`r, z, p, q`), so an outer
Picard iteration on ``\eta = \eta(u, v)`` can reuse the solver across
viscosity updates without re-allocating.

## Chronopoulos–Gear PCG

On parallel machines the dominant cost of CG is rarely the matrix-vector
product — it is the **global reductions** that serialise across all ranks.
The Chronopoulos–Gear variant rearranges the recurrences so that the two
inner products

```math
\gamma_k = (r_k, u_k), \qquad
\delta_k = (A u_k, u_k), \quad \text{with } u_k = M^{-1} r_k,
```

can be combined into a **single** reduction per iteration.  The step lengths
are recovered without an explicit ``(p_k, A p_k)`` inner product via the
recurrences

```math
\beta_k = \frac{\gamma_k}{\gamma_{k-1}}, \qquad
\alpha_k = \frac{\gamma_k}{\delta_k - \beta_k\,\gamma_k / \alpha_{k-1}}.
```

The full algorithm (initial step plus the loop body):

```math
\begin{aligned}
r_0 &= b - A x_0,\quad u_0 = M^{-1} r_0,\quad w_0 = A u_0,\\
\gamma_0 &= (r_0, u_0),\quad \delta_0 = (w_0, u_0),\\
\alpha_0 &= \gamma_0 / \delta_0,\\
p_0 &= u_0,\quad s_0 = w_0,\\
x_1 &= x_0 + \alpha_0 p_0,\quad r_1 = r_0 - \alpha_0 s_0,\\[4pt]
\text{for } k &= 1, 2, \dots\\
  u_k &= M^{-1} r_k,\quad w_k = A u_k,\\
  \gamma_k &= (r_k, u_k),\quad \delta_k = (w_k, u_k),\qquad\!\!\text{(single reduction)}\\
  \beta_k &= \gamma_k / \gamma_{k-1},\\
  \alpha_k &= \gamma_k / \!\left(\delta_k - \beta_k \gamma_k / \alpha_{k-1}\right),\\
  p_k &= u_k + \beta_k p_{k-1},\\
  s_k &= w_k + \beta_k s_{k-1},\\
  x_{k+1} &= x_k + \alpha_k p_k,\\
  r_{k+1} &= r_k - \alpha_k s_k.
\end{aligned}
```

The local arithmetic cost is slightly higher (two extra AXPY operations to
maintain the auxiliary ``p_k, s_k`` recurrences), but on a distributed
machine the saving of one global reduction per iteration typically dominates.
On a single workstation the runtime difference is negligible — what matters
for this package is that the two variants are *algebraically equivalent* and
produce the same iterates up to associativity of floating-point reductions.

[`ChronopoulosGearPCGSolver`](@ref) implements this with the same interface
as [`PCGSolver`](@ref):

```julia
sol  = ChronopoulosGearPCGSolver(A; tol = 1e-10)
info = solve!(x, sol, A, b)
```

The test suite checks (i) both solvers converge to the direct-solve answer to
machine precision, (ii) their iteration counts agree within a tiny tolerance,
and (iii) the SPD precondition actually holds for the operator that PCG sees
(``x^\top A x > 0`` for random nonzero ``x``).

## Preconditioner interface

The PCG solvers are written against an abstract type so that any operator
that knows how to apply ``M^{-1}`` slots in.  The minimal interface is

```julia
abstract type AbstractPreconditioner end
apply!(z::AbstractVector, P::AbstractPreconditioner, r::AbstractVector)  #  z ← M⁻¹ r
```

Two concrete preconditioners are shipped:

- [`IdentityPreconditioner`](@ref) — the no-op (``M^{-1} = I``); reduces PCG to plain CG.
  Useful as a baseline / sanity check.
- [`JacobiPreconditioner`](@ref) — the diagonal preconditioner
  ``M = \mathrm{diag}(A)``.  Cheap to apply (a single elementwise multiply),
  needs only ``2N`` of extra storage, and is the default choice — matching
  CISM's `which_ho_precond = 'diagonal'` option.

```julia
P = JacobiPreconditioner(A)            # extracts and inverts diag(A) once
sol = PCGSolver(A; P, tol = 1e-10)     # other preconditioners drop in identically
```

ILU(0) or algebraic-multigrid preconditioners would be added by defining a
new subtype of `AbstractPreconditioner` plus its `apply!` method — neither
the assembly nor the solver loops need to change.

## Why this matters for the package

The energy-functional formulation, the `spd = true` assembly option, and the
PCG solvers compose to give a complete iterative-solver path for the SSA
system that is

- **derived, not transcribed** — symmetry comes from the energy Hessian, not
  from careful per-row bookkeeping;
- **preconditioner-agnostic** — Jacobi today, ILU/AMG tomorrow, with no
  changes outside `preconditioners.jl`;
- **scalable in principle** — the Chronopoulos–Gear variant is the same
  algorithm used by CISM (Lipscomb et al. 2019) for production-scale runs,
  and is the natural target for a future MPI-parallel assembly.

## References

- Hestenes, M.R. and Stiefel, E. (1952) *Methods of conjugate gradients for
  solving linear systems*. **Journal of Research of the National Bureau of
  Standards**, 49(6), 409–436.
- D'Azevedo, E. and Romine, C.H. (1992) *Reducing communication costs in the
  conjugate gradient algorithm on distributed-memory multiprocessors*.
  Technical report ORNL/TM-12192, Oak Ridge National Laboratory.
- Chronopoulos, A.T. and Gear, C.W. (1989) *s-step iterative methods for
  symmetric linear systems*. **Journal of Computational and Applied
  Mathematics**, 25(2), 153–168.
- Lipscomb, W.H. *et al.* (2019) *Description and evaluation of the Community
  Ice Sheet Model (CISM) v2.1*. **Geoscientific Model Development**, 12,
  387–424. [doi:10.5194/gmd-12-387-2019](https://doi.org/10.5194/gmd-12-387-2019).
