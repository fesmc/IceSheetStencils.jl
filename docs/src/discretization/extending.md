```@meta
CurrentModule = IceSheetStencils
```

# Extending to a realistic solver

!!! note "Design notes"
    Nothing on this page is shipped yet. The
    [C-grid (SSA)](c_grid.md) stencil currently supports periodic boundary
    conditions only. This page collects the design we'd follow to grow it
    into a solver that handles realistic ice-sheet geometry — Dirichlet
    inflows, calving fronts, and time-varying ice / no-ice masks — without
    redesigning what already works.

## What's needed for a realistic SSA solver

Beyond what the package ships today, a production-grade SSA solver has to
handle three things:

1. **Lateral domain BCs** — Dirichlet (no-slip walls, prescribed inflow)
   or Neumann (free-slip / prescribed traction).
2. **Calving / ice-front traction** — at the seaward boundary of an ice
   shelf, a depth-integrated water-pressure stress balances the
   longitudinal stress.
3. **Internal ice / no-ice masks** — most ice sheets sit inside a
   computational domain that contains open ocean and bare bedrock cells.
   The momentum balance is degenerate where ``H = 0`` and has to be
   excluded from the linear system.

These split into two tiers: the first two slot into the existing pipeline
with little disturbance; the third is the one that needs actual design
work.

## Tier 1 — natural extensions

### Dirichlet BCs

Replace the row at a constrained `u`- or `v`-point with the trivial
equation ``u = u_{\mathrm{BC}}``: zero the row, set the diagonal to one,
and put ``u_{\mathrm{BC}}`` into the corresponding entry of `b`.
Symmetry is preserved by **also** zeroing the column and folding the
displaced product ``A_{ij}\,u_{\mathrm{BC}}`` into `b` at every row `i`
that coupled to the Dirichlet point. This is a few dozen lines added to
[`assemble_ssa`](@ref) once a `dirichlet_mask` is plumbed through.

The current C-grid stencil itself is untouched.

### Calving-front traction

At a vertical ice front the boundary integral

```math
\int_{-d}^{s} \bigl( \sigma_{xx}\,n_x + \sigma_{xy}\,n_y \bigr)\,\mathrm{d}z
\;=\;
\tfrac{1}{2} \rho_i g H^2 - \tfrac{1}{2} \rho_w g d^2
```

(with `d` the submerged thickness and ``n`` the outward normal) is the
hydrostatic-difference stress that must balance the longitudinal stress
at the front. Crucially, this term is

- **linear in the geometry ``(H, d)``** — i.e. it depends on the *fields*,
  not on the velocity unknowns, so it contributes to ``b`` and never to
  ``A``;
- **localised to faces between an ice cell and an ice-free cell** — so
  it's a per-face addition, not a per-row one.

In practice we'd derive it once symbolically alongside the existing
driving stress (it has the same structure: a discrete divergence of an
``H``-quadratic), compile it via `build_function`, and add its
contribution into ``b`` for every u-/v-point flanked by an ocean cell.

The LHS — and therefore everything in [`compile_stencil`](@ref) — is
unchanged.

## Tier 2 — internal ice / no-ice masks

This is the only piece that requires real design work, because the
current stencil's averaging assumptions break down at the ice/no-ice
interface:

- The cross-shear corner average ``\eta H`` mixes the four cells around a
  corner; if one of them is ice-free, the average no longer represents a
  meaningful stress.
- The membrane stress ``\partial_x[2\eta H(2 u_x + v_y)]`` is
  ill-posed where ``H = 0`` on one side of the face.

There are two clean ways to handle this.

### Option (a): pruned-domain solve

Build the global matrix as today (treating no-ice cells as if they had
some nominal ``H, \eta``), then **delete the rows and columns** that
correspond to inactive `(u, v)` degrees of freedom and add the
calving-front traction (Tier 1) to ``b`` for the surviving rows whose
neighbour was masked out. Concretely this means:

1. add an `ice_mask::Matrix{Bool}` field to [`SSAFields`](@ref);
2. derive an `active_dofs(ice_mask)` index that selects only u-/v-points
   between two ice cells;
3. run today's [`assemble_ssa`](@ref) at full size, then take
   `A_act = A[active, active]` and `b_act = b[active]`;
4. add the front-traction contribution to `b_act` for u-/v-points whose
   *would-have-been* neighbour was inactive.

This treats every ice/ocean face automatically as a free boundary in the
discrete sense, which is exactly the calving-front condition. It re-uses
**the entire current pipeline** — one stencil, one assembly loop, one
sparse matrix — and is the approach taken by ISSM, Úa, and CISM.

```
ice_mask         build full A,b      restrict to active     add front traction
   ┌───┐  ────►  (today)        ────► A_act, b_act      ────► to b_act on faces
   │ ▮ │                                                       touching no-ice
   └───┘
```

### Option (b): per-pattern stencils

Pre-derive a [`StencilSpec`](@ref) for each distinct mask configuration
around the centre cell and dispatch by reading the local 3×3 mask in the
assembly loop. With ice/ocean only, the 9 surrounding cells give ``2^9 =
512`` possible patterns, but symmetry collapses this to roughly 8–16
distinct stencils in practice (interior, +x face, −x face, NE corner,
etc.). Each one would re-use [`compile_stencil`](@ref) verbatim.

This is more rigorous than (a) — every stencil is the *correct* discrete
operator for its local geometry — but it multiplies the code, the test
matrix, and the symbolic-derivation time accordingly.

### Recommendation

Take **Option (a)**. It composes with what's already shipped, mirrors how
the established ice-sheet codes do it, and the pruning step is a few
lines on top of the existing assembly. Option (b) is worth keeping in
mind if a downstream user finds (a) unacceptably diffusive across a
particular boundary — but that's a much later concern.

## What changes to the public API

| Today | After Tier 1 + Option (a) |
|---|---|
| `SSAFields(H, η, β, s)` | `SSAFields(H, η, β, s; ice_mask = trues(...), dirichlet = nothing)` |
| `assemble_ssa(fields; dx, dy, ρ, g)` returns `(A, b)` of size ``2N \times 2N`` | same call, returns the *pruned* `(A, b)` of size ``2 N_{\mathrm{active}}`` plus an index map so the caller can lift the solution back onto the full grid |
| `derive_x_residual` / `derive_y_residual` | unchanged |
| `compile_stencil`, `ssa_stencils`     | unchanged |

The Symbolics-driven core — derivation, Jacobian extraction,
`build_function`-compiled stencils — is exactly what we have now. The
extension is concentrated in `assemble_ssa` and in two small new helpers
(boundary-traction RHS, mask-to-active-dofs).

## What stays the same

- The single-stencil [`compile_stencil`](@ref) machinery.
- The [validation suite](../index.md#Validation): every test stays
  meaningful — periodic-BC consistency, discrete adjointness, MMS
  convergence — and the new boundary code adds its own small tests on top
  (e.g. an MMS in a rectangular domain with prescribed inflow).
- The contract that the LHS is **linear in `(u, v)` for fixed
  `(\eta, H, \beta)`**, with the outer Picard / Newton iteration on
  ``\eta = \eta(u, v)`` living above the package.
