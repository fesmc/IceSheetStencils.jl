```@meta
CurrentModule = IceSheetStencils
```

# 2-D level-set method for ice-sheet mass conservation

The package supplies symbolic stencil derivations for the two equations
introduced on the [theory page](../theory.md):

- the depth-integrated **thickness equation** for ``H`` (in two
  variants — first-order upwind and WENO5), and
- the **2-D margin level-set equation** for ``\varphi``, with a
  Godunov upwind for the advective term and an Osher–Sethian formula
  for the calving (Hamilton–Jacobi) term.

Both equations are discretised on the same Arakawa C-grid layout used
by [`derive_x_residual`](@ref) and [`derive_y_residual`](@ref): scalar
fields (``H``, ``\varphi``, ``\dot M_s``, ``\dot M_b``, ``c``) live at
cell centres, while velocity unknowns sit on cell faces — ``\bar u`` on
east faces, ``\bar v`` on north faces.

## C-grid layout reminder

```
                v_{i, j+1/2}
                     ▲
                     │
   ─────────────●────┼────●─────────────
                     │
   u_{i-1/2, j} ──▶  ●  ──▶ u_{i+1/2, j}
                     │
   ─────────────●────┼────●─────────────
                     │
                     ▼
                v_{i, j-1/2}
```

`H`, `φ`, and the source fields live at the cell centres `●`. The
existing [SSA solver](../../ssa/discretization/c_grid.md) produces
exactly this velocity layout, so the mass-conservation step consumes
its output directly.

## Thickness equation — first-order upwind

Discretising ``\partial H / \partial t + \nabla\!\cdot(\bar{\mathbf{u}}\,H)
= \dot M_s - \dot M_b`` in flux form on the C-grid gives, at cell ``(i,j)``,

```math
\frac{\partial H_{i,j}}{\partial t}
\;=\;
- \frac{F^x_{i+1/2,j} - F^x_{i-1/2,j}}{\Delta x}
- \frac{F^y_{i,j+1/2} - F^y_{i,j-1/2}}{\Delta y}
+ \dot M_{s,i,j} - \dot M_{b,i,j}.
```

The face fluxes use an upwind selection on the face velocity:

```math
F^x_{i+1/2,j}
\;=\; \begin{cases}
\bar u_{i+1/2,j}\,H_{i,j}     & \bar u_{i+1/2,j} \ge 0 \\
\bar u_{i+1/2,j}\,H_{i+1,j}   & \bar u_{i+1/2,j} <    0
\end{cases}
```

and analogously for ``F^y``. This is the formulation produced by
[`derive_thickness_residual`](@ref): a 5-point cross stencil on ``H``,
the four surrounding face velocities, and the two source fields at the
cell centre.

## Thickness equation — WENO5

For smoother solutions and lower numerical diffusion, the same flux
form is used with a fifth-order WENO reconstruction of ``H`` on each
face:

```math
H^{\mathrm{WENO5}}_{i+1/2,j}
\;=\; \omega_0\,p_0(H) + \omega_1\,p_1(H) + \omega_2\,p_2(H),
```

where ``p_k`` are three candidate third-order stencils and ``\omega_k``
are nonlinear weights driven by smoothness indicators ``\beta_k``. The
upwind side (``i-2 \dots i+2`` for ``\bar u \ge 0``,
``i-1 \dots i+3`` for ``\bar u < 0``) is selected by an `ifelse` on
``\bar u_{i+1/2,j}``. The full symbolic expression is large but
mechanical to derive — exactly the kind of expression Symbolics.jl is
designed to handle. [`derive_thickness_residual_weno5`](@ref) returns
the corresponding `Num` and a compiled in-place RHS function.

## Level-set equation — Godunov upwind plus calving

The margin equation
``\partial \varphi / \partial t + \bar{\mathbf{u}}\!\cdot\!\nabla\varphi
= c\,|\nabla\varphi|`` is discretised at each cell centre using the four
one-sided gradients

```math
D^{-x}\varphi_{i,j} = \frac{\varphi_{i,j} - \varphi_{i-1,j}}{\Delta x},\quad
D^{+x}\varphi_{i,j} = \frac{\varphi_{i+1,j} - \varphi_{i,j}}{\Delta x},
```

(and analogously in ``y``). The advective term is upwinded on the
sign of the local cell-centred velocity ``\bar u``, ``\bar v``
(averaged from the two flanking C-grid faces):

```math
\bar{\mathbf{u}} \!\cdot\! \nabla\varphi
\;\approx\;
\max(\bar u, 0)\,D^{-x}\varphi + \min(\bar u, 0)\,D^{+x}\varphi
+ \max(\bar v, 0)\,D^{-y}\varphi + \min(\bar v, 0)\,D^{+y}\varphi.
```

The calving term uses the Osher–Sethian Godunov formula for an
inward-propagating front (normal speed ``-c``, ``c \ge 0``):

```math
|\nabla\varphi|^{\mathrm{God}}_{i,j}
\;=\; \sqrt{\;
\bigl(\min(D^{-x}\varphi,0)\bigr)^2 + \bigl(\max(D^{+x}\varphi,0)\bigr)^2
+ \bigl(\min(D^{-y}\varphi,0)\bigr)^2 + \bigl(\max(D^{+y}\varphi,0)\bigr)^2
\;}.
```

Together,

```math
\frac{\partial \varphi_{i,j}}{\partial t}
\;=\;
- \bigl[\bar{\mathbf{u}}\!\cdot\!\nabla\varphi\bigr]_{i,j}
+ c_{i,j}\,|\nabla\varphi|^{\mathrm{God}}_{i,j}.
```

[`derive_levelset_residual`](@ref) assembles this expression
symbolically and returns it together with the offset metadata that
identifies which cells the stencil reads from.

## Stencil shapes

| Stencil | Scalar reads | Velocity reads (face) | Source reads |
|---|---|---|---|
| `derive_thickness_residual`        | 5-point cross on ``H`` | ``\bar u`` on east/west faces, ``\bar v`` on north/south faces | ``\dot M_s, \dot M_b`` at centre |
| `derive_thickness_residual_weno5`  | 13-point cross on ``H`` (5 in ``x``, 5 in ``y``, with ``H_{0,0}`` shared)            | as above | as above |
| `derive_levelset_residual`         | 5-point cross on ``\varphi`` | as above | ``c`` at centre |

The first-order and WENO5 thickness stencils share the same metadata
*shape* — only the symbolic flux differs — so they plug into the same
assembly loop with no other changes. The level-set stencil exposes the
same offset-list interface but operates on its own scalar field
``\varphi``.

## Reinitialisation

Reinitialisation of ``\varphi`` to a signed-distance function is done by
iterating

```math
\frac{\partial \varphi}{\partial \tau}
+ \mathrm{sign}(\varphi_0)\bigl(|\nabla\varphi|_{\mathrm{God}} - 1\bigr) = 0
```

to steady state in pseudo-time ``\tau``. The same Godunov ``|\nabla\varphi|``
formula above is used (with the sign of ``\mathrm{sign}(\varphi_0)``
selecting the appropriate one-sided gradients). A driver-level
implementation is provided in
[`examples/levelset_margin_demo.jl`](https://github.com/fesmc/IceSheetStencils.jl/blob/main/examples/levelset_margin_demo.jl);
it uses the same symbolic ``|\nabla\varphi|`` expression that
[`derive_levelset_residual`](@ref) produces, ensuring that reinit and
advection are consistent by construction.
