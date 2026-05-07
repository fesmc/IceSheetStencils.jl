# Theory

## Mass conservation for an ice sheet

Vertically integrating the incompressible mass-conservation equation from
the bed ``z = b(x,y,t)`` to the surface ``z = s(x,y,t)`` and applying the
kinematic boundary conditions at top and bottom yields the
**depth-integrated thickness equation**

```math
\frac{\partial H}{\partial t} \;+\; \nabla\!\cdot\!\bigl(\bar{\mathbf{u}}\,H\bigr)
\;=\; \dot M_s \;-\; \dot M_b ,
```

where

| Symbol | Meaning | Units |
|---|---|---|
| ``H``                 | ice thickness                                | m |
| ``\bar{\mathbf{u}} = (\bar u, \bar v)`` | depth-averaged horizontal velocity | m s⁻¹ |
| ``\dot M_s``          | surface mass balance (accumulation − ablation) | m s⁻¹ |
| ``\dot M_b``          | basal melt rate (positive = melting)         | m s⁻¹ |

The kinematic boundary conditions at the upper and lower surfaces have
already been used to eliminate the vertical velocity ``w`` in favour of
the surface and basal mass-balance source terms. The full 3-D
englacial velocity is therefore *not* required to advance the geometry —
only its depth average ``\bar{\mathbf{u}}``, plus the two scalar source
fields ``\dot M_s`` and ``\dot M_b``.

## The lateral margin

The thickness equation above is valid wherever ice exists (``H > 0``).
The set of points where ``H = 0`` defines the **lateral margin** of the
ice sheet. Physical processes that act *only* at the margin —
calving of marine fronts, cliff collapse, lake-terminating retreat — are
naturally written as a normal speed of the margin curve, not as a
volumetric source in the thickness equation.

Several margin representations are in use in ice-sheet modelling:

- **Mask-based**: ``H``-cells flagged as ice or ocean; the margin is the
  jump between flags. Simple but introduces step-function artefacts and
  makes margin speed difficult to prescribe smoothly.
- **Volume-of-fluid / partial cells**: a cell-fraction field tracks
  partially-filled cells. Conserves mass but has well-known
  reconstruction issues at sharp corners.
- **Level set**: a smooth signed function ``\varphi(x,y,t)`` whose zero
  contour is the margin. Margin speed is prescribed directly along the
  outward normal and topology changes (calving an iceberg, joining two
  ice masses) are handled automatically.

`IceSheetStencils.jl` adopts the **2-D level-set** approach: a single
horizontal scalar ``\varphi`` tracks the planform geometry, while ``H``
is integrated only where ``\varphi \le 0`` (inside the ice mask). The
margin moves under the combined action of (i) advection by the
depth-averaged flow and (ii) an inward calving speed prescribed at the
margin.

## Margin level-set equation

Convention: ``\varphi < 0`` *inside* the ice, ``\varphi > 0`` outside,
``\varphi = 0`` at the lateral margin. The unit outward normal is
``\hat{\mathbf{n}} = \nabla\varphi / |\nabla\varphi|``.

Two effects move the margin:

1. **Advection** by the depth-averaged horizontal flow
   ``\bar{\mathbf{u}}``: ``\partial \varphi / \partial t + \bar{\mathbf{u}}\!\cdot\!\nabla\varphi = 0``.
2. **Calving** as a normal speed ``c \ge 0`` directed *inward*
   (i.e.\ in the ``-\hat{\mathbf{n}}`` direction). For an inward
   normal speed the level-set equation reads
   ``\partial \varphi / \partial t = c\,|\nabla\varphi|``.

Combining the two,

```math
\boxed{\;
\frac{\partial \varphi}{\partial t}
\;+\; \bar{\mathbf{u}}\!\cdot\!\nabla\varphi
\;=\; c\,|\nabla\varphi| .
\;}
```

The calving rate ``c`` is supplied by the caller as a cell-centred
field. Whether ``c`` comes from an analytical retreat law, an
empirical fit to observations, or a coupled hydro-fracture model is
left out of scope — the discretisation only requires its values on the
grid.

To keep ``\varphi`` close to a true signed-distance function (so that
``|\nabla\varphi|\approx 1`` near the zero set), it is periodically
**reinitialised** by iterating to steady state of

```math
\frac{\partial \varphi}{\partial \tau}
\;+\; \mathrm{sign}(\varphi_0)\,\bigl(|\nabla\varphi| - 1\bigr) \;=\; 0
```

in pseudo-time ``\tau``, with ``\varphi_0`` the un-reinitialised
field. Reinitialisation is the standard remedy for the gradient-norm
drift that level-set advection produces.

## Coupling to the thickness equation

The two scalar fields are advanced in lock-step:

1. Update ``\varphi`` by one step of the level-set equation above.
2. Update ``H`` by one step of the thickness equation, restricted to
   cells where ``\varphi \le 0``. Cells that have just been uncovered
   by a retreating margin (``\varphi`` crossed zero from negative to
   positive) are projected to ``H = 0``; cells just covered are
   initialised from upwind neighbours.
3. Every few steps, reinitialise ``\varphi`` to a signed distance.

The momentum balance (e.g.\ [SSA](../ssa/theory.md)) is solved on the
masked domain to produce ``\bar{\mathbf{u}}`` for the next step.

## Why this split is symbolic-friendly

The thickness equation is a scalar conservation law and the level-set
equation is a Hamilton–Jacobi equation. Both reduce, on a structured
grid, to **upwind flux differences** — short symbolic expressions with
conditional branches on the sign of the local velocity (or, for the
level set, on the local one-sided gradients). Symbolics.jl is well
suited to deriving these expressions and to generating the corresponding
fast Julia kernels via `build_function`. See the
[discretization page](discretization/levelset.md) for the concrete
stencils.
