# Theory

## The Shallow Shelf Approximation

The Shallow Shelf Approximation (SSA, also called the *shelfy-stream
approximation*) is the standard depth-integrated momentum balance for
ice that slides freely along its bed (ice shelves, ice streams). It is
obtained from the Stokes equations by

1. retaining only the leading-order terms in the small aspect ratio
   `ε = [thickness] / [horizontal extent]`,
2. assuming horizontal velocities are independent of depth (plug flow),
   and
3. integrating the resulting balance from bed to surface.

Two scalar equations remain — one per horizontal direction — for the
depth-averaged velocity ``\mathbf{u} = (u, v)``:

```math
\frac{\partial}{\partial x}\!\left[\,2\,\eta\, H \left(2\,\frac{\partial u}{\partial x}
                                              + \frac{\partial v}{\partial y}\right)\right]
+ \frac{\partial}{\partial y}\!\left[\,\eta\, H \left(\frac{\partial u}{\partial y}
                                              + \frac{\partial v}{\partial x}\right)\right]
- \beta\,u
\;=\; \rho\,g\,H\,\frac{\partial s}{\partial x}
```

```math
\frac{\partial}{\partial y}\!\left[\,2\,\eta\, H \left(2\,\frac{\partial v}{\partial y}
                                              + \frac{\partial u}{\partial x}\right)\right]
+ \frac{\partial}{\partial x}\!\left[\,\eta\, H \left(\frac{\partial u}{\partial y}
                                              + \frac{\partial v}{\partial x}\right)\right]
- \beta\,v
\;=\; \rho\,g\,H\,\frac{\partial s}{\partial y}
```

| Symbol | Meaning | Units |
|---|---|---|
| ``u, v``     | depth-averaged horizontal velocities       | m s⁻¹ |
| ``H``        | ice thickness                              | m |
| ``s``        | ice-surface elevation                      | m |
| ``\eta``     | effective viscosity                        | Pa·s |
| ``\beta``    | basal-friction coefficient (linear sliding)| Pa·s·m⁻¹ |
| ``\rho``     | ice density                                | kg m⁻³ |
| ``g``        | gravitational acceleration                 | m s⁻² |

## Structure of the operator

For *fixed* ``(\eta, H, \beta)`` the left-hand side is a linear, symmetric,
elliptic operator on ``(u, v)``. The right-hand side (the *driving stress*)
depends only on the geometry ``(H, s)`` and the constants ``(\rho, g)``.

- The first term in each equation is a **longitudinal (membrane) stress
  divergence** — second derivatives of the velocity component aligned with
  the equation.
- The second term is a **lateral-shear (cross-shear) stress divergence** —
  mixed second derivatives coupling ``u`` and ``v``.
- ``-\beta u``, ``-\beta v`` are the **basal-drag** terms (a linear sliding
  law). For floating ice, set ``\beta = 0``.
- The right-hand side is the **driving stress** from the surface-elevation
  gradient.

Self-adjointness of the operator implies the assembled discrete matrix is
symmetric — a property the package's tests verify on a periodic grid (see
[Validation](index.md#Validation)).

## Nonlinearity

In a full ice-flow model the effective viscosity itself depends on the
strain rate via Glen's flow law,

```math
\eta = \tfrac{1}{2}\,A^{-1/n}\,\bigl(\dot\varepsilon_e\bigr)^{(1-n)/n},
\qquad
\dot\varepsilon_e^2 = u_x^2 + v_y^2 + u_x v_y + \tfrac{1}{4}(u_y + v_x)^2,
```

with ``n \approx 3``. This makes the SSA momentum balance **nonlinear** in
``(u, v)``.

`IceSheetStencils.jl` builds the **linear-in-``(u, v)``** stencil only:
``(\eta, H, \beta)`` are treated as known fields supplied by the caller.
Outer Picard or Newton iteration on ``\eta = \eta(u, v)`` lives outside the
package; this keeps the symbolic derivation tractable and the assembled
matrix a single sparse linear operator.

## Boundary conditions

Only **periodic** boundary conditions are implemented in the current
release. The natural extensions — Dirichlet (no-slip / prescribed inflow)
and Neumann (calving-front stress balance) — will share the same
`StencilSpec` / `compile_stencil` / `assemble_ssa` machinery and are
planned for a future release. See the
[discretization overview](discretization/index.md) for the design.
