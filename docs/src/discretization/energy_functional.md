```@meta
CurrentModule = IceSheetStencils
```

# Variational (energy-functional) formulation

The SSA momentum balance is the Euler–Lagrange equation of a **quadratic energy
functional** ``E[u, v]``.  Deriving the discrete stencil from ``E`` rather than
directly from the strong form has two structural advantages.

- **Symmetry is guaranteed.** The stiffness matrix equals the (negated) Hessian
  of a quadratic, so ``A = A^\top`` follows from algebra, not from careful
  bookkeeping.
- **Boundary conditions are local energy modifications.** Each BC type — Dirichlet
  inflow, calving-front traction, or stress-free wall — is an additive term in
  ``E`` at the relevant boundary faces.  The interior stencil is untouched.

## The continuous functional

Define the depth-integrated energy density

```math
W(u,v) \;=\;
  \eta H \!\left[2u_x^2 + 2v_y^2 + 2\,u_x v_y + \tfrac{1}{2}(u_y + v_x)^2\right]
  + \tfrac{1}{2}\beta\!\left(u^2 + v^2\right)
  + \rho\,g\,H\!\left(u\,\frac{\partial s}{\partial x}
                     + v\,\frac{\partial s}{\partial y}\right)
```

and the total energy

```math
E[u,v] \;=\; \iint_\Omega W(u,v)\;\mathrm{d}A.
```

!!! note "Sign of the driving-stress term"
    The driving-stress contribution enters with a **positive** sign in ``W``.
    ``E`` is convex (positive-definite quadratic in the strain rates plus the
    basal drag) and is *minimised* at the physical solution.  The positive sign
    is what makes ``\delta E / \delta u = 0`` produce the correct SSA equation
    rather than its negative — see the derivation below.

## Euler–Lagrange equations

Setting the first variation of ``E`` to zero gives (after integration by parts
over a domain with periodic or zero-flux boundary conditions)

```math
\frac{\delta E}{\delta u} \;=\;
\beta u + \rho g H \frac{\partial s}{\partial x}
- \frac{\partial}{\partial x}\!\Bigl[2\eta H\!\left(2u_x + v_y\right)\Bigr]
- \frac{\partial}{\partial y}\!\Bigl[\eta H\!\left(u_y + v_x\right)\Bigr]
\;=\; 0.
```

Multiplying by ``-1`` recovers the SSA x-momentum equation:

```math
\frac{\partial}{\partial x}\!\Bigl[2\eta H\!\left(2u_x + v_y\right)\Bigr]
+ \frac{\partial}{\partial y}\!\Bigl[\eta H\!\left(u_y + v_x\right)\Bigr]
- \beta u \;=\; \rho g H \frac{\partial s}{\partial x}.
```

The y-momentum equation follows from ``\delta E / \delta v = 0`` by the
``x \leftrightarrow y``, ``u \leftrightarrow v`` symmetry of ``W``.

## The discrete energy on the C-grid

The continuous integral is replaced by a sum over the four types of contribution
that arise from the C-grid layout.

### Membrane terms (at cell centres)

At each cell centre ``(i,j)`` the longitudinal strain rates are approximated
with the flanking face velocities

```math
\dot\varepsilon_{xx}^{i,j} = \frac{u_{i,j} - u_{i-1,j}}{dx},
\qquad
\dot\varepsilon_{yy}^{i,j} = \frac{v_{i,j} - v_{i,j-1}}{dy},
```

and the membrane energy contribution is

```math
E_{\mathrm{mem}}^{i,j} \;=\;
  \eta H\big|_{i,j}
  \!\left(2\bigl(\dot\varepsilon_{xx}^{i,j}\bigr)^{\!2}
         + 2\bigl(\dot\varepsilon_{yy}^{i,j}\bigr)^{\!2}
         + 2\,\dot\varepsilon_{xx}^{i,j}\,\dot\varepsilon_{yy}^{i,j}\right)
  dx\,dy.
```

### Shear terms (at corners)

At each cell corner ``(i+\tfrac12, j+\tfrac12)`` the shear strain rate is

```math
\dot\varepsilon_{xy}\big|_{i+\frac12,j+\frac12} \;=\;
  \tfrac{1}{2}\!\left(\frac{u_{i,j+1} - u_{i,j}}{dy}
                     + \frac{v_{i+1,j} - v_{i,j}}{dx}\right),
```

with ``\eta H`` averaged from the four surrounding cell centres,

```math
\widetilde{\eta H}\big|_{i+\frac12,j+\frac12} \;=\;
  \tfrac{1}{4}\!\left(\eta H\big|_{i,j} + \eta H\big|_{i+1,j}
                    + \eta H\big|_{i,j+1} + \eta H\big|_{i+1,j+1}\right),
```

giving

```math
E_{\mathrm{shear}}^{i+\frac12,j+\frac12} \;=\;
  2\,\widetilde{\eta H}\big|_{i+\frac12,j+\frac12}
  \bigl(\dot\varepsilon_{xy}\big|_{i+\frac12,j+\frac12}\bigr)^{\!2}
  \;dx\,dy.
```

### Drag and driving-stress terms (at velocity DOFs)

At each u-point ``u_{i+\frac12,j}``, with ``\beta`` and ``H`` averaged from the
two flanking cells ``(i,j)`` and ``(i+1,j)``:

```math
E_{\mathrm{drag},u}^{i+\frac12,j} \;=\;
  \tfrac{1}{2}\bar\beta\; u_{i+\frac12,j}^2\; dx\,dy,
\qquad
E_{\mathrm{drive},u}^{i+\frac12,j} \;=\;
  \rho g \bar H\; u_{i+\frac12,j}\;\frac{s_{i+1,j} - s_{i,j}}{dx}\; dx\,dy.
```

The v-point contributions are the ``x \leftrightarrow y``, ``u \leftrightarrow v``
mirror.

## Gradient identity

Taking the partial derivative of ``E_{\mathrm{disc}}`` with respect to the
u-point ``u_{i+\frac12,j}`` gives

```math
\frac{\partial E_{\mathrm{disc}}}{\partial u_{i+\frac12,j}}
\;=\; -dx\,dy \;\times\; R_x(i,j),
```

where ``R_x(i,j)`` is the strong-form residual from [`derive_x_residual`](@ref)
evaluated at the local field values.  Setting ``\partial E / \partial u = 0`` is
therefore exactly equivalent to setting the residual to zero: **both formulations
produce the same linear system** ``A\,x = b``.

The stiffness matrix satisfies ``A = -\partial^2 E_{\mathrm{disc}} / \partial
x^2 \,/\, (dx\,dy)``, i.e. it is the negated (and area-normalised) Hessian of a
convex quadratic.  Symmetry follows immediately.

The code below reconstructs the local energy patch with the same Symbolics.jl
variables used by [`derive_x_residual`](@ref) and verifies the gradient identity
numerically.

```@example enfunc
using IceSheetStencils, Symbolics

# Symbols shared with derive_x_residual — same names → same Symbolics objects.
@variables dx dy ρ g
@variables u_m1_0 u_0_0 u_p1_0 u_0_m1 u_0_p1
@variables v_w_s v_e_s v_w_n v_e_n
@variables H_w_s H_e_s H_w_c H_e_c H_w_n H_e_n
@variables η_w_s η_e_s η_w_c η_e_c η_w_n η_e_n
@variables β_w_c β_e_c
@variables s_w_c s_e_c

# ── Membrane: west cell (i,j) and east cell (i+1,j) ──────────────────────────
dudx_W = (u_0_0 - u_m1_0) / dx;   dvdy_W = (v_w_n - v_w_s) / dy
dudx_E = (u_p1_0 - u_0_0) / dx;   dvdy_E = (v_e_n - v_e_s) / dy

E_mem = η_w_c*H_w_c*(2*dudx_W^2 + 2*dvdy_W^2 + 2*dudx_W*dvdy_W)*dx*dy +
        η_e_c*H_e_c*(2*dudx_E^2 + 2*dvdy_E^2 + 2*dudx_E*dvdy_E)*dx*dy

# ── Shear: south corner (i+½,j-½) and north corner (i+½,j+½) ────────────────
ηH_S = (η_w_c*H_w_c + η_e_c*H_e_c + η_w_s*H_w_s + η_e_s*H_e_s) / 4
ηH_N = (η_w_c*H_w_c + η_e_c*H_e_c + η_w_n*H_w_n + η_e_n*H_e_n) / 4
dudy_S = (u_0_0 - u_0_m1) / dy;   dvdx_S = (v_e_s - v_w_s) / dx
dudy_N = (u_0_p1 - u_0_0) / dy;   dvdx_N = (v_e_n - v_w_n) / dx

E_shear = ηH_S*(1//2)*(dudy_S + dvdx_S)^2*dx*dy +
          ηH_N*(1//2)*(dudy_N + dvdx_N)^2*dx*dy

# ── Drag and drive: averaged to the u-point ───────────────────────────────────
E_drag  = (1//2)*((β_w_c + β_e_c)/2)*u_0_0^2*dx*dy
E_drive = ρ*g*((H_w_c + H_e_c)/2)*u_0_0*((s_e_c - s_w_c)/dx)*dx*dy

E_patch = E_mem + E_shear + E_drag + E_drive
nothing # hide
```

```@example enfunc
# Gradient ∂E/∂u₀
∂E_∂u = expand_derivatives(Differential(u_0_0)(E_patch))

# Reference residual from the strong-form derivation
spec = derive_x_residual()

# Verify numerically: substitute a concrete set of field values
vals = Dict(dx => 1000.0, dy => 1000.0, ρ => 910.0, g => 9.81,
            u_m1_0 => 100.0, u_0_0 => 120.0, u_p1_0 => 145.0,
            u_0_m1 => 108.0, u_0_p1 => 133.0,
            v_w_s => -8.0,  v_e_s => -6.0,  v_w_n => 11.0, v_e_n => 13.0,
            H_w_s => 950.0, H_e_s => 970.0,
            H_w_c => 1000.0, H_e_c => 1020.0,
            H_w_n => 1050.0, H_e_n => 1080.0,
            η_w_s => 1e14, η_e_s => 1.1e14,
            η_w_c => 1.2e14, η_e_c => 1.3e14,
            η_w_n => 1.4e14, η_e_n => 1.5e14,
            β_w_c => 1e9,   β_e_c => 1.2e9,
            s_w_c => 2010.0, s_e_c => 1990.0)

sub(ex) = Float64(Symbolics.value(Symbolics.substitute(ex, vals)))

grad_norm = sub(∂E_∂u) / (1000.0 * 1000.0)
res_val   = sub(spec.residual)
println("∂E / (dx·dy·∂u) = ", grad_norm)
println("      −residual  = ", -res_val)
println("   relative diff = ", abs(grad_norm + res_val) / abs(res_val))
```

## Boundary conditions via energy modification

The key insight is that every BC type is an **additive modification to ``E``**
at the boundary faces — the interior energy expression and the stencil coefficients
are left intact.

### Stress-free (natural Neumann) BC

No modification needed.  When a domain boundary is present and no boundary term
is added, integration by parts leaves the natural condition

```math
\bigl[\,2\eta H(2u_x + v_y)\,\bigr]_{\partial\Omega}\cdot n_x = 0,
```

i.e. zero longitudinal stress at the open face.  This is the correct default for
a freely calving or free-slip wall.

### Prescribed traction (calving-front BC)

At an ice–ocean front with outward normal ``\hat{x}``, the depth-integrated
stress balance gives the traction

```math
T_{\mathrm{front}} \;=\; \tfrac{1}{2}\rho_i g H^2 - \tfrac{1}{2}\rho_w g d^2,
```

where ``d`` is the submerged thickness.  Imposing this traction is equivalent to
adding the boundary work integral

```math
E_{\mathrm{front}} \;=\;
  -\int_{\partial\Omega_{\mathrm{front}}} u\; T_{\mathrm{front}}\;\mathrm{d}\ell
```

to the total energy.  Because ``E_{\mathrm{front}}`` is **linear** in ``u``, its
gradient adds only to the right-hand side vector ``b`` — the stiffness matrix
``A`` is unchanged.

```@example enfunc
@variables ρ_w d

T_front = (1//2)*ρ*g*H_e_c^2 - (1//2)*ρ_w*g*d^2   # traction at east face
E_front = -u_0_0 * T_front * dy                    # boundary work (dy = face length)

rhs_addition = expand_derivatives(-Differential(u_0_0)(E_front))
println("Calving-front traction added to b[row]:")
println("  ", rhs_addition)
```

### Dirichlet BC (penalty method)

At a u-point constrained to ``u_{\mathrm{BC}}``, add the penalty energy

```math
E_{\mathrm{pen}} \;=\; \tfrac{1}{2}\kappa\,
  \bigl(u - u_{\mathrm{BC}}\bigr)^2\; dx\,dy, \qquad
  \kappa \;\gg\; \max_i |A_{ii}|.
```

Its gradient ``\partial E_{\mathrm{pen}} / \partial u = \kappa(u - u_{\mathrm{BC}})``
contributes ``+\kappa`` to the diagonal entry of ``A`` and ``+\kappa\,u_{\mathrm{BC}}``
to the corresponding entry of ``b``.  As ``\kappa \to \infty`` the solution
converges to ``u = u_{\mathrm{BC}}``.

```@example enfunc
@variables u_BC κ

E_pen = (1//2)*κ*(u_0_0 - u_BC)^2*dx*dy
∂E_pen = expand_derivatives(Differential(u_0_0)(E_pen))
println("∂E_pen / (dx·dy·∂u₀) = ", expand(∂E_pen / (dx*dy)))
println("→ adds κ to A[row,row] and κ·u_BC to b[row]")
```

## Stiffness matrix from the Hessian

For a quadratic ``E_{\mathrm{disc}}``, the stiffness matrix is exactly the
negated, area-normalised Hessian:

```math
A_{kl} \;=\;
  -\frac{1}{dx\,dy}\;
  \frac{\partial^2 E_{\mathrm{disc}}}{\partial x_k \,\partial x_l}.
```

Because `Symbolics.jacobian` is used twice (once by the energy formulation to
get the residual gradient, once more to get the stiffness), both the
energy-based and the existing strong-form pipelines produce **identical**
``A``.  The Hessian route also makes it straightforward to assemble only the
upper triangle and exploit the guaranteed symmetry.

See [`examples/energy_functional_demo.jl`](https://github.com/fesmc/IceSheetStencils.jl/blob/main/examples/energy_functional_demo.jl)
for a self-contained script that carries out all of the above steps,
including the Hessian extraction and the symbolic verification.
