# Theory audit for the full-cylinder WNL code

This note records the mathematical convention implemented by the current
`weakly_nonlinear` module. It supersedes the parts of the downloaded v17.1
manuscript that describe every retained mode as a complex neutral mode and
every calculation as a common-neutral-point Landau expansion.

## 1. Governing residual and sign convention

After mapping the moving layers to fixed reference domains, the code uses the
descriptor residual

```text
B0*q_t - L(t;a)*q - C(q,q) - D(q,q,q) + O(||q||^4) = 0,
```

where

```text
C(x,y)   = -(1/2) R_qq(0)[x,y],
D(x,y,z) = -(1/6) R_qqq(0)[x,y,z].
```

Thus the evolution form is `B0*q_t=L*q+C(q,q)+D(q,q,q)+...`. The signs in
the coefficient routines are consistent with this convention.

`B0` is singular by design. Momentum and interface-kinematic variables carry
time derivatives; pressure, incompressibility, traction, boundary, volume,
and gauge rows are algebraic. Direct and adjoint vectors are normalized by

```text
y_j^H * Bslow * phi_j = 1.
```

## 2. Floquet modes and amplitude coordinates

The shifted Floquet block used by the code is

```text
A_{m,s}(Lambda) =
  [Lambda + i*omega*(n+s)]*B0*delta_nn' - L_{n-n'}.
```

At the operating point, a retained mode satisfies
`A_{m_j,s_j}(lambda_j)*phi_j=0`. Product fields use the sum of the input
exponents, not an unshifted neutral block.

There are two different physical coordinate types.

1. A distinct azimuthal/conjugate pair is represented by a complex internal
   coordinate `a_j`:

   ```text
   q_1 = a_j*phi_j + conj(a_j*phi_j).
   ```

2. A self-conjugate axisymmetric harmonic or subharmonic mode with a real
   Floquet exponent is represented by one real signed coordinate `x_j`:

   ```text
   q_1 = x_j*phi_j,   x_j in R.
   ```

The direct interface carrier is normalized to unit complex peak. Therefore
the user-facing peak displacement `A_j=(peak zeta_j)/h` is related to the
internal coordinate by

```text
A_j = S_j*a_j,   S_j=2 for a distinct complex pair,
                     1 for a self-conjugate real mode.
```

Consequently a cubic coefficient expressed in user-facing peak-amplitude
units is scaled column by column:

```text
g_peak(j,k) = g_internal(j,k)/S_k^2.
```

For a real mode, an initial phase is not a continuous degree of freedom. A
phase of zero gives a positive signed amplitude and a phase of pi gives a
negative one. Other phases are rejected.

## 3. Two related, but different, reductions

### Common-neutral multiple scales

Near a simple neutral instability, the distinguished balance is

```text
q = epsilon*q1 + epsilon^2*q2 + epsilon^3*q3 + ...,
a-ac = epsilon^2*Deltaa,
T = epsilon^2*t.
```

This balance makes detuning, slow evolution, and cubic feedback enter at the
same order. In physical peak-amplitude coordinates it gives

```text
dA_j/dt = (a-ac)*mu_j*A_j + cubic terms + O(A^5).
```

This is the controlled Landau interpretation. Two retained modes require a
common neutral point (a codimension-two condition), unless one mode is
strictly slaved rather than retained.

The exact linear Floquet problem does not require this scaling. It can remain
accurate far from `ac` at small disturbance amplitude because it discards
nonlinear powers of the disturbance, not powers of `a-ac`.

### Operating-point local reduction

With `reference='analysisAmplitude'`, the code instead evaluates
`lambda_j`, the direct/adjoint modes, the nonlinear actions, and all shifted
homological equations at the requested acceleration. Its reduced equation is
an adjoint-projected, slaved-field cubic model in modal coordinates:

```text
dA_j/dt = lambda_j*A_j
        + A_j*sum_k g_peak(j,k)*|A_k|^2
        + conj(A_j)*sum_k h_peak(j,k)*A_k^2
        + O(|A|^4).
```

The last term is present only when azimuthal and Floquet selection rules allow
it. This operating-point equation does not require `lambda_j` to be small in
order to define a finite-time cubic correction. Unlike the neutral-point
Landau coefficient, however, an off-neutral cubic coefficient is not
coordinate invariant: nonresonant terms can be redistributed between the
modal ODE and higher-order spatial embedding by a near-identity coordinate
change. The code fixes one adjoint-projection/slaved-field convention, but it
does not compute the third-order embedding. Treating its indefinite
integration as a universal saturation law additionally requires slow rates,
spectral isolation, coefficient convergence, and external validation.

For a controlled early-time comparison, the recommended option is
`transientModel='smallAmplitudeCorrection'`. It evaluates the cubic forcing
on the exact linear trajectory. With

```text
I(z,t) = (exp(z*t)-1)/z,   I(0,t)=t,
A_L,j(t) = A_j(0)*exp(lambda_j*t),
```

the implemented first correction is

```text
delta A_j^g = A_L,j * sum_k g_jk*|A_k(0)|^2
                         * I(2*Re(lambda_k),t),

delta A_j^h = sum_k h_jk*conj(A_j(0))*A_k(0)^2*exp(lambda_j*t)
                         * I(conj(lambda_j)+2*lambda_k-lambda_j,t).
```

This directly measures departure from exact linear growth. The code truncates
the returned interval when the requested amplitude, relative-correction, or
nonlinear-to-linear rate limit is reached.

## 4. Correct second- and third-order combinatorics

Let `A(Lambda)` denote the shifted product block and let `N_j=1` be the
descriptor normalization.

### Distinct complex pair

For one complex mode,

```text
A(2*lambda_j) q_jj = C(phi_j,phi_j),
A(lambda_j+conj(lambda_j)) q_jbarj = 2*C(phi_j,bar(phi_j)),

g_jj = y_j^H [2*C(phi_j,q_jbarj)
              +2*C(bar(phi_j),q_jj)
              +3*D(phi_j,phi_j,bar(phi_j))].
```

For a complex source mode `k` in target equation `j`,

```text
g_jk = y_j^H [2*C(phi_j,q_kbark)
              +2*C(phi_k,q_jbark)
              +2*C(bar(phi_k),q_jk)
              +6*D(phi_j,phi_k,bar(phi_k))].
```

When both coordinates are complex and

```text
2*m_k-m_j = m_j,
wrap(2*s_k-s_j) = s_j,
```

the symmetry-complete equation also contains
`h_jk*conj(A_j)*A_k^2`, with

```text
h_jk = y_j^H [2*C(phi_k,q_barj,k)
              +2*C(bar(phi_j),q_kk)
              +3*D(bar(phi_j),phi_k,phi_k)].
```

This term is required, for example, for complex harmonic/subharmonic modes
with the same azimuthal wavenumber.

### Self-conjugate real mode

For a real signed mode there is no independent conjugate mean field:

```text
A(2*lambda_j) q_jj = C(phi_j,phi_j),
g_jj = y_j^H [2*C(phi_j,q_jj)+D(phi_j,phi_j,phi_j)].
```

For a real source `k`,

```text
A(2*lambda_k) q_kk = C(phi_k,phi_k),
A(lambda_j+lambda_k) q_jk = 2*C(phi_j,phi_k),

g_jk = y_j^H [2*C(phi_j,q_kk)
              +2*C(phi_k,q_jk)
              +3*D(phi_j,phi_k,phi_k)].
```

If the source is real but the target is a distinct nonaxisymmetric pair,
`q_jk` and the three target-frequency feedback terms remain complex fields.
Only `q_kk` is self-conjugate. Applying a real-field projection to `q_jk`
would impose a false `m=0` symmetry and is therefore forbidden.

Applying the complex-pair formulas to a real axisymmetric mode overcounts a
scalar self coefficient by a factor of three in the simplest polynomial test
and overcounts a real-source cross coefficient by a factor of two. The code
now dispatches these cases separately.

### Axisymmetric pressure quotient

The equal-order primitive discretization can contain pressure-only right null
vectors `Z` satisfying `A*Z approximately 0` and `Bslow*Z=0`. Their matching
left nullspace represents redundant discrete compatibility directions, not a
physical velocity or interface degree of freedom. For an `m=0` forced field,
the optional pressure-quotient gate first certifies this algebraic right
nullspace and then decomposes

```text
r = A*q-f = r_physical + r_compatibility.
```

The field is accepted only if `norm(r_physical)/norm(f)` passes the ordinary
forced-field tolerance and `norm(r_compatibility)/norm(f)` passes a separate,
looser ceiling. The unprojected residual is always retained. This procedure
is never applied to nonaxisymmetric blocks or to a nullspace that fails the
zero-mass/algebraic certification.

The compact equation is not complete whenever symmetry admits cubic
monomials outside the implemented `g/h` families. Two independent radial
branches with identical `(m,s)` permit transfers such as
`A_k*|A_j|^2`. A 1:3 azimuthal pair can permit `A_1^3` and
`conj(A_1)^2*A_2`. The code enumerates the allowed cubic monomials and rejects
such a mode set before an expensive recovery. A full cubic coefficient tensor
is needed to support it.

## 5. Physical assumptions

The cylinder implementation assumes:

- two immiscible, incompressible, Newtonian fluids with constant density,
  viscosity, and interfacial tension;
- a flat, quiescent base interface with no imposed circulation or mean flow;
- equal nondimensional layer depths, `-1 <= z <= 0` and `0 <= z <= 1`;
- sinusoidal vertical acceleration represented through the hydrostatic
  normal-traction term;
- rigid no-slip horizontal walls;
- a selectable impermeable sidewall model. The preferred stress-free model is
  `u_r=0`, `d_r w=0`, and `R*d_r(u_theta)-u_theta=0`; the retained
  `legacyFreeSlide` option instead imposes `d_r(r*u_theta)=0`. Neither is a
  no-slip sidewall;
- either an ideal free contact line, `d_r zeta=0`, or an ideal pinned contact
  line. The bundled reduced Bessel mode is an exact separated radial carrier
  only for the compatible legacy/free-wall idealization. With the true
  stress-free wall and `m>0`, it is used only as a tracking seed and the full
  reconstructed interface harmonics are allowed to acquire radial corrections
  `c_n(r)`. For the pinned case the nonlinear sidewall rows retain the ALE
  radial-derivative correction generated by a nonzero contact-angle
  perturbation;
- no dynamic contact angle, contact-line hysteresis, wetting film, surfactant,
  Marangoni stress, compressibility, thermal effect, or turbulence;
- a single-valued interface graph and a nonsingular ALE map;
- fixed volume for the axisymmetric interface component.

The ideal stress-free or free-slide sidewall and ideal contact-line models are
likely the largest physical simplifications for comparison with a laboratory
cylinder. A physical no-slip wall would additionally require resolution of
the sidewall viscous boundary layer.

For each temporal harmonic, the implementation reports the cylindrical
weighted projection

```text
zeta_n(r) = b_n J_m(beta r) + c_n(r),
<J_m,c_n>_(r dr) = 0.
```

This is a diagnostic decomposition of the full eigenmode, not an additional
closure assumption. The WNL coefficients use the full mode regardless of the
plotting policy. The optional `auto` policy omits `c_n` from a reconstructed
surface only when its combined time-radius relative L2 norm is below the user
threshold; this omission should be checked under radial and vertical grid
refinement.

The reconstructed second-order fields are the particular slaved responses
of the modal parameterization and are therefore present at the initial plotted
time. This assumes the initial full fluid state lies on that local invariant
manifold. A laboratory initial condition containing only the primary
interface displacement also excites homogeneous velocity/pressure transients;
capturing their startup requires a time-domain Volterra/DNS calculation, not
the instantaneous slaved-field reconstruction used here.

## 6. What must be checked before using a coefficient quantitatively

Numerical convergence and asymptotic validity are separate requirements.

1. Direct and adjoint residuals must pass their reported physical gates, and
   `y^H*Bslow*phi` must remain one.
2. The full and reduced operating-point exponents must agree under refinement.
3. Every second-order forced field must pass the unequilibrated
   forcing-relative residual gate. An exploratory field is not a converged
   field.
4. Quadratic resonances and near-resonant large slaved responses must be
   excluded or retained in an enlarged amplitude system.
5. `N`, `Nr`, both vertical grids, `Ntheta`, and directional-difference steps
   must be varied until `lambda`, `g`, and `h` stabilize.
6. The reconstructed interface amplitude and slope must stay small, and the
   ALE Jacobian must stay positive.
7. For the finite-time correction, both
   `|delta A_j|/|A_L,j|` and `|N_j(A_L)|/|lambda_j A_L,j|` must stay below the
   chosen perturbative limits.
8. For a long cubic-envelope trajectory away from onset, compare against DNS
   or experiment; small residuals alone do not establish physical accuracy.

The reconstruction includes primary fields and available second-order slaved
fields. It does not include the third-order spatial embedding, so its modal
ODE is cubic while its displayed interface shape is accurate only through
second order.

## 7. Audit changes implemented

- separated real self-conjugate and complex-pair coefficient formulas;
- made peak-displacement normalization explicit and mode dependent;
- rejected nonphysical continuous phases for real modes, including saved-data
  postprocessing overrides;
- added phase-sensitive cubic terms for compatible complex H/SH pairs;
- rejected every retained mode set whose symmetries require cubic monomials
  outside the implemented `g/h` families, including same-`(m,s)` and 1:3
  azimuthal pairs;
- corrected physical reconstruction multiplicities, including mixed
  real/complex mode pairs;
- retained the nonlinear ALE sidewall corrections required by a pinned
  contact line;
- replaced trapezoidal integration on Chebyshev nodes by Clenshaw-Curtis
  quadrature for the axisymmetric volume constraint;
- made old caches fail explicitly when their coefficient convention omits one
  of these required terms. Their direct/adjoint mode vectors may still be
  reused when recomputing coefficients.
