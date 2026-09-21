# Linear Floquet stability and interface dynamics

This module computes viscous two-fluid Floquet growth rates, multipliers, and
small-amplitude interface motion from physical inputs. It supports 2D Cartesian
(x,z), 3D Cartesian (x,y,z), and 3D cylindrical (r,theta,z) geometries using the
repository's existing separable interface equation.

Requires MATLAB and Optimization Toolbox (`fsolve`). The plotting driver uses
`tiledlayout` (MATLAB R2019b or later); verification was run in R2023b.

## Start here

From the repository root:

```matlab
run('linear_floquet/examples/run_linear_floquet.m')
```

Edit the **User inputs** section of that script. Select the geometry in the
`vi_linear_defaults(...)` call before editing its mode fields. The driver defaults to **cartesian2d** and scans 81 effective wavenumbers
`kStarRange=linspace(0.2,20,81)`, where `kh=k_eff*lowerDepth`. It computes
Floquet growth rates for every sample, then reconstructs interface dynamics
**only for the sample with the largest growth rate**. If all rates are negative,
this is the least damped mode. Refine the range and spacing near the maximum;
the selected point is a sampled maximum, not a continuous optimization.

The Cartesian scan uses unbounded Fourier modes; `Lx` and `Ly` set the plotted
window. For a periodic box, select `horizontalBoundary='periodic'` and supply
compatible discrete Fourier modes programmatically. Horizontal wall conditions
are unchanged. For cylindrical geometry, the driver instead scans the requested
`azimuthalOrders` and `radialIndices` at the fixed radius, respecting the discrete
Bessel spectrum. In 3D Cartesian geometry, `directionRad` chooses wave direction.

The driver checks temporal resolution, plots growth versus `kh` with the selected
point marked, plots that mode's dynamics, and saves a new folder for each run.
`growthRates` and `result.sweep.table` contain the entire scan; `result.sweep.modes`
preserves its eigenvalues, harmonics and convergence diagnostics.
`result.modes`, `result.config.modes`, and the interface arrays contain only the
selected mode. `growth_rates.csv` contains all samples; `modal_dynamics.csv`
contains only the selected mode. The MAT file preserves both.

Use `vi_linear_growth(cfg)` for growth rates alone, `vi_linear_most_unstable(cfg)`
for a scan followed by the selected dynamics, or `vi_linear_floquet(cfg)` for
superposed dynamics of every supplied mode. `vi_linear_dynamics(growthResult)`
reconstructs a previously computed solution without repeating root searches.

For programmatic use:

```matlab
addpath('linear_floquet')
cfg = vi_linear_defaults('cartesian2d');
cfg.geometry.Lx = 0.08;
cfg.modes.kx = 2*pi/cfg.geometry.Lx;
cfg.forcing.kind = 'acceleration_g';
cfg.forcing.amplitude = 0.8;
cfg.forcing.frequencyHz = 20;
cfg.forcing.phaseRad = pi/4;
cfg.sampling.timeSeconds = linspace(0,10/cfg.forcing.frequencyHz,401);
result = vi_linear_floquet(cfg);
```

The defaults are water/air with two 22 mm layers, a 35 mm cylinder radius or
70 mm Cartesian periods, 15 Hz forcing, acceleration amplitude 0.3 g, and an
initial modal displacement scale of 10 micrometres. Changing the forcing
frequency does not automatically change `sampling.timeSeconds`; set both if
you want a fixed number of forcing cycles.

## Required physical and modal inputs

Start with `vi_linear_defaults` to obtain a complete configuration. Missing
fields, misspelled fields, unsupported boundary conditions, zero wavenumbers,
and inconsistent periodic wavenumbers are rejected before solving.

| Input | Meaning / units |
|---|---|
| `fluids.rhoLower`, `rhoUpper` | Densities, kg/m^3 |
| `fluids.muLower`, `muUpper` | **Dynamic** viscosities, Pa s (mu = rho*nu) |
| `fluids.surfaceTension` | Interfacial tension, N/m |
| `geometry.lowerDepth`, `upperDepth` | Positive distances from the interface to the two horizontal walls, m; unequal depths supported |
| `geometry.radius` | Cylinder radius, m |
| `geometry.Lx`, `Ly` | Cartesian periods, or plotting-window lengths for unbounded modes, m |
| `gravity.magnitude` | Positive gravity/reference acceleration, m/s^2 |
| `gravity.sign` | +1 or -1, multiplying only the mean gravitational term |
| `forcing.frequencyHz` | Positive frequency in Hz; omega = 2*pi*f |
| `forcing.amplitude`, `kind` | Nonnegative amplitude and its units, described below |
| `forcing.phaseRad` | Initial phase of the acceleration cosine, radians |
| `modes` | Nonempty struct array of requested spatial modes |
| `modes(i).amplitude_m`, `phaseRad` | Magnitude and phase of the initial **complex modal interface displacement** |
| `sampling.timeSeconds` | Increasing physical times, s |

The other important inputs are the boundary condition and the spatial mode
indices/wavenumbers. Geometry and material properties alone do not select a
unique disturbance. Multiple requested modes are superposed linearly.

Lower and upper label the two sides of the reference interface, at z<0 and z>0.
Either density ordering is allowed, but this density-jump formulation rejects
nearly equal densities (`abs(At)<1e-8`). Positive mean gravity stabilizes a denser
lower fluid. Reversing `gravity.sign` reverses that gravitational contribution;
RT ordering corresponds to `(rhoLower-rhoUpper)*gravity.sign < 0`.

### Vibration amplitude and phase

The acceleration convention, shared with the existing reduced/full operators,
is

```text
g_eff(t) = g*gravity.sign + a*cos(2*pi*f*t + forcing.phaseRad).
```

`forcing.kind` is one of:

- `acceleration_g`: amplitude is a/g.
- `acceleration_m_s2`: amplitude is a in m/s^2.
- `displacement_m`: amplitude is A in m, with a = A*(2*pi*f)^2.

**The phase always describes acceleration**, including when the amplitude is
entered as displacement. A base displacement `A*cos(omega*t+phi_z)` has
acceleration phase `phi_z+pi`. Changing the forcing phase changes the periodic
waveform/time origin, not the exact Floquet growth rate.

`modes(i).phaseRad` is a separate disturbance coefficient phase. It does not
change the applied vibration. With the normalization below, at t=0 the initial
interface is `real(sum(A_i*Phi_i))`, where
`A_i=amplitude_m*exp(1i*phaseRad)`. Thus, for a real axisymmetric shape and a phase
of pi/2, the initial real displacement can vanish even though the complex
modal amplitude is nonzero.

### Geometry and boundary conditions

All cases retain no slip and no penetration at the two horizontal walls,
velocity continuity and tangential stress continuity at the interface, and
capillary/normal stress balance. There is no base shear flow.

| Geometry selector | Mode fields | Spatial factor before taking the real part |
|---|---|---|
| `cartesian2d` | `kx` in rad/m | exp(i*kx*x) |
| `cartesian3d` | `kx`, `ky` in rad/m | exp(i*(kx*x+ky*y)); k=hypot(kx,ky) |
| `cylindrical3d` | integer `m>=0`, positive integer `radialIndex` | J_m(k*r)*exp(i*m*theta), with J'_m(k*R)=0 |

Cartesian `geometry.horizontalBoundary='periodic'` is the default. Then kx*Lx
and ky*Ly must be integer multiples of 2*pi. For continuous Fourier waves,
select `'unbounded'`; Lx and Ly then define only the displayed window. Lateral
no-slip walls or Cartesian pinned edges are not imposed by Fourier modes.

The cylinder uses `'bessel-neumann'`, the repository's free-contact-line,
separable radial model. Its wavenumbers are the **positive** derivative roots
`k=j'_(m,j)/R`; the spatially uniform volume-changing mode is excluded. This is
not the full no-slip sidewall problem. Pinned contact lines are not supported
by this driver; the existing pinned solvers remain separate legacy tools.

Cartesian complex modes represent traveling-wave components. Combine opposite
wavevectors with matching temporal amplitudes to form Cartesian standing
waves. Cylindrical modes are normalized so their continuous radial maximum is
one; Cartesian spatial factors already have unit modulus.

Examples of mode selection:

```matlab
% One oblique Fourier mode in a doubly periodic Cartesian domain
cfg = vi_linear_defaults('cartesian3d');
cfg.modes.kx = 2*pi/cfg.geometry.Lx;
cfg.modes.ky = 4*pi/cfg.geometry.Ly;

% Two cylindrical modes
cfg = vi_linear_defaults('cylindrical3d');
cfg.modes(1).m = 0; cfg.modes(1).radialIndex = 1;
cfg.modes(2) = cfg.modes(1);
cfg.modes(2).label = 'm=2, j=1';
cfg.modes(2).m = 2;
cfg.modes(2).amplitude_m = 5e-6;
cfg.modes(2).phaseRad = pi/3;
```

## What is solved and reconstructed

Let h be `lowerDepth`, tc=sqrt(h/g), kStar=k*h, and omegaStar=omega*tc. The code
calculates

```text
At  = (rhoLower-rhoUpper)/(rhoLower+rhoUpper)
eta = muUpper/muLower
C   = (muLower/rhoLower)/sqrt(g*h^3)
Bd  = rhoLower*g*h^2/surfaceTension
Ac  = a/g.
```

It uses `vi_reduced_cylinder_coefficients` and
`vi_floquet_acceleration_matrix` to solve

```text
M(sStar)*z = [diag(An(sStar)) - Ac*B(phase)]*z = 0.
```

The same vertical equation applies to each geometry's horizontal Laplacian
mode; no Mathieu approximation replaces the viscous fluid calculation.
`vi_reduced_cylinder_coefficients` now accepts optional dimensionless layer
depths, defaulting to [1 1], so existing equal-depth callers remain compatible.
Its wall-shifted exponential basis avoids exponentially growing basis values.
An exactly static harmonic is evaluated analytically.

A multistart search solves for a **complex** exponent with indices -N:N. The
returned exponent already contains the quasi-frequency: no extra omega/2 is
added for a subharmonic response. Harmonic and subharmonic responses occur at
quasi-frequencies equivalent to 0 and omega/2; away from these, oscillatory
complex branches can occur. Integer-frequency aliases are re-solved in the
principal frequency interval with their matching harmonic vectors.

The reconstruction is

```text
p_i(t) = sum_n z_(i,n)*exp(i*n*omega*t), normalized to p_i(0)=1
zeta(x,t) = real(sum_i A_i*Phi_i(x)*exp(sStar_i*t/tc)*p_i(t)).
```

The physical growth rate is real(sStar)/tc in 1/s, and the multiplier for one
forcing period is exp(2*pi*sStar/omegaStar). A positive rate means growth for
that retained mode. The full interface is real; complex modal traces are also
saved so normalization and phase can be inspected. If p(0) is too close to a
node for this normalization, the code asks for a different forcing phase.

This is an **eigenmode-seeded linear evolution**, not the solution for an
arbitrary prescribed initial interface and fluid velocity. Each retained
Floquet eigenmode determines its compatible initial velocity. Specifying only
an arbitrary interface shape cannot determine a unique viscous initial-value
solution. The calculation has no nonlinear saturation and applies while
perturbations remain small.

## Numerical checks and scope of stability conclusions

- `numerics.harmonics`: N for the temporal indices -N:N.
- `checkConvergence=true`: repeat at N+`harmonicIncrement` (default 4), retaining
  the higher-resolution result only if the dominant exponent agrees.
- `initialGuesses`: optional complex dimensionless exponent guesses. Empty
  uses a grid of real/quasi-frequency seeds plus finite-depth inviscid seeds.
- `rootTolerance`: scaled matrix/null-vector residual tolerance.
- `edgeTolerance`: maximum energy fraction in the outer temporal coefficients.
- `convergenceTolerance`: relative exponent-change tolerance.

The root finder accepts only converged nonlinear solves with a small matrix
residual and small harmonic-edge energy. It reports the largest growth rate
**among accepted roots for each requested spatial mode**. A finite multistart
search is not a proof that every viscous Floquet branch was found. Broaden the
initial guesses and increase N when surveying new regimes. Stable requested
modes do not establish stability against unrequested spatial modes. The
returned `rootSearch` records accepted candidates and rejected trial messages.
A reported truncation agreement checks the exponent, not every point of the
reconstructed waveform; refine N and the output samples for quantitative
waveform studies.

## Results, output files, and folder organization

```text
linear_floquet/
  examples/run_linear_floquet.m   editable driver
  vi_linear_defaults.m           complete input configuration
  vi_linear_floquet.m            public analysis function
  vi_linear_plot.m               public plotting function
  private/                       validation, root search, reconstruction, saving
  tests/run_linear_tests.m       regression/physics checks
  results/run_*/                 generated outputs (git-ignored)
```

Shared root-level solvers were retained in place to avoid breaking the existing
validation scripts and weakly nonlinear module. New linear workflow code is
contained in this folder.

Important return fields:

- `result.modes(i).growthRatePerSecond`, `exponentStar`, `multiplier`.
- `harmonicIndices`, `harmonicCoefficients`, `harmonicConvergenceError`.
- `result.modalDisplacement_m`: complex matrix [number of modes, Nt].
- `result.interface_m`: real displacement in m. Array sizes are [1,Nx,Nt],
  [Ny,Nx,Nt], or [Ntheta,Nr,Nt] for the three geometries respectively.
- `result.grid.X`, `.Y`: physical Cartesian coordinates for plotting, including
  the cylindrical surface; cylindrical `.r` and `.theta` are also retained.
- `result.timeSeconds`, `config`, `parameters`: reproducible inputs and units.

Set `output.reconstruct=false` for a growth-rate/modal-trace sweep without
allocating the full spatial-time array. Set `output.plot=false` for batch work.
When `output.save=true`, a unique folder contains:

- `growth_rates.csv`: labels, k, kh, growth rates, quasi-frequencies, multipliers.
- `modal_dynamics.csv`: time and real/imaginary modal displacements.
- `linear_floquet.mat`: the complete result, configuration, diagnostics, and
  reconstructed interface if requested.

To replot a saved result:

```matlab
saved = load('path/to/linear_floquet.mat');
vi_linear_plot(saved.result);
```

## Verification

```matlab
addpath('linear_floquet/tests')
run_linear_tests
```

The checks cover equivalent wavenumbers in all three geometries, forcing-phase
invariance and coefficient rotation, amplitude-unit conversion, one-period
multipliers, compatibility with the existing cylindrical growth solver, the
unforced low-viscosity finite-depth dispersion relation at unequal depths,
RT growth, multiple-mode reconstruction, saved outputs, and rejected inputs.

The Floquet formulation follows the approach described by
[Kumar and Tuckerman (1994), *Parametric instability of the interface between
two fluids*](https://doi.org/10.1017/S0022112094003812). The repository's actual
shared viscous coefficient implementation is the numerical model used here.
