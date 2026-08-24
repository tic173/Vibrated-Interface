# Vibrated-Interface
Multiphase interfacial dynamics in vibrated domains

## Weakly nonlinear analysis

The `weakly_nonlinear` folder contains the MATLAB coefficient engine for
Floquet modes, including adjoints, quadratic-resonance tests, forced
second-order fields, self/cross Landau coefficients, and the full-state
`cylinder_wnl_operators.m` ALE discretization. Start with
`weakly_nonlinear/README.md` and run
`weakly_nonlinear/tests/run_wnl_tests.m`.

The existing reduced interface solver supplies an independent threshold and
mode-validation target. Converge the full-state linear recovery before using
the physical nonlinear coefficients quantitatively.

The supported editable driver is
`weakly_nonlinear/examples/vi_wnl_user_run_full_cylinder.m`. It accepts one or
two Floquet modes, a user-defined magnitude and phase for each initial complex
amplitude, and any requested vibration acceleration. With
`input.weaklyNonlinear.reference='analysisAmplitude'`, both modes and all
nonlinear coefficients are evaluated at that physical operating point; the
modes do not need a common critical acceleration. The full implementation
uses Bessel-enriched radial differentiation, matched multi-domain vertical
Chebyshev grids, full primitive-variable eigenpair refinement, and separate
linear, coefficient, and forced-field residual gates.

The editable driver defaults to `input.execution.profile='balanced'`. Use
`development` only for a quick feasibility/recovery pass and `final` for a
resolution study and reportable coefficients. Full-eigenpair corrections
anchor their phase and scale to the interface displacement and preserve the
physical `Bslow*phi` descriptor branch, so an algebraic pressure/gauge null
vector cannot masquerade as a converged interface mode.
A descriptor-aware residual normalization excludes zero-mass pressure/gauge
variables from the state norm while retaining every primitive equation in the
residual numerator. The LSQR correction fallback also works in physical
coordinates through a regularized objective assembled before column scaling,
preventing that scaling from favoring a pressure-inflated candidate. Poorly
solved bordered corrections are rejected before line search.
A promising GMRES update whose descriptor direction is correct but whose norm
is dominated by zero-mass pressure/gauge variables receives a descriptor-DAE
algebraic completion. The velocity/interface variables are held fixed while
only algebraic variables are recomputed in a regularized minimum-norm solve;
the completed state must then pass the original physical residual and branch
gates.
The same pressure-safe DAE completion now treats the nonzero right-hand sides
of quadratic mean/combination/second-harmonic fields. Forced corrections are
physically regularized, while acceptance still uses the original
forcing-relative full-operator residual. Two-mode cross branches fail fast
after an invalid required field, avoiding downstream solves that cannot
restore that coefficient.
The cylinder forced-field solver now exploits that vibration couples adjacent
Floquet harmonics only through interface-displacement columns. It factors the
primitive spatial blocks, solves a small temporal interface Schur system, and
reconstructs the full state before checking the original equation residual.
The column-equilibrated adaptive-rank minimum-norm calculation remains the
fallback; the older temporal-block GMRES path remains optional.
The cylinder nonlinear action also removes the exact flat-interface linear
remainder before directional differentiation and checks a factor-of-two step
sequence with Richardson extrapolation. This prevents linear cancellation
error from appearing as a spurious incompatible quadratic forcing.
A near-converged full eigenpair uses temporal-block preconditioned corrections
whose GMRES candidate is physically tested before any LSQR fallback. It can
continue algebraically without repeating its expensive interface lift;
coefficient evaluation still requires both strict residual convergence and
agreement between the reduced and full operating-point exponents.
Coefficient runs can also reuse coefficient-ready direct/adjoint modes from a
saved recovery MAT-file. Reuse requires an exact operating-point/grid/branch
signature match and rechecks every current full-operator residual and
reduced/full exponent gate. By default, a coefficient run stops immediately
on a cache miss instead of silently repeating an hours-long mode recovery.
For two modes, numerical coefficient validity and slow-envelope validity are
reported separately. An arbitrary-acceleration operating-point reduction may
be computed for diagnostics, but the driver withholds a quantitative WNL
transient when either retained exponent is not slow relative to the forcing.
The conjugate cylinder adjoint also applies the required swap of the two
axis-regularity equation rows under `m -> -m`; direct states and algebraic
adjoints therefore use their mathematically distinct conjugation maps.
