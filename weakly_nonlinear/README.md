# Weakly nonlinear Floquet coefficient module

This folder implements the numerical reduction described in
`Weakly_Nonlinear_RT_Faraday_Cylinder.tex`.

## Scope

The following parts are implemented:

- Fourier–Floquet descriptor operator assembly;
- direct and algebraic-adjoint null vectors;
- normalization with the physical slow-time matrix `B0`;
- correct conjugation of harmonic, subharmonic, and general Floquet fields;
- azimuthal and temporal selection rules;
- bordered second-order solves;
- detection and projection of quadratic resonances;
- self-interaction coefficients;
- cross-interaction coefficients between different azimuthal modes;
- shared two-mode mean/sum/difference fields without duplicate forced solves;
- shifted Floquet blocks for operating-point modal reductions;
- nonresonant coupled Landau equations;
- a primitive-variable Bessel-enriched radial, multi-domain
  Chebyshev--Fourier cylinder discretization;
- full ALE bulk/interface residuals and field-level directional actions;
- free and pinned contact-line rows at every perturbation order;
- cancellation-safe directional actions of the exact nonlinear residual
  remainder, with adaptive step-doubling and Richardson checks.

The existing repository eliminates the bulk fields and returns a reduced
linear equation for `zeta_all`. That is sufficient for linear thresholds but
not for a physical Landau coefficient. `cylinder_wnl_operators.m` now adds an
independent full-state discretization returning `B0`, `Lhat`, `C`, `D`, and
`P`. It must recover the reduced solver's neutral point under grid refinement
before its nonlinear coefficients are used quantitatively.

No empirical or guessed nonlinear terms are inserted by this module.

## State ordering

For one temporal harmonic, the cylinder factory uses

```text
q_n = [u_r^d; u_theta^d; w^d; p^d;
       u_r^l; u_theta^l; w^l; p^l; zeta]
```

including all constraint and boundary unknowns needed by the chosen
discretization. A Floquet field stores the coefficients as

```matlab
field.coeff(:,j)  % spatial state at temporal index field.spec.n(j)
```

and the corresponding matrix vector is `field.coeff(:)`.

The block specification contains:

- `m`: azimuthal wavenumber;
- `s`: Floquet quasifrequency in `(-1/2,1/2]`;
- `lambda`: continuous modal exponent after removing the periodic carrier;
- `n`: retained integer temporal indices;
- `ndof`: spatial state size for one temporal harmonic.

For `s=1/2`, conjugation maps `n` to `-n-1`. The helper functions perform
this shift automatically.

## Quick verification

In MATLAB:

```matlab
cd Vibrated-Interface
addpath('weakly_nonlinear')
addpath('weakly_nonlinear/examples')
addpath('weakly_nonlinear/tests')
run_wnl_tests
```

## Editable cylinder input file

Start from `examples/vi_wnl_user_run_full_cylinder.m` when
changing physical parameters, the gravity sign, forcing, selected `(m,l,s)`
mode, Floquet cutoff, or spatial resolution. The script calculates
`omegaStar`, `R0`, `C`, `Bd`, `At`, and `eta`, outputs the linear growth rate
and Floquet multiplier, plots the reconstructed interface transient, and
compares it with the Landau prediction from the included WNL operators.
From the repository root, run

```matlab
run(['weakly_nonlinear/examples/', ...
    'vi_wnl_user_run_full_cylinder.m'])
```

Choose a transparent runtime profile near the top of that file:

```matlab
input.execution.profile = 'balanced';    % default intermediate verification
% input.execution.profile = 'development'; % fast feasibility check
% input.execution.profile = 'final';     % use entered resolution unchanged
```

The development and balanced profiles cap `N`, `Nr`, automatic vertical
element sizes, nonlinear temporal oversampling, and the most expensive Krylov
iteration counts; they are not coefficient-convergence evidence. Development
is deliberately a quick feasibility/recovery pass and may stop at a residual
gate. The final profile preserves the entered discretization. Forced fields
use the cylinder's low-rank temporal Schur solve by default; it eliminates
the complete spatial blocks and solves only for the interface coordinates
that couple neighboring Floquet harmonics. The column-equilibrated rank-aware
minimum-norm path remains a checked fallback. The optional temporal-block
GMRES path is disabled for the large two-mode cylinder because it was slower
than these paths. Stage timings are saved with the result and printed when
`input.execution.reportTiming=true`.

Each full-eigenpair Newton correction is first solved by temporal-block
preconditioned GMRES and immediately tested with the unequilibrated physical
operator. A short LSQR solve is invoked only when that GMRES candidate fails
to reduce the physical residual on the tracked descriptor branch. This
replaces the repeated 20,000-iteration LSQR corrections that dominated earlier
run times:

```matlab
input.options.eigenpairUseBlockGmres = true;
input.options.eigenpairGmresRestart = 30;
input.options.eigenpairGmresMaxCycles = 80;
input.options.eigenpairLsqrFallbackMaxIterations = 3000;
input.options.eigenpairCorrectionRegularization = 1.0e-12;
input.options.eigenpairCorrectionEquationTolerance = 1.0e-1;
input.options.eigenpairUseAlgebraicCompletion = true;
input.options.eigenpairAlgebraicCompletionRegularization = 1.0e-14;
input.options.eigenpairAlgebraicCompletionMaxIterations = 1000;
input.options.eigenpairCorrectionConstraintTolerance = 1.0e-10;
input.options.eigenpairMinimumDescriptorRetention = 1.0e-2;
input.options.eigenpairMinimumDescriptorOverlap = 0.10;
input.options.eigenpairMaximumSeedResidualRatio = 5.0;
```

The LSQR fallback minimizes the row-equilibrated bordered residual plus a
small Tikhonov penalty on the physical velocity--pressure correction. The
complete objective is then column-equilibrated. Because the physical norm is
inserted before this variable transformation, column scaling improves Krylov
conditioning without redefining a huge pressure correction as cheap. A
rejected GMRES correction is supplied as the LSQR initial iterate only if its
best trial retained the required descriptor norm and overlap; a
pressure-inflated correction restarts LSQR from zero.

An inexact correction is tested in two stages. First, the full bordered
correction equation must satisfy `eigenpairCorrectionEquationTolerance`;
second, its line-search trial must reduce the physical mode residual while
retaining the descriptor branch. A correction-equation residual near one is
not a Newton direction, even if it happens to lower the mode residual by a few
percent. Such a correction is rejected immediately and automatic continuation
is disabled, preventing repeated costly nonconvergent steps.

There is one important descriptor-system exception. A GMRES trial can already
have a small physical equation residual and nearly unit overlap with the
reference `Bslow` direction, but still have near-zero descriptor retention
because pressure or another zero-mass variable has grown enormously. For that
case, `wnl_complete_algebraic_state` holds all descriptor-supported variables
and prescribed interface coefficients fixed and solves only for the algebraic
variables. The completion is accepted only when the resulting state reduces
the actual full-operator residual and passes the same phase, descriptor, and
eigenvalue checks. It cannot turn an inaccurate velocity/interface trial into
an accepted mode by changing a tolerance.

The full primitive-variable descriptor is singular because pressure and
constraint variables do not carry slow-time derivatives. Consequently,
`A*phi=0` alone is not sufficient to identify the physical interfacial mode:
an iterative correction can otherwise collapse onto a pressure/gauge null
direction with `Bslow*phi=0`. For cylinder modes, the eigenvector phase and
scale are therefore anchored to the prescribed interface-displacement
coefficients, rather than to the pressure-dominated norm of the full primitive
state. Models without explicit interface metadata use the `Bslow` variables as
the anchor. Each accepted correction must also retain both the norm and
direction of the reference descriptor vector. A residual reduction that fails
either test is rejected and the last physical iterate is retained.

The reported direct residual is

```text
||A phi|| / (||A||scale ||phi on nonzero Bslow columns||),
```

and the adjoint uses the corresponding nonzero `Bslow` rows. All pressure,
gauge, wall, interface, and constraint equations remain in the numerator.
Only zero-mass algebraic variables are excluded from the denominator. Thus a
candidate cannot make its residual arbitrarily small by adding an enormous
pressure/null component. The saved `directResidualDetails` and
`leftResidualDetails` structures include `fullToPhysicalNormRatio` for this
check.

The safeguarded line search now permits twelve halvings. This matters for a
poorly conditioned descriptor pencil: a large inexact correction can leave the
physical branch at unit step even though a smaller interface-anchored step is a
valid residual-decreasing direction. Output reports the GMRES and LSQR physical
trials separately.

At a nonzero operating-point exponent, the prescribed-interface solve is only
the starting vector for this correction. It may therefore stop at
`directLiftSeedResidualTolerance=5e-6`; this does not relax the final
`modeResidualTolerance` or `coefficientModeResidualTolerance` gates.
The correction is skipped when that starting residual remains above
`eigenpairMaximumSeedResidualRatio*modeResidualTolerance`, because a costly
bordered solve is then unlikely to remain on the requested physical branch.

Preliminary runtime profiles cap direct-correction work aggressively. The
development profile performs only the short GMRES probe and no LSQR fallback;
balanced permits 15 GMRES cycles and 750 fallback iterations. The final profile
retains the entered solver limits.

By default, a profile-limited eigenpair pass that reduces but narrowly misses
the strict coefficient gate continues from its best state. Continuation neither
repeats the prescribed-interface lift nor changes `N`, `Nr`, or the vertical
grid:

```matlab
input.options.autoEscalateEigenpairRefinement = true;
input.options.eigenpairEscalationMaxSteps = 8; % total, not additional
input.options.eigenpairEscalationMaximumResidualRatio = 1.0e4;
```

Set the first option to `false` for a strictly capped timing benchmark.
Continuation is skipped when the first-pass residual is farther from the
coefficient gate than the configured ratio, avoiding a long algebraic solve
on a clearly under-resolved grid.
Refinement also stops once two improving corrections produce a linear-ready
mode whose full exponent is clearly outside the allowed reduced/full interval.
At that point more iterations on the same grid cannot repair the discretization
mismatch.
Recovery-only runs print and save two distinct states. `linearReady` uses
`modeResidualTolerance`; `coefficientReady` applies the stricter
`coefficientModeResidualTolerance` to the primary direct/adjoint pair and the
physical conjugate direct field. Do not start coefficient evaluation from a
result labeled `LINEAR-READY ONLY`.

Coefficient readiness also requires the refined full-state exponent to agree
with the reduced operating-point exponent:

```matlab
input.options.enforceReducedFullEigenvalueAgreement = true;
input.options.maximumReducedFullEigenvalueRelativeMismatch = 0.05;
input.options.maximumReducedFullEigenvalueAbsoluteMismatch = 1.0e-3;
```

The allowed absolute shift is the larger of the absolute tolerance and the
relative tolerance times the exponent magnitude. This prevents a coarse-grid
full eigenpair with a small algebraic residual but a materially different
growth rate from entering the Landau coefficient calculation.

Specify the physical vibration amplitude independently of the neutral point:

```matlab
input.forcing.analysisAmplitude = 0.42;  % actual a/g0 to analyze
input.forcing.findCriticalAcceleration = true;
input.forcing.aCritical = 0.30;          % ignored when the line above is true
input.weaklyNonlinear.reference = 'analysisAmplitude';
```

The exact linear Floquet growth rate and transient are evaluated at
`analysisAmplitude`. If automatic threshold calculation is enabled, the code
calculates the first mode's `aCritical` separately as an onset diagnostic.
With the default `reference='analysisAmplitude'`, each retained mode is solved
at the requested physical acceleration and written as

```text
q_j(t) = exp(lambda_j t) sum_n q_jn exp[i(n+s_j) omega t].
```

The full cylinder operators and the nonlinear coefficients are assembled at
that same acceleration. The coupled amplitude equation therefore uses the
exact operating-point exponents,

```text
dA_j/dt = lambda_j A_j + A_j sum_k g(j,k)|A_k|^2.
```

The second-order mean, sum, and difference fields use the summed exponents
`lambda_j+lambda_k`; this is the homological equation required away from a
neutral point. The two modes need not have the same individual critical
acceleration.

For a finite-time, small-amplitude comparison away from onset, select

```matlab
input.weaklyNonlinear.transientModel = 'smallAmplitudeCorrection';
input.comparison.maximumRelativeCubicCorrection = 0.10;
input.comparison.maximumNonlinearToLinearRateRatio = 0.10;
```

This does not integrate the cubic equation as a saturation model. It forms
the exact operating-point linear amplitudes `A_L` and the first perturbative
correction

```text
delta A_j^(3)(t) = A_L,j(t) sum_k g(j,k)|A_k(0)|^2
                   integral_0^t exp(2 Re(lambda_k)s) ds.
```

The comparison is truncated before either `|delta A^(3)|/|A_L|` or the
instantaneous nonlinear-to-linear rate ratio exceeds its configured limit.
Use `transientModel='cubicEnvelope'` only when nonlinear feedback and a
near-onset saturation trajectory are intended.

An existing completed coefficient run can be postprocessed without repeating
the expensive cylinder solves. The fast postprocessor honors
`input.weaklyNonlinear.transientModel`: `smallAmplitudeCorrection` evaluates
the first correction above, while `cubicEnvelope` integrates the fully coupled
system using the saved `lambda_j` and `g(j,k)`. The cache path is generated
from the same user input structure as the main driver:

```matlab
result = vi_saved_small_amplitude_transient(input);
```

This produces corresponding mode-by-mode linear and nonlinear trajectories,
probe displacements, and the actual right-hand-side terms `lambda_j*A_j` and
`A_j*sum_k(g(j,k)*abs(A_k)^2)`. In `smallAmplitudeCorrection` mode it also
returns both perturbative validity diagnostics.
Initial conditions can be changed without recomputing those coefficients:

```matlab
settings.amplitudesOverH = [1e-3; 1e-6];
settings.phases = [0; pi/2];
result = vi_saved_small_amplitude_transient( ...
    input,settings);
```

An explicit MAT filename remains accepted when postprocessing outside the
driver input workflow.

Choose the number of retained amplitudes explicitly and assign each initial
complex amplitude through a magnitude and phase:

```matlab
input.numberOfModes = 1;  % or 2
input.initialConditions.amplitudesOverH = [1.0e-3; 1.0e-6];
input.initialConditions.phases = [0.0; pi/4];
```

The candidate definitions `input.modes(1)` and `input.modes(2)` remain in the
input section; only the first `input.numberOfModes` entries are retained. Thus
the example above gives `A_1(0)=10^-3` and
`A_2(0)=10^-6 exp(i*pi/4)`. The driver reports the exact `lambda_j`, wrapped
`s_j`, and reduced/full residuals for both modes at `analysisAmplitude`.

The editable full-cylinder script now defaults to an exploratory
operating-point calculation with exactly the modes entered in
`input.modes(1:input.numberOfModes)`. In that default configuration
`autoSelectSlowModes=false`, `requireSlowModesForComparison=false`, and
`stopCoefficientRunWhenSlowEnvelopeFails=false`; therefore the
`maximumSlowRateRatio=0.10` threshold is reported as a warning/validity
diagnostic instead of stopping or replacing the requested pair.

The same user driver can produce a clearly labeled exploratory trajectory when
a forced field narrowly misses the strict residual gate.  Its default
`allowExploratoryTrajectoryWithUnconvergedForcedFields=true` retains all
finite fields below `forcedExploratoryResidualTolerance=1e-5`, computes the
complete `g(j,k)` matrix, and compares the reconstructed interface with the
exact linear result.  The original `forcedSolveResidualTolerance=1e-6` is not
relaxed: `forcedSolvesValid`, `numericalCoefficientValidity`, and
`quantitativelyValid` remain false.  The exploratory ceiling also stops
residual correction once reached, avoiding repeated Krylov iterations that do
not improve the physical residual.  Quadratic resonances, nonfinite fields, or
fields above the exploratory ceiling still withhold the trajectory.

For the editable full-cylinder script, `autoSelectSlowModes=true` can keep the
first requested branch and replace a fast second branch before the expensive
full WNL coefficient pass. The default search scans `m` from the requested
mode set, radial indices `1:10`, and harmonic/subharmonic branches. For the
default inputs at `a/g0=0.70`, this avoids the fast `m2_l2_harmonic` branch
and selects the slowest admissible replacement found by the configured scan.
The selected set is saved in `output.slowModeSelection`.

For two retained modes, the complete nonresonant system is integrated:

```text
dA_j/dt = lambda_j A_j + A_j sum_k g(j,k)|A_k|^2,
j,k=1,2.
```

The default two-mode coefficient path does not solve identical second-order
problems twice. It reuses each mode's self-generated mean, the common
`q_AB=q_BA` sum field, and the conjugate difference field after evaluating its
residual with the full conjugate operator. This reduces the usual ten forced
solves to six for a nonresonant pair. The saved
`output.weaklyNonlinear.optimization` structure reports the actual reuse. For
an explicit duplicate-solve verification, set

```matlab
input.options.reuseTwoModeForcedFields = false;
```

Consequently, a mode assigned exactly zero remains zero: this nonresonant
cubic system has cross-saturation but no additive transfer term. Use a small
nonzero seed when you want to observe competition between the two modes.
The histories are returned as columns of
`output.comparison.modeAmplitudesWnl`; the full coupling matrix is
`output.comparison.gMatrixInZetaUnits`. The parsed complex starting values are
also saved in `output.initialConditions`, even when the comparison is skipped.
The modal comparison plots signed coordinates rather than magnitudes:
`Re(A_j)` from the exact linear evolution and coupled WNL evolution share the
same axes. The corresponding complex linear trajectory is saved as
`output.comparison.modeAmplitudesLinear`; the real trajectories,
linear-minus-WNL differences, RMSE, relative L2 errors, and maximum absolute
differences are saved in the `realModeAmplitude*` fields. Magnitudes remain in
the output for weak-amplitude cutoff and backward-compatible diagnostics.
For an operating-point two-mode run, the reconstructed linear reference also
contains both exact modal carriers and both exact exponential amplitudes; it is
not a first-mode-only reference.

This removes the codimension-two software restriction, but it does not make a
cubic truncation uniformly accurate at arbitrarily large growth rates. The
driver records a separate slow-envelope validity flag when
`abs(real(lambda_j))/omegaStar` exceeds
`input.weaklyNonlinear.maximumSlowRateRatio`. With
`requireSlowModesForComparison=true`, coefficients can still be evaluated as
an exploratory operating-point reduction by setting
`input.weaklyNonlinear.stopCoefficientRunWhenSlowEnvelopeFails=false`, but
the comparison is withheld.  With the default
`requireSlowModesForComparison=false`, the exploratory plot is produced and
marked non-quantitative. Always
check
spectral isolation, small interface amplitude, residual convergence, and
truncation convergence.

Set `input.weaklyNonlinear.reference='commonNeutral'` to reproduce the former
neutral-point formulation. Only in that compatibility mode must both modes be
neutral at one `aCritical`, and the linear coefficients are
`(analysisAmplitude-aCritical)*mu_j`.

The file runs the existing H/SH linear neutral-threshold preview by default.
The current code uses one shared finite-depth reduced Stokes operator for both the threshold
and growth-rate calculations. It enforces no slip at both `z=-1` and `z=+1`
with a scaled exponential basis, even when the Stokes wavenumber is large.
The reduced forcing matrix now represents the same
`cos(tau+input.forcing.phase)` convention as the primitive WNL operator.
For operating-point WNL, each full-state seed is evaluated at
`analysisAmplitude` and at its exact `gamma_j=lambda_j+i*s_j*omegaStar`.
The script prints the reduced fixed-exponent residual before attempting the
full lift.
Its `output.linear` structure contains `growthRateStar`,
`growthRatePerSecond`, `floquetMultiplier`, `harmonicIndices`,
`zetaCoefficients`, `periodicPart`, and the reconstructed
`displacementOverH`, `analysisAmplitude`, and
`accelerationOffsetFromCritical`. Set
`input.forcing.analysisAmplitude` to the desired physical value to plot a
growing or decaying transient away from neutral onset. It uses
`vi_dominant_floquet_root` to run several initial guesses, reject roots with
large singular residuals, remove Floquet aliases, and retain the accepted root
with largest real part.

The full-state null solve is branch-tracked. The driver forms an interface
seed from the requested radial shape `J_m(betaStar*r)` and, for every retained
harmonic or subharmonic mode, its reduced fixed-exponent Floquet interface
vector at the common requested acceleration. A primitive
incompressible matrix can have many smaller pressure/bulk singular directions
whose interface displacement is exactly zero, so a finite list of globally
smallest singular vectors is not used. Instead, `wnl_compute_mode` minimizes
the residual of the complete primitive-variable operator subject to a nonzero
projection onto the interface seed. The seed labels the branch; it does not
replace the full solve. Consequently the reconstructed velocity, pressure,
and displacement satisfy the full bulk, interface, wall, pressure-gauge,
volume, and chosen contact-line rows to the reported residual tolerance. It
obtains the paired adjoint from an exactly normalized bordered system containing
`A'`, the recovered direct vector, and the `Bslow*direct` normalization row;
this removes the null direction without a large penalty weight. By default,
the solver factors each uncoupled temporal-harmonic Stokes block. It combines
those factors with the direct-mode border column and descriptor-normalization
row using a scalar Schur complement. This bordered block-Jacobi preconditioner
improves the physical adjoint equation and the normalization equation in the
same GMRES step. Restartable bordered LSQR remains an automatic fallback. The
controls are

```matlab
input.options.adjointUseBlockGmres = true;
input.options.adjointGmresRestart = 30;
input.options.adjointGmresMaxCycles = 120;
input.options.adjointBlockRegularization = 1.0e-6;
input.options.adjointBlockRegularizationGrowth = 10.0;
input.options.adjointBlockRegularizationAttempts = 5;
```

For the default subharmonic cutoff `N=10`, the log should report that 22
temporal blocks are being factored and should identify the final method as
`block-GMRES` (or `block-GMRES + LSQR refinement`).
The block regularization is applied only to the preconditioner. It does not
modify the full bordered adjoint system or the reported physical residual.
Zero, nonfinite, or incorrectly normalized adjoint candidates are rejected,
including GMRES flag-2 returns with zero iterations.
A GMRES flag-1 candidate is accepted only when independent evaluation of the
full physical adjoint residual and normalization-row residual passes both
tolerances; the solver flag alone neither accepts nor rejects a mode.
The effective singular value, scaled residual, iterative-solver flag, and
iteration count are printed and stored in
`output.weaklyNonlinear.mode.tracking`.
The physical conjugate direct field uses the usual `m -> -m` and Floquet-
index reversal. Its algebraic adjoint also transforms equation rows. In the
cylinder operator, the two axis-regularity rows containing
`u_r+i*u_theta` and `u_r-i*u_theta` swap under conjugation; treating the
adjoint like a state vector omits that swap and gives a false conjugate-
adjoint residual even when the primary adjoint is converged. The factory now
provides this row map through `conjugateAdjointRows`, and the conjugate-mode
log states whether the model-specific map was applied.
For the free-contact-line branch, the reduced solver's interface coefficients
are the Schur-complement branch seed. The code first lifts the complete
velocity-pressure state and then releases the fixed interface coefficients in
a bordered full-state eigenpair correction for both the eigenvector and its
Floquet exponent. The reconstruction and adjoint systems use row and
column equilibration and a deterministic matrix-norm bound; they do not call
`normest`. Nonlinear products are not evaluated when either full-state mode
residual exceeds `input.options.modeResidualTolerance`. A difficult
prescribed-interface lift is restarted from its latest LSQR iterate up to
`input.options.directLiftMaxRestarts` times. Each attempt is then judged with
the unequilibrated full-cylinder residual used by the WNL gate. The solver
retains the best physical candidate and stops immediately when that gate is
satisfied. This matters because LSQR minimizes a scaled, regularized augmented
objective: a later restart can improve that internal objective while worsening
the physical mode residual. The complete attempt history is stored in
`output.weaklyNonlinear.mode.tracking.direct.attemptScaledResiduals`, and
`selectedAttempt` identifies the retained candidate. The convergence tolerance
is unchanged. If the mode still fails, the driver preserves
and saves the exact linear result and records the WNL failure in
`output.weaklyNonlinearFailure` instead of terminating the script.

Result saving also does not depend on `pwd`. A relative setting such as

```matlab
input.run.outputFile = 'vi_wnl_user_result.mat';
```

is resolved from the repository root. Relative subdirectories are created
automatically. The absolute destination is printed and stored in
`output.save.savedFile`; if the requested location cannot be written, the
driver attempts a repository `weakly_nonlinear/results` directory and then
MATLAB's temporary directory.

With `input.run.outputFile='auto'`, the compact filename identifies the
operating point and retained modes, for example

```text
vi_wnl_ag0-3_fHz-30_modes-m0l6-m2l2.mat
```

Initial amplitudes and phases are intentionally omitted because they do not
change the operating-point `lambda_j` or `g(j,k)`. They remain stored in
`output.userInput`, `output.input`, and `output.initialConditions`, and may be
overridden during `vi_saved_small_amplitude_transient` postprocessing.

### Reusing a coefficient-ready recovery

Mode recovery is usually more expensive than evaluating one additional set
of nonlinear coefficients. The supported driver can reuse the modes in its
previous MAT-file:

```matlab
input.run.modeRecoveryOnly = false;
input.run.reuseRecoveredModes = true;
input.run.recoveredModeFile = []; % automatically use outputFile
input.run.requireRecoveredModes = false;
input.run.allowModeRecoveryDuringCoefficientRun = false;
```

The cache is accepted only when the physical parameters, requested
acceleration, contact-line condition, operator factory, complete numerical
grid, azimuthal/radial/Floquet branches, and vector dimensions match. Its
saved direct and adjoint residuals must already pass the strict coefficient
gate. The code then rebuilds the current blocks and independently recomputes
the normalization, direct/adjoint residuals, conjugate modes, and reduced/full
exponent agreement. A recovery-only run may recover normally after a rejected
cache. A coefficient run stops immediately when no cache is accepted, because
silently rebuilding two modes can take hours. Set
`allowModeRecoveryDuringCoefficientRun=true` only when intentionally running
recovery and coefficients in one uninterrupted calculation. Cache use and its
source file are saved in `output.recoveryCache`.

### Pressure-safe forced fields

For a nonresonant quadratic field, the code solves `A*q=f`. In the cylinder,
the only coupling between adjacent Floquet harmonics comes from `Lplus` and
`Lminus`, and those matrices act only on interface-displacement columns. The
preferred solver factors each uncoupled velocity-pressure temporal block,
eliminates it, solves the resulting small block-tridiagonal interface Schur
system, and reconstructs every primitive variable. This is an exact block
elimination of the same matrix `A`: it does not weight wall rows or project
the forcing. Every reconstructed candidate is checked with the original
unequilibrated residual
`norm(A*q-f)/norm(f)`. Thus
`input.options.forcedSolveResidualTolerance` remains the acceptance criterion.

If a temporal diagonal cannot be eliminated, the fallback equilibrates
columns only, which preserves the original physical least-squares objective,
and evaluates progressively smaller QR-rank tolerances. If that also fails,
the field remains invalid and the existing descriptor-DAE completion and
regularized correction diagnostics are applied.
V38 additionally rejects a completed refinement seed that increases this
physical residual unless its full-state norm is reduced below
`input.options.forcedAlgebraicCompletionMaximumFullNormRatio` times the raw
state norm. This prevents a well-scaled forced field from being discarded for
a worse algebraic completion while retaining completion for genuinely
pressure/gauge-inflated states.

The principal controls are

```matlab
input.options.forcedUseCylinderSchur = true;
input.options.forcedCylinderSchurBlockRefinementSteps = 2;
input.options.forcedCylinderSchurReducedRefinementSteps = 2;
input.options.forcedUseRankAwareMinimumNorm = true;
input.options.forcedRankToleranceFactors = [1.0e-2,1.0e-4];
input.options.forcedTryDefaultRankTolerance = false;
input.options.forcedColumnEquilibrationFloorRatio = 1.0e-14;
input.options.forcedUseBlockGmres = false; % optional legacy seed path
input.options.forcedGmresRestart = 30;
input.options.forcedGmresMaxCycles = 120;
input.options.forcedGmresStagnationCycles = 3;
input.options.forcedGmresMinimumCycleImprovement = 1.0e-2;
input.options.forcedBlockRegularization = 1.0e-6;
input.options.forcedBlockRegularizationGrowth = 10.0;
input.options.forcedBlockRegularizationAttempts = 5;
input.options.forcedBlockMaximumInverseGain = 1.0e10;
```

The Schur diagnostics report the number of interface coordinates, reduced
dimension, any temporal-block regularization used to obtain a finite seed,
block-response and reduced-system residuals, and the final complete equation
residual. The regularized Schur path is only a seed; acceptance is still based
on the original unequilibrated forced residual. With `solveTolerance=1e-10`,
the fallback rank factors
above try `1e-12` and then `1e-14`. The first candidate meeting the unchanged
physical gate is retained.
The `development` and `balanced` profiles cap optional GMRES at 8 and 30
cycles, respectively; `final` retains the entered limit. A failed,
unavailable, or stagnant block preconditioner returns to the rank-aware
minimum-norm path.

In a two-mode calculation, a cross coefficient stops immediately when a
required mean, difference, or sum field fails this gate; later fields cannot
repair that coefficient. The opposite cross coefficient is still evaluated
when it has an independent valid path. `wnlResult.optimization` reports the
numbers of executed, shared/reused, and fail-fast-skipped forced solves.

When the WNL calculation finishes, `output.comparison` reports the
linear and WNL growth rates, their relative difference, the cubic coefficient
expressed in `zeta/h`, the predicted saturation amplitude, and both transient
envelopes. This requires `input.run.modeRecoveryOnly=false`; recovery-only
mode deliberately does not calculate `mu` or `g`. When a comparison is not
available, `output.comparison.available`, `status`, `skipReason`, and
`requiredAction` explain why. For a nontrivial selected transient near onset,
set `input.forcing.analysisAmplitude` slightly above or below the computed
critical value. The
standalone linear transient retains 20 forcing periods, but the default WNL
comparison, its error norms, and its plots use only periods 0 through 4:

```matlab
input.linear.numberOfPeriods = 20;
input.comparison.endForcingPeriod = 4.0;
input.comparison.snapshotForcingPeriod = 4.0;
```

Set `endForcingPeriod=[]` to compare through the complete linear record. The
output stores both the requested and actual sampled endpoint, as well as the
linear, WNL, and linear-minus-WNL values at that endpoint. The
comparison assumes `operators.P` is the derivative with
respect to one unit of `a/g0`. Cubic amplitudes also require the direct mode
normalization specified by
`input.comparison.zetaOverHPerUnitAmplitude` or returned as
`operatorMetadata.zetaOverHPerUnitAmplitude`.
With `input.comparison.runGrowthRateSweep=true`, the script additionally
recomputes the dominant linear growth rate over the specified acceleration
offsets and compares that curve with the WNL tangent
`Re(mu)*(a-aCritical)`. The sweep RMSE and relative L2 error are returned in
`output.comparison`. The code also stores signed differences using an explicit
`LinearMinusWnl` suffix. It prints the growth-sweep table and error norms and
creates a separate three-panel difference figure for the growth-rate sweep,
transient envelope, and reconstructed interface signal. Use
`input.comparison.printDifferenceTable` and
`input.comparison.plotDifference` to control these outputs.

The full reconstructed interface comparison is enabled by

```matlab
input.comparison.plotInterfaceDynamics = true;
input.comparison.includeSlavedHarmonics = true;
input.comparison.probeRadiusOverR = [];  % automatic modal maximum
input.comparison.probeTheta = 0.0;
```

At a fixed azimuth, it compares the same real displacement field in the two
models. The linear field is the reduced Floquet mode with its exact linear
growth rate. For one mode, the WNL field is reconstructed as

```text
zeta_WNL/h = Re{A phi + A^2 q_AA + |A|^2 q_AbarA},
```

where `A` obeys the cubic Landau equation, `q_AA` is the slaved second
azimuthal harmonic (`m_out=2m`), and `q_AbarA` is the slaved mean field
(`m_out=0`). The cubic coefficient changes the amplitude evolution; a separate
third-order spatial-shape correction is not reconstructed. Before differencing,
the arbitrary complex phase of the WNL eigenvector is aligned with the linear
Floquet carrier over the first forcing period; subsequent nonlinear phase
drift is retained. The code produces
a probe-history comparison, radial-time maps, a linear-minus-WNL map, and a
full circular-interface snapshot. Numerical arrays and error measures are in
`output.comparison.interfaceDynamics`, including `linearField`,
`wnlPrimaryField`, `wnlSecondOrderField`, `wnlTotalField`, `fieldRmse`,
`fieldRelativeL2`, and `fieldMaximumAbsoluteDifference`.

For two modes, the reconstruction adds both primary fields, both pairs of
self-generated mean/second-harmonic fields, and the cross-generated
sum/difference fields:

```text
zeta_WNL/h = Re{sum_j A_j phi_j
              + sum_j [A_j^2 q_jj + |A_j|^2 q_jbarj]
              + A_1 A_2 q_12 + A_1 conj(A_2) q_1bar2}.
```

At an operating point, the two-mode interface comparison uses both exact
linear modal carriers and both initial amplitudes as its linear reference.

The default `input.run.operatorFactory='cylinder_wnl_operators'` is included.

### Bessel-enriched radial differentiation

The Bessel-enriched implementation retains radial point values, so the complete pointwise ALE residual and
nonlinear products need no reformulation, but changes how radial derivatives
are constructed. For every selected `(m,radialIndex)` branch, the adapted
space contains `J_m(beta*r)` and the `J_|m-1|(beta*r)` and
`J_|m+1|(beta*r)` horizontal-velocity structures. Analytically differentiated
quadratic Bessel products and regular powers fill the remaining slots.

```matlab
input.numerics.Nr = 9;
input.numerics.radialGrid.type = 'besselEnriched';
input.numerics.radialGrid.maximumProductOrder = 2;
input.numerics.radialGrid.maximumConditionNumber = 1e10;
input.numerics.radialGrid.fallbackToChebyshev = false;
```

The reported basis condition number and derivative reproduction residuals must
be acceptable. The selected linear Bessel branches are differentiated exactly
within roundoff even when their radial index is high, but nonlinear products
and forced secondary fields still require an `Nr` convergence study. Compare
at least `Nr=7,9,11`; two-mode or cubic-product calculations often need more
slots. Use `radialGrid.type='chebyshev'` for the legacy radial operator.

### Multi-domain vertical differentiation

The multidomain implementation uses a three-element vertical Chebyshev grid in each fluid: a wall-layer
element, a bulk element, and an interface-layer element. The automatic thin
element widths are based on the oscillatory Stokes thickness and are capped at
20 percent of the fluid depth. Configure it with

```matlab
input.numerics.verticalGrid.type = 'multidomain';
input.numerics.verticalGrid.boundaryLayerWidthFactor = 4;
input.numerics.verticalGrid.maximumBoundaryLayerFraction = 0.20;
input.numerics.verticalGrid.pointsPerBoundaryLayer = 11;
input.numerics.verticalGrid.pointsInBulk = 15;
```

For full control, provide `lowerBreaks/lowerPoints` and
`upperBreaks/upperPoints`. For example, `lowerBreaks=[-1,-.98,-.02,0]` and
`lowerPoints=[11,15,11]`. Adjacent elements keep separate copies of their
common endpoint. At every interior radial point the operator imposes
continuity of velocity and pressure and continuity of the three vertical
velocity derivatives. The remaining pressure row retains incompressibility.
Axis and sidewall rows take precedence at the radial corners. These artificial
interface rows are linear under the shared per-fluid ALE map, so their
quadratic and cubic residual entries are zero.

The driver prints every break point, point count, Stokes thickness, first wall
spacing, and number of points inside the wall and interface Stokes layers.
`NzLower` and `NzUpper` are used only when
`verticalGrid.type='single'`. Set `input.run.modeRecoveryOnly=true` while
refining `N`, `Nr`, the element break points, and the element point counts;
this reconstructs and validates the direct/adjoint modes while skipping
expensive nonlinear products. Once mode residuals converge, set it back to
`false`, then refine `Ntheta`, `quadraticStep`, and `cubicStep` before
reporting `g`. The default nonlinear action evaluates a factor-of-two step
sequence and compares neighboring Richardson extrapolants:

```matlab
input.execution.adaptiveDirectionalSteps = true;
input.execution.quadraticStepMultipliers = [1,2,4,8];
input.execution.cubicStepMultipliers = [1,2,4];
```

The cylinder residual first removes its exact flat-interface linear part, so
the finite-difference stencil acts only on terms that can contribute to `C`
or `D`. The printed Richardson disagreement is a directional-step convergence
diagnostic; it is not a replacement for the forced-field residual gate. Set
`adaptiveDirectionalSteps=false` only for a controlled one-step comparison.
These fields are execution settings, so changing them does not invalidate a
direct/adjoint recovery cache on the same physical and numerical grid.

If recovery still fails, inspect the printed row-family table before changing
the grid again. A dominant normal-traction family points to a neutral-point or
forcing-phase mismatch; dominant bulk momentum/divergence families point to
spatial resolution or an iterative-lift failure. Do not relax
`modeResidualTolerance` to make a coefficient run continue.

## Boundary conditions used by the cylinder factory

The full operator preserves the conditions implicit in the reduced linear
branch:

- no slip at `z=-1` and `z=+1`;
- at `r=R0`, `u_r=0`, `d_r(w)=0`, and `d_r(r*u_theta)=0`;
- regular Fourier fields at `r=0`;
- continuity of all three velocity components at the interface;
- both tangential-traction conditions and the full normal-traction balance;
- the exact kinematic condition on the moving interface;
- `d_r(zeta)=0` for `boundary.contactLine='free'`, matching the
  `J'_m(beta*R0)=0` branches;
- `zeta=0` for `boundary.contactLine='pinned'` at every order.

For `m=0`, one kinematic collocation row is replaced by
`integral(zeta*r dr)=0`; this enforces fixed layer volumes and removes the
constant radial displacement excluded by the reduced positive-root Bessel
basis.

The nonlinear residual uses
`z=Z+(1+Z)*zeta` below and `z=Z+(1-Z)*zeta` above. It includes convection,
ALE time/spatial metrics, incompressibility, stress evaluation on the moving
surface, the exact graph normal and curvature, capillarity, and vibration
acting on nonlinear traction geometry. The nonlinear action is evaluated on
azimuthal and fast-time pseudospectral grids and projected onto the requested
`(m,n+s)` block.

The Stuart–Landau test has an analytic coefficient. For

```math
q_t = Jq + \mu q + g_0(q^Tq)q,
```

the selected normalization gives `mu=1` and `g=4*g0`.

## Constructing the cylinder model

Suppose `operators` contains the callbacks documented in
`vi_wnl_model.m`:

```matlab
config.omega = omega_c;
config.N = 10;
config.ndof = number_of_unknowns_per_temporal_harmonic;

model = vi_wnl_model(config,operators);
spec = model.makeSpec(m,s,'selected_mode');

opts.verbose = true;
opts.detuning = parameter_direction;
result = wnl_analyze_single_mode(model,spec,opts);
```

For several retained modes:

```matlab
specA = model.makeSpec(mA,sA,'A');
specB = model.makeSpec(mB,sB,'B');
result = wnl_analyze_mode_set(model,{specA,specB},opts);

mu = result.mu;
g = result.g;
rhs = @(T,A) wnl_rhs_landau(T,A,mu,g);
```

The reduced equations are

```math
\frac{dA_j}{dT}
=\mu_jA_j+\sum_k g_{jk}A_j|A_k|^2.
```

If `result.quadraticResonances` is nonempty, this cubic-only system is not
the correct scaling. The reported entries identify the target and the two
input modes, together with the projected quadratic coefficient.

## Required cylinder operators

### `operators.B0(spec)`

Return the singular descriptor matrix. It has identity blocks in the two
momentum rows and the kinematic `zeta_t` row. Incompressibility, pressure,
traction, velocity-continuity, wall, and contact-line rows receive zero.

Do not pass the arrays called `B`, `B1`, or `B0` in the present linear
MATLAB files. Those arrays shift temporal harmonics under vibration.

### `operators.Lhat(spec,k)`

Return the `k`th fast-time Fourier coefficient of the complete linear
operator. With

```math
L(\tau)=L_0+g_{\rm sgn}G
+a_c\cos(\tau+\phi_0)G,
```

the only nonzero temporal coefficients are

```math
L_0+g_{\rm sgn}G,\qquad
\frac{a_c}{2}e^{i\phi_0}G,\qquad
\frac{a_c}{2}e^{-i\phi_0}G.
```

`wnl_fourier_model` assembles

```math
A_{n,n'}=
i\omega_c(n+s)B_0\delta_{n,n'}-L_{n-n'}.
```

### `operators.C`, `operators.D`, `operators.CField`, and `operators.DField`

These return one spatially projected bilinear or trilinear forcing. They
must include:

- bulk convection;
- fixed-domain/ALE metric terms;
- evaluation of velocity and stress at the displaced interface;
- nonlinear kinematic terms;
- normal-vector and traction geometry;
- capillary curvature;
- the contact-line condition at every perturbation order.

The cylinder factory uses the field-level callbacks to process all temporal
harmonics together. The coefficient-level callbacks remain available for
tests. For production coefficients, compare multiple differentiation steps;
the defaults use symmetric finite-difference polarization.

### `operators.P`

This optional callback returns parameter-detuning blocks. When frequency is
varied, include

```math
-\omega_2 B_0\partial_\tau
```

in addition to explicit derivatives of the spatial operator.

## Relation to the existing linear solver

The existing functions remain useful for:

- neutral-point searches;
- Bessel roots and radial branch labels;
- `zeta_all` as an interface-mode validation target;
- analytic vertical-velocity reconstruction;
- finite-difference checks of linear growth-rate sensitivities.

The new full operator should reproduce those thresholds and interface
eigenvectors before nonlinear coefficients are trusted.

`bessel_derivative_root.m` has been added at repository level because that
dependency was missing from the public repository.

## Required validation order

1. Reproduce the selected reduced-solver neutral acceleration with the full
   `B0`/`Lhat` block.
2. Verify direct and left residuals and the acceleration sensitivity `mu`.
3. Inspect the forced mean and second-harmonic fields and their boundary
   residuals.
4. Vary the quadratic and cubic differentiation steps.
5. Converge `g` independently in temporal, radial, azimuthal, and vertical
   resolution.
6. Compare the reconstructed branch with early nonlinear DNS data.
