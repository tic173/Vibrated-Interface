# Current clean release (V42)

The distributed package now contains one supported cylinder input driver:

`examples/vi_wnl_user_run_full_cylinder.m`

Obsolete version-numbered user drivers and their incremental release notes are
not included in the clean archive.

## Low-rank temporal Schur solve for cylinder forced fields

- The V41 diagnostics showed that the nonlinear directional actions were
  step-converged, while the remaining forced residual was concentrated in
  homogeneous sidewall and interface-continuity equations. The global
  minimum-norm solve was trading small boundary violations against the much
  larger bulk right-hand side; further directional-step or rank-tolerance
  changes could not correct that behavior.
- The cylinder now uses the exact low-rank structure of vertical vibration.
  `Lplus` and `Lminus` couple neighboring Floquet harmonics only through the
  interface-displacement columns. Each uncoupled primitive-variable temporal
  block is factored once, eliminated, and a small block-tridiagonal Schur
  system is solved for the coupled interface coordinates. The complete
  velocity-pressure field is then reconstructed.
- This elimination is algebraically equivalent to the original full Floquet
  equations. It does not project the nonlinear forcing, weight boundary rows,
  relax a tolerance, or discard pressure equations. The reconstructed field
  must still pass `norm(A*q-f)/norm(f) <= 1e-6` in the original,
  unequilibrated operator before it can enter a cubic coefficient.
- The console and saved forced-field result report the Schur dimension,
  number of active coupling columns, block-response residual, reduced-system
  residual, factorization time, and final complete-equation residual. If a
  temporal diagonal cannot be eliminated, the existing rank-aware global
  solver remains a checked fallback.
- A coefficient run now stops immediately on a missing or incompatible mode
  cache unless `allowModeRecoveryDuringCoefficientRun=true` is explicitly
  selected. This prevents an accidental cache-path change from silently
  repeating several hours of two-mode recovery.

The separate asymptotic warning remains unchanged: the supplied
`m=2,l=2,s=0` operating-point mode has a rate that is not slow relative to
the forcing frequency. Numerically converged coefficients make an exploratory
two-mode reduced system available, but do not by themselves make that branch
a controlled weakly nonlinear amplitude.

## Exact nonlinear remainder and adaptive directional steps

- V40 established that changing the minimum-norm rank tolerance did not
  reduce the harmonic and cross-field residual floors: the same candidates
  remained at forcing-relative residuals between about `1.9e-6` and
  `5.9e-6`. The defect was therefore in the assembled nonlinear right-hand
  side, not in the number of forced-solver iterations.
- The cylinder action now subtracts every analytically linear flat-interface
  contribution before taking a quadratic or cubic directional derivative:
  flat time, pressure, viscosity, incompressibility, velocity continuity,
  kinematics, stress, gravity, capillarity, and artificial element-matching
  rows. These terms have exactly zero second and third Frechet derivatives;
  removing them changes neither physical nonlinear operator but prevents
  cancellation of large `O(h)` values from contaminating the desired
  `O(h^2)` or `O(h^3)` action.
- Field-level actions evaluate a factor-of-two step sequence and use
  neighboring Richardson estimates to select a cancellation-safe result.
  The chosen step range and disagreement are printed. The base steps remain
  editable, and setting `adaptiveDirectionalSteps=false` restores one-step
  evaluation.
- Failed forced fields now print their largest individual cylinder equations,
  including radial/vertical indices and whether each row is bulk momentum,
  incompressibility, a physical boundary, an interface condition, or a
  multi-domain matching equation.
- Directional-step controls are execution settings and are deliberately
  excluded from the recovered-mode signature. A coefficient-ready V40 mode
  cache on the same physical parameters, acceleration, grid, and branches can
  therefore be reused without repeating direct/adjoint recovery.

Numerical convergence of the forced fields and asymptotic validity remain
separate requirements. In the supplied `a/g0=0.7` example the
`m=2,l=2,s=0` exponent has `|Re(lambda)|/omegaStar ~= 0.194`, above the
configured `0.10` slow-envelope limit. Numerical solver improvements can
remove a coefficient obstruction, but they cannot make that fast branch a controlled
weakly nonlinear slow amplitude.

## Rank-aware two-mode validity and forced fields

- Nonresonant `A*q=f` fields now use a column-equilibrated rank-aware
  minimum-norm solve by default. Column scaling changes variables but not the
  physical residual norm. The solver tries successively smaller QR-rank
  tolerances and stops only when the original unequilibrated
  `norm(A*q-f)/norm(f)` gate passes.
- The generic `solveTolerance` is no longer used as the only absolute
  `lsqminnorm` rank decision. This targets the V39 harmonic-mode residual
  floors of `1.9e-6`--`5.9e-6`, where weak physical directions could be
  discarded just above the required `1e-6` gate. No forcing projection or
  tolerance relaxation is performed; failure still produces `NaN`
  coefficients.
- Temporal-block GMRES remains available but is disabled in the editable
  large-cylinder input. When enabled, it reports `improved-seed` separately
  from `passed-gate` and terminates after several stagnant restart cycles.
  This removes the misleading V39 `accepted=1` wording and prevents a failed
  seed calculation from consuming thousands of iterations before fallback.
- The recovery cache defaults to the selected output MAT-file when
  `recoveredModeFile=[]`, preventing an output rename from silently causing
  another multi-hour direct/adjoint recovery.
- Numerical coefficient validity and slow-envelope validity are now saved
  separately. Modes at any requested acceleration can still be evaluated as
  an exploratory operating-point reduction, but a linear/WNL transient is
  withheld when `|Re(lambda_j)|/omegaStar` exceeds the configured slow-rate
  limit.

## V39 scaled temporal-block solver for forced fields

- Nonresonant mean, sum, difference, and second-harmonic systems now try a
  row- and column-equilibrated temporal-block GMRES solve before the generic
  minimum-norm/LSQR path. Each diagonal primitive-variable Stokes block is
  factored with the same adaptive regularization and finite-inverse checks
  used by the successful direct/adjoint mode solvers.
- When no model-specific direct solve is supplied, GMRES starts from zero.
  The expensive minimum-norm solve is skipped only when the resulting field
  passes the unchanged unequilibrated forcing-relative residual gate. If it
  does not, the previous minimum-norm, DAE-completion, and regularized-LSQR
  paths remain active and the best physical-residual candidate is retained.
- The console and MAT output now distinguish scaled and physical GMRES
  residuals, factorization and Krylov time, block regularization range,
  inverse-gain probes, and whether the minimum-norm fallback ran.
- This targets the V38 two-mode diagnostic in which the harmonic-mode LSQR
  correction equations had residual ratios between 38 and 46: those
  directions were not approximate solutions of the forced equations, so
  adding iterations could not repair the fields.

## V38 selective forced-field DAE completion

- A nonzero-right-hand-side algebraic completion is no longer adopted as the
  next forced-field refinement seed merely because its residual remains below
  the loose seed-growth bound. If completion increases the physical residual,
  it must also reduce the full-state norm by at least a factor of four,
  demonstrating that it actually removed a pressure/gauge-inflated state.
- This fixes the two-mode failure in which ordinary harmonic-mode fields with
  state norms of order 1e2 were changed from residuals of order 1e-3 to
  residuals of order one without any norm reduction. V38 keeps the better raw
  field and starts pressure-safe correction from it.
- The console and saved completion diagnostics now state whether selection was
  caused by residual improvement or substantial algebraic compaction. The
  unequilibrated forcing-relative acceptance gate remains 1e-6.
- A regression test covers both rejection of a noncompacting harmful
  completion and acceptance of a genuinely pressure-inflated completion.

## Pressure-safe two-mode forced fields

- The descriptor-DAE algebraic completion now accepts a nonzero right-hand
  side and is applied when an ordinary nonresonant mean, sum, difference, or
  second-harmonic correction misses its residual gate. Already-valid fields
  keep the faster legacy path. Completion holds the velocity/interface
  descriptor variables fixed while recomputing zero-mass pressure, gauge,
  boundary, and constraint variables at the actual quadratic forcing.
- Forced residual corrections now use row/column equilibration with a
  Tikhonov norm defined in physical state coordinates. A line search and the
  original unequilibrated forcing-relative residual select the returned
  field; the `1e-6` forced-field gate is unchanged.
- Each forced field reports its full-state norm, descriptor norm, DAE
  completion history, correction residuals, and wall-clock time. This makes
  pressure cancellation distinguishable from spatial/truncation error.
- Cross-coupling branches stop as soon as a required mean, difference, or sum
  field fails. The reverse coefficient is still attempted when it has a
  viable independent path. Output records executed, reused, and fail-fast
  skipped solves instead of assuming all ten baseline solves ran.

## Recovery restart

- A coefficient run can reuse direct/adjoint modes stored by an earlier
  recovery run through `input.run.reuseRecoveredModes` and
  `input.run.recoveredModeFile`.
- Reuse requires matching physical parameters, acceleration, boundary model,
  operator factory, complete numerical grid, mode branches, vector sizes, and
  strict saved residuals. The current operator then recomputes normalization,
  direct/adjoint residuals, conjugates, and reduced/full exponent agreement.
  A cache therefore saves recovery time without bypassing an acceptance gate.

## Conjugate-adjoint correction

- The physical conjugate of a cylinder direct mode still uses the usual
  Floquet map `m -> -m`, `s -> wrap(-s)`, and, for `s=1/2`, `n -> -n-1`.
  Its algebraic adjoint now additionally uses the cylinder equation-row map.
- At the cylinder axis, the regularity equations `u_r+i*u_theta` and
  `u_r-i*u_theta` exchange their stored `ur`/`ut` row locations under
  complex conjugation. Earlier releases conjugated the adjoint as though it
  were a state vector. This left the direct conjugate accurate but produced
  a spurious conjugate-adjoint residual near the ordinary mode gate.
- `cylinder_wnl_operators` now supplies the missing axis-row permutation,
  `vi_wnl_model` passes it to the generic WNL model, and
  `wnl_conjugate_mode` records whether the model-specific row map was used.
  An operator-level regression test checks both the direct and adjoint
  conjugacy identities.

## Runtime improvements

- A descriptor-DAE algebraic completion now repairs the specific failure
  observed in V34: a block-GMRES trial can have a very small full equation
  residual and the correct `Bslow` direction, yet be rejected because an
  enormous pressure/gauge component drives descriptor retention toward zero.
- The completion holds every nonzero-`Bslow` column and every prescribed
  interface coefficient fixed. It recomputes only pressure, gauge, boundary,
  and constraint variables through a regularized minimum-norm least-squares
  solve, then reapplies the unchanged phase, descriptor, residual, and
  reduced/full eigenvalue gates.
- Development disables this extra solve, balanced permits 500 iterations,
  and final permits the entered 1000 iterations and one restart. The saved
  eigenpair history records whether completion was attempted and accepted.

- The pressure-safe eigenpair fallback now adds a physical Tikhonov norm to
  the row-equilibrated bordered equations before column equilibration. This
  retains the intended physical minimum-norm selection while restoring the
  conditioning advantage that was lost in V33's unscaled fallback.
- A Newton correction is no longer accepted solely because a line search
  finds a small decrease. Its full bordered equation must first pass
  `eigenpairCorrectionEquationTolerance` (default `0.1`). The reported V33
  corrections had residuals near one and therefore consumed eight expensive
  steps without approaching the mode gate; V34 stops that work immediately.
- A failed correction-equation solve has its own stop reason and cannot
  trigger automatic Newton continuation. GMRES/LSQR logs distinguish an
  unavailable algebraic correction from a branch-preserving physical trial.
- The editable l=6 example now requests `Nr=16`. The balanced profile still
  caps it at 12, while selecting `final` now performs a genuine radial
  refinement together with the entered `N=11` and 37-point vertical grids.

- Direct and adjoint residuals are now normalized by the descriptor-supported
  physical state, not by the complete velocity-pressure vector. Pressure,
  gauge, and other zero-mass variables remain in every residual equation but
  cannot inflate the denominator and manufacture false convergence.
- Saved residual diagnostics report the full-to-physical state-norm ratio, so
  an algebraically inflated candidate is explicit rather than hidden behind a
  tiny normalized residual.
- The eigenpair LSQR fallback minimizes a regularized physical correction
  norm. Column equilibration is applied to the complete objective, so it
  accelerates convergence without rewarding enormous pressure corrections.
- A rejected GMRES correction is used as the LSQR warm start only when its
  best trial retained the physical descriptor norm and overlap. A
  pressure-inflated GMRES direction now restarts physical-coordinate LSQR
  from zero.
- Relative MAT-file output paths are resolved from the repository containing
  the driver rather than MATLAB's current directory. Missing subdirectories
  are created, and repository-results/temporary fallbacks preserve the output
  if the requested destination is unavailable.

- `development`, `balanced`, and `final` runtime profiles make cost/accuracy
  choices explicit.
- The editable driver now defaults to `balanced`; high-radial-index modes were
  frequently only linear-ready on the development grid.
- Development substantially reduces the Floquet, radial, automatic vertical,
  and iterative-solver sizes; it is a quick feasibility pass and can stop at
  a residual gate before evaluating coefficients.
- Earlier releases disabled the unscaled experimental forced-field GMRES
  path. V39 replaces it with the equilibrated, adaptively regularized method
  described above.
- Mode lift, full eigenpair refinement, adjoint, forced solves, self/cross
  coefficients, nonlinear actions, and total wall-clock time are reported.
- Independent nonlinear temporal samples can use an existing Parallel
  Computing Toolbox pool. Pool creation remains an explicit user choice.
- Development uses the minimum de-aliased nonlinear temporal sampling,
  balanced uses 1.5x oversampling, and final retains the original 2x
  oversampling.
- A recovery-only run skips every O(A^2) solve and cubic action.
- A profile-limited eigenpair solve that improves but misses the coefficient
  gate now continues from its best state with separate escalation limits. It
  does not repeat the prescribed-interface lift or change the discretization.
- The escalation step count is now a total Newton-step budget. The earlier
  interpretation could run five balanced-profile steps followed by eight
  additional steps.
- Full-eigenpair corrections now use temporal-block preconditioned GMRES with
  a short LSQR fallback and residual-decreasing line search.
- GMRES and LSQR are now staged rather than ranked only by their global
  bordered residuals. The GMRES correction is first tested against the full
  physical operator; LSQR is run only if that descriptor-preserving trial
  fails to decrease the mode residual.
- Development and balanced profiles substantially cap direct-correction
  Krylov work. Development disables the LSQR fallback; balanced uses at most
  15 GMRES cycles and 750 LSQR iterations per Newton step.
- The default prescribed-interface seed target is tightened from `1e-5` to
  `5e-6`. Full-eigenpair refinement is skipped when the lift is still more
  than five times above the ordinary physical-mode residual gate, avoiding a
  long correction from a demonstrably unsafe starting state.
- Rejected corrections now report the physical trial residual and descriptor
  retention/overlap, making the reason for rejection observable.
- The eigenvector phase/scale border is now constructed from the prescribed
  interface-displacement coefficients. The previous full-state inner product
  was dominated by pressure and could make the bordered correction point away
  from the requested interfacial branch.
- When explicit interface DOFs are unavailable, the phase border uses only the
  `Bslow` descriptor variables, excluding pressure and algebraic constraints.
- The safeguarded line search now allows twelve halvings instead of four, so a
  large inexact correction can be damped back onto the physical branch.
- GMRES and LSQR physical trial outcomes are printed separately, including the
  best attempted residual and its descriptor retention/overlap.
- The eigenpair phase/scale row is now imposed algebraically after an
  iterative correction. This prevents a globally small LSQR residual from
  hiding an order-one error in the single normalization row.
- Residual reduction is accepted only while `Bslow*phi` retains sufficient
  norm and overlap with the lifted interfacial branch. Pressure/gauge null
  vectors with zero slow-time descriptor content are rejected and cannot be
  passed to the Landau adjoint normalization.
- A zero-vector candidate is assigned an infinite normalized residual rather
  than the misleading value zero.
- A nonzero-exponent prescribed-interface lift can stop once it reaches the
  configurable Newton-seed tolerance. Final mode and coefficient residual
  gates are unchanged.
- Refinement stops early when a linear-ready corrected mode already proves a
  reduced/full exponent mismatch. Further algebraic refinement on that same
  grid would not make its nonlinear coefficients valid.
- Recovery-only output now distinguishes `LINEAR-READY ONLY` from
  `COEFFICIENT-READY RECOVERY` and saves the per-mode readiness diagnostics.
- The primary adjoint is solved directly to the coefficient-level tolerance.
  A mapped conjugate adjoint is required to pass the ordinary linear gate,
  while its conjugate direct mode retains the strict coefficient gate in a
  nonresonant cubic reduction.
- The supported editable driver starts in recovery-only mode; coefficient and
  interface-comparison work must be explicitly enabled after readiness passes.
- The two-mode coefficient path reuses both self-generated mean fields, the
  common sum field, and (after a full-operator residual check) the conjugate
  difference field. A nonresonant pair therefore uses six forced solves
  instead of ten.

## Accuracy safeguards retained

- Full primitive-variable eigenvector/Floquet-exponent refinement.
- Explicit reduced/full operating-point eigenvalue-agreement gate.
- Separate linear and coefficient-level mode residual gates.
- Forcing-relative residual gates for every slaved field.
- Separate diagnostics for quadratic resonance and forced-solve failure.
- Weak-amplitude stopping for divergent cubic envelope models.

Development and balanced results are preliminary. Reported coefficients must
be repeated with the final profile and checked for convergence in temporal,
radial, vertical, azimuthal, and directional-differentiation resolution.
