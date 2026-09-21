function opts = wnl_options(userOpts)
%WNL_OPTIONS Default numerical options for weakly nonlinear analysis.
%
% opts = WNL_OPTIONS()
% opts = WNL_OPTIONS(userOpts)
%
% Unknown user fields are retained so model-specific solvers can share the
% same option structure.

opts = struct();
opts.nullTolerance = 1.0e-10;
opts.resonanceTolerance = 1.0e-8;
opts.solveTolerance = 1.0e-10;
opts.fullSvdMax = 2500;
opts.modeTrackingConstraintWeight = 1.0e3;
opts.modeTrackingRegularization = 1.0e-12;
opts.modeTrackingSolveTolerance = 1.0e-9;
opts.modeTrackingMaxIterations = 10000;
opts.directLiftMaxRestarts = 3;
% At a nonzero operating-point exponent the prescribed-interface lift is
% only a Newton seed. Its own solve may stop at this residual; the corrected
% full mode must still pass the unchanged linear/coefficient gates below.
opts.directLiftSeedResidualTolerance = 5.0e-6;
opts.modeTrackingMinimumConstraint = 1.0e-8;
opts.adjointBorderedMaxRestarts = 3;
% Prefer a model-specific adjoint solve when one is available. The cylinder
% implementation eliminates full temporal blocks through their low-rank
% interface coupling; wnl_compute_mode still checks A'*left in full.
opts.adjointUseModelSolver = true;
% Equal-order primitive collocation can leave pressure-only null vectors in
% every uncoupled temporal block. The cylinder Schur adjoint removes only
% those directions with a minimum-norm SVD border before eliminating the
% block. The complete, unmodified Floquet adjoint residual is still checked.
opts.adjointCylinderSchurUseNullspaceBordering = true;
opts.adjointCylinderSchurNullTolerance = 1.0e-12;
opts.adjointCylinderSchurMaximumNullity = 8;
opts.adjointCylinderSchurBlockRefinementSteps = 2;
opts.adjointCylinderSchurFactorResidualTolerance = 1.0e-8;
opts.adjointUseBlockGmres = true;
opts.adjointGmresRestart = 30;
opts.adjointGmresMaxCycles = 120;
opts.adjointBlockRegularization = 1.0e-6;
opts.adjointBlockRegularizationGrowth = 10.0;
opts.adjointBlockRegularizationAttempts = 5;
opts.adjointBlockMaximumInverseGain = 1.0e10;
opts.modeResidualTolerance = 1.0e-6;
opts.stopOnUnconvergedMode = true;
opts.stopOnQuadraticResonance = true;
opts.forcedSolveResidualTolerance = 1.0e-6;
opts.forcedSolveConstraintTolerance = 1.0e-8;
opts.stopOnUnconvergedForcedSolve = true;
% Optional ceiling for forming explicitly exploratory coefficients from a
% finite forced field that misses the strict residual gate.  Empty preserves
% the strict behavior.  It never changes solution.valid; it only labels a
% field exploratoryUsable and, when stopping is disabled, avoids futile
% refinement once this looser ceiling has been reached.
opts.forcedExploratoryResidualTolerance = [];
opts.forcedSolveRefinementTolerance = 1.0e-9;
opts.forcedSolveMaxIterations = 10000;
opts.forcedSolveMaxRestarts = 5;
% A primitive-variable forced solution may hide a large pressure/gauge
% component in the zero-mass part of the descriptor system.  Recompute
% those algebraic variables at the nonzero quadratic right-hand side before
% iterative refinement.  This is the forced-field counterpart of the DAE
% completion used by the full eigenpair refinement.
opts.forcedUseAlgebraicCompletion = true;
opts.forcedAlgebraicCompletionRegularization = 1.0e-14;
opts.forcedAlgebraicCompletionSolveTolerance = 1.0e-10;
opts.forcedAlgebraicCompletionEquationTolerance = 1.0e-4;
opts.forcedAlgebraicCompletionMaxIterations = 1000;
opts.forcedAlgebraicCompletionMaxRestarts = 1;
% A completed state is allowed to be a refinement seed even if its raw
% equation residual temporarily grows, provided that growth is bounded.
% Every reported/accepted forced field is still selected solely by the
% unequilibrated forcing-relative residual gate above.
opts.forcedAlgebraicCompletionMaximumResidualGrowth = 100.0;
% A completion whose residual is temporarily worse is useful only when it
% actually removes a pressure/gauge-inflated state. Require at least this
% much full-state compaction (final norm / initial norm) before replacing a
% better raw forced-field seed. Residual-improving completions do not need
% to satisfy this norm-ratio test.
opts.forcedAlgebraicCompletionMaximumFullNormRatio = 0.25;
% Iterative refinement solves a row/column-equilibrated Tikhonov problem,
% with regularization defined in physical state coordinates.  This avoids
% recreating the enormous algebraic cancellation removed by completion.
opts.forcedCorrectionRegularization = 1.0e-14;
opts.forcedLineSearchMaxCuts = 8;
% Stop a cross-coefficient branch as soon as one required mean,
% difference, or sum field is invalid.  Later fields cannot rescue that
% coefficient, so computing them only wastes the dominant run time.
opts.forcedFailFast = true;
% The recommended large-cylinder forced solve is a column-equilibrated,
% rank-aware minimum-norm calculation.  Column equilibration preserves the
% physical least-squares norm, while progressively smaller QR-rank
% tolerances recover weak velocity/interface directions that a single
% absolute LSQMINNORM tolerance can incorrectly discard.
opts.forcedUseRankAwareMinimumNorm = true;
opts.forcedRankToleranceFactors = [1.0e-2,1.0e-4];
opts.forcedTryDefaultRankTolerance = false;
opts.forcedColumnEquilibrationFloorRatio = 1.0e-14;
% Direct sparse QR inside LSQMINNORM can require a dense-sized workspace
% for very large Floquet blocks. Above this dimension use iterative LSQR
% with the same column equilibration and physical-residual selection.
opts.forcedDirectMinimumNormMaxDimension = 1.0e5;
% The cylinder's vibration blocks couple neighboring temporal harmonics only
% through a small set of interface-displacement columns. Eliminate the full
% primitive-variable blocks and solve the resulting small temporal Schur
% system before falling back to a global minimum-norm calculation.
opts.forcedUseCylinderSchur = true;
% Select the primitive temporal-block factor. Column-only QR preserves the
% original least-squares equation norm and is robust to pressure-dominated
% numerical rank loss; equilibrated LU is faster on well-conditioned grids.
opts.forcedCylinderSchurFactorMethod = 'equilibratedLu';
% Equilibrate each exact primitive-variable temporal block before LU. This
% changes coordinates only; reconstructed fields are always tested in the
% original unequilibrated equations.
opts.forcedCylinderSchurUseEquilibratedBlocks = true;
opts.forcedCylinderSchurBlockRefinementSteps = 2;
opts.forcedCylinderSchurReducedRefinementSteps = 2;
% If an exact temporal block is singular, the model-specific Schur solve
% may still provide a useful seed by using a tiny shifted block inverse.
% The unshifted complete residual gate below remains the only acceptance
% criterion; an invalid regularized seed falls back to the global solve.
opts.forcedCylinderSchurBlockRegularization = 1.0e-12;
opts.forcedCylinderSchurBlockRegularizationGrowth = 10.0;
opts.forcedCylinderSchurBlockRegularizationAttempts = 5;
opts.forcedCylinderSchurBlockMaximumInverseGain = 1.0e12;
% Equal-order primitive-variable collocation can leave pressure-only null
% vectors in an otherwise valid temporal block. Complete only those vectors
% with minimum-norm borders before considering a diagonal regularization.
opts.forcedCylinderSchurUseNullspaceBordering = true;
opts.forcedCylinderSchurNullTolerance = 1.0e-12;
opts.forcedCylinderSchurMaximumNullity = 8;
% Equal-order m=0 blocks contain certified pressure-only null vectors. A
% least-squares field can therefore miss the unprojected equations solely
% in the corresponding redundant compatibility directions. When enabled,
% accept against the residual in the pressure quotient while retaining and
% bounding the removed residual separately. Nonaxisymmetric blocks and
% uncertified nullspaces always use the complete residual.
opts.forcedCylinderSchurUsePressureCompatibilityProjection = false;
opts.forcedCylinderSchurPressureCompatibilityTolerance = 1.0e-5;
% The older temporal-block GMRES path remains available as an optional
% seed generator, but is disabled by default: for the full two-mode
% cylinder it can spend thousands of Krylov iterations before the same
% minimum-norm fallback is needed.  When enabled, stagnation detection
% stops it after a few unproductive restart cycles.
opts.forcedUseBlockGmres = false;
opts.forcedGmresRestart = 30;
opts.forcedGmresMaxCycles = 120;
opts.forcedGmresStagnationCycles = 3;
opts.forcedGmresMinimumCycleImprovement = 1.0e-2;
opts.forcedBlockRegularization = 1.0e-6;
opts.forcedBlockRegularizationGrowth = 10.0;
opts.forcedBlockRegularizationAttempts = 5;
opts.forcedBlockMaximumInverseGain = 1.0e10;
% Diagnostic/production escape hatch for very large models: retain only the
% model-specific forced candidate and skip global minimum-norm/refinement
% fallbacks. A missed full residual gate remains invalid and triggers the
% ordinary coefficient fail-fast policy.
opts.forcedModelSolverOnly = false;
% Long LSQR corrections often sit at the pressure compatibility floor for
% all 10,000 iterations. Work in bounded chunks and stop only after repeated
% chunks fail to improve the original physical residual materially.
opts.forcedRefinementChunkIterations = 1000;
opts.forcedRefinementStagnationChunks = 2;
opts.forcedRefinementMinimumChunkImprovement = 1.0e-3;
opts.reuseTwoModeForcedFields = true;
% Empty evaluates the complete coefficient matrix. An N-by-2 array of
% [target,source] indices evaluates only those entries, which is useful for
% convergence and invariance checks of one suspicious coefficient.
opts.coefficientPairs = [];
opts.refineOperatingPointEigenpair = true;
opts.eigenpairRefinementTolerance = 1.0e-9;
opts.eigenpairRefinementMaxSteps = 8;
opts.eigenpairRefinementSolveTolerance = 1.0e-10;
% Solve each bordered Newton correction with the same temporal-block
% strategy used successfully by the adjoint.  LSQR is retained only as a
% short fallback when the preconditioned correction is unavailable.
opts.eigenpairUseBlockGmres = true;
opts.eigenpairGmresRestart = 30;
opts.eigenpairGmresMaxCycles = 80;
opts.eigenpairGmresAcceptanceTolerance = 1.0e-5;
opts.eigenpairLsqrFallbackMaxIterations = 3000;
opts.eigenpairLsqrFallbackMaxRestarts = 0;
% The fallback solves a Tikhonov-regularized least-squares problem in the
% physical correction norm, then column-equilibrates that complete
% objective. This both accelerates LSQR and suppresses algebraic pressure
% directions without changing the physical norm being minimized.
opts.eigenpairCorrectionRegularization = 1.0e-12;
% An inexact Newton direction is not useful merely because a line search
% finds a tiny residual decrease. Require the bordered correction equation
% itself to be solved to this relative tolerance before testing a step.
opts.eigenpairCorrectionEquationTolerance = 1.0e-1;
% A promising Newton trial can have an accurate descriptor state but an
% enormous pressure/gauge component. Recompute only the zero-mass algebraic
% variables while holding Bslow-supported and interface variables fixed.
opts.eigenpairUseAlgebraicCompletion = true;
opts.eigenpairAlgebraicCompletionRegularization = 1.0e-14;
opts.eigenpairAlgebraicCompletionSolveTolerance = 1.0e-10;
opts.eigenpairAlgebraicCompletionEquationTolerance = 1.0e-4;
opts.eigenpairAlgebraicCompletionMaxIterations = 1000;
opts.eigenpairAlgebraicCompletionMaxRestarts = 1;
opts.eigenpairLineSearchMaxCuts = 12;
opts.eigenpairStagnationRelativeImprovement = 5.0e-3;
opts.eigenpairStagnationSteps = 2;
% A residual-only Newton acceptance test can mistake an algebraic
% pressure/gauge null vector for the interfacial Floquet mode. Anchor phase
% and scale to the prescribed interface displacement (or Bslow variables
% when no interface metadata exists), and require every accepted correction
% to retain both the norm and direction of Bslow*phi.
opts.eigenpairCorrectionConstraintTolerance = 1.0e-10;
opts.eigenpairMinimumDescriptorRetention = 1.0e-2;
opts.eigenpairMinimumDescriptorOverlap = 0.10;
% Do not launch a costly full-state Newton correction from a lift that is
% still far above the ordinary physical-mode gate.
opts.eigenpairMaximumSeedResidualRatio = 5.0;
% Once at least two improving Newton corrections give a linear-ready state
% whose exponent is clearly outside the reduced/full agreement interval,
% more algebraic work cannot make that discretization coefficient-ready.
opts.stopEigenpairOnReferenceMismatch = true;
opts.eigenpairReferenceMismatchAbortFactor = 1.25;
opts.eigenpairReferenceMismatchMinimumSteps = 2;
opts.eigenpairReferenceMismatchMaximumLastStepRatio = 0.25;
% Internal continuation bookkeeping; ordinary callers leave these values.
opts.eigenpairReferenceLambda = [];
opts.eigenpairPriorAcceptedSteps = 0;
% A runtime profile may deliberately cap the first eigenpair pass.  When
% that pass improves the residual but does not reach the coefficient gate,
% continue from its best eigenpair with these larger algebraic limits.  No
% spatial or Floquet resolution is changed by this continuation.
opts.autoEscalateEigenpairRefinement = true;
opts.eigenpairEscalationMaxSteps = 8;
opts.eigenpairEscalationMaximumResidualRatio = 1.0e4;
opts.coefficientModeResidualTolerance = 1.0e-8;
% Axisymmetric harmonic/subharmonic branches with a real Floquet exponent
% are one-dimensional real eigenspaces. Numerical recovery may leave a tiny
% imaginary exponent and phase defect; project such a mode onto its physical
% self-conjugate representative before evaluating nonlinear products.
opts.realModeEigenvalueTolerance = 1.0e-4;
opts.realModeMinimumConjugateOverlap = 0.99;
% Real quadratic and cubic responses are projected onto the exact fixed
% point of physical conjugation. Reject the coefficient if the unprojected
% field is too far from that subspace; projection must remain a numerical
% cleanup, not a repair of a physically inconsistent solve.
opts.realFieldConjugacyTolerance = 1.0e-5;
% A coefficient in a real modal equation must itself be real. This gate is
% applied after explicit self-conjugate projection and catches broken
% conjugation bookkeeping or an inadequately resolved nonlinear action.
opts.realCoefficientImaginaryTolerance = 1.0e-7;
opts.enforceReducedFullEigenvalueAgreement = true;
opts.maximumReducedFullEigenvalueRelativeMismatch = 0.05;
opts.maximumReducedFullEigenvalueAbsoluteMismatch = 1.0e-3;
opts.checkDirectModeResidual = true;
opts.checkAdjointModeResidual = true;
opts.verbose = true;
opts.direct = [];
opts.left = [];
opts.normalizeDirect = [];
opts.detuning = [];

if nargin == 0 || isempty(userOpts)
    return;
end

names = fieldnames(userOpts);
for j = 1:numel(names)
    opts.(names{j}) = userOpts.(names{j});
end
end
