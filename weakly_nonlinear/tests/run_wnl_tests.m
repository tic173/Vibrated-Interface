function run_wnl_tests()
%RUN_WNL_TESTS Basic regression tests for the weakly nonlinear module.

thisFile = mfilename('fullpath');
moduleRoot = fileparts(fileparts(thisFile));
repositoryRoot = fileparts(moduleRoot);
addpath(moduleRoot);
addpath(fullfile(moduleRoot, 'examples'));
addpath(repositoryRoot);

wnl_demo_stuart_landau();
test_subharmonic_conjugation();
test_subharmonic_product();
test_shifted_floquet_block();
test_shifted_spec_algebra();
test_periodic_nonlinear_shift();
test_seeded_mode_tracking();
test_prescribed_lift_keeps_best_physical_residual();
test_descriptor_residual_rejects_algebraic_inflation();
test_minimum_norm_algebraic_completion();
test_nonzero_rhs_algebraic_completion();
test_forced_completion_seed_selection();
test_rank_aware_forced_minimum_norm();
test_exploratory_forced_field_gate();
test_cylinder_temporal_schur_forced();
test_temporal_block_gmres_forced();
test_bordered_adjoint();
test_temporal_block_gmres_adjoint();
test_temporal_block_gmres_eigenpair();
test_physical_coordinate_lsqr_eigenpair();
test_inaccurate_seed_skips_eigenpair();
test_nonfinite_mode_gate();
test_reference_eigenvalue_gate();
test_bessel_derivative_roots();
test_floquet_forcing_phase();
test_reduced_neutral_mode_consistency();
test_multidomain_chebyshev_grid();
test_bessel_enriched_radial_grid();
test_cylinder_operator_smoke();
test_interface_dynamics_reconstruction();
test_two_mode_initial_conditions_and_reconstruction();
test_analysis_amplitude_resolution();
test_comparison_time_window();
test_small_amplitude_cubic_transient();
test_runtime_profiles();
test_two_mode_shared_forced_fields();
test_cross_forced_fail_fast();
test_recovered_mode_cache();
test_output_save_path_resolution();
test_automatic_output_filename();
test_automatic_mode_labels();
fprintf('All WNL tests passed.\n');
end

function test_small_amplitude_cubic_transient()
time = linspace(0,0.5,51).';
A0 = [1.0e-2;2.0e-2i];
lambda = [0.4;-0.2+0.1i];
g = [-3,0.5;0.2,-1];
limits = struct('maximumAmplitude',1, ...
    'maximumRelativeCorrection',1, ...
    'maximumNonlinearToLinearRateRatio',1, ...
    'linearRateFloor',1.0e-12);
result = vi_cubic_transient_correction(time,A0,lambda,g,limits);
expectedLinear = exp(time*lambda.') .* A0.';
integrals = [(expm1(2*real(lambda(1))*time)/(2*real(lambda(1)))), ...
    (expm1(2*real(lambda(2))*time)/(2*real(lambda(2))))];
expectedMultiplier = (integrals.*abs(A0.').^2)*g.';
expectedCorrection = expectedLinear.*expectedMultiplier;
assert(result.completedRequestedWindow);
assert(norm(result.linearAmplitudes-expectedLinear,'fro') < 1.0e-13);
assert(norm(result.cubicCorrections-expectedCorrection,'fro') < 1.0e-13);
assert(norm(result.correctedAmplitudes- ...
    (expectedLinear+expectedCorrection),'fro') < 1.0e-13);

limits.maximumRelativeCorrection = 1.0e-5;
limited = vi_cubic_transient_correction(time,A0,lambda,g,limits);
assert(~limited.completedRequestedWindow);
assert(numel(limited.time) < numel(time));
assert(any(strcmp(limited.stopReasons, ...
    'relative cubic-correction limit')));
end

function test_automatic_mode_labels()
assert(strcmp(vi_wnl_mode_label( ...
    struct('m',0,'radialIndex',6,'s',0.5)), ...
    'm0_l6_subharmonic'));
assert(strcmp(vi_wnl_mode_label( ...
    struct('m',2,'radialIndex',2,'s',0)), ...
    'm2_l2_harmonic'));
assert(strcmp(vi_wnl_mode_label( ...
    struct('m',3,'radialIndex',4,'s',0.25)), ...
    'm3_l4_s0p25'));
assert(strcmp(vi_wnl_mode_label( ...
    struct('m',1,'radialIndex',3,'s',-0.5)), ...
    'm1_l3_subharmonic'));
end

function test_automatic_output_filename()
input.dimensional.frequencyHz = 30;
input.dimensional.h = 22.0e-3;
input.dimensional.R = 35.0e-3;
input.forcing.analysisAmplitude = 2;
input.numberOfModes = 2;
input.modes(1) = struct('m',2,'radialIndex',6,'s',0.5);
input.modes(2) = struct('m',2,'radialIndex',2,'s',0);
input.initialConditions.amplitudesOverH = [1.0e-3;1.0e-6];
input.initialConditions.phases = [0;pi/2];
input.execution.profile = 'final';
input.numerics = struct('N',11,'Nr',16,'NzLower',37, ...
    'NzUpper',37,'Ntheta',32);
input.boundary.contactLine = 'free';
input.comparison.endForcingPeriod = 4;
filename = vi_wnl_output_filename(input);
expected = 'vi_wnl_ag0-2_fHz-30_modes-m2l6-m2l2.mat';
assert(strcmp(filename,expected));
assert(strcmp(filename,vi_wnl_output_filename(input)));
input.initialConditions.amplitudesOverH = [0.2;0.3];
input.initialConditions.phases = [pi/3;-pi/4];
assert(strcmp(filename,vi_wnl_output_filename(input)));
input.forcing.analysisAmplitude = 2.5;
assert(~strcmp(filename,vi_wnl_output_filename(input)));
end

function test_rank_aware_forced_minimum_norm()
% A single absolute QR-rank tolerance discards the weak direction of this
% nearly dependent system.  The second rank-aware attempt must recover it,
% and validity must still be based on the original equation residual.
delta = 1.0e-6;
A = sparse([1,1;0,delta]);
% This right-hand side has a material component in the weak singular
% direction. A rank-one truncation leaves an O(delta) residual, whereas
% the tighter second rank decision recovers the complete solution.
exact = [2;-1];
forcing = A*exact;
model.block = @(spec) struct('A',A,'Bslow',speye(2)); %#ok<NASGU>
spec = wnl_spec(0,0,0,2,'rank_aware_forced_test');
opts = struct('verbose',false,'solveTolerance',1.0e-4, ...
    'forcedUseBlockGmres',false, ...
    'forcedUseRankAwareMinimumNorm',true, ...
    'forcedRankToleranceFactors',[1,1.0e-6], ...
    'forcedTryDefaultRankTolerance',false, ...
    'forcedSolveResidualTolerance',1.0e-10, ...
    'forcedSolveMaxRestarts',0, ...
    'forcedUseAlgebraicCompletion',false);
solution = wnl_solve_forced(model,spec,forcing,{},opts);
diagnostic = solution.rankAwareDiagnostics;
assert(diagnostic.attempted && diagnostic.numberOfAttempts == 2);
assert(diagnostic.selectedAttempt == 2);
assert(diagnostic.passedPhysicalGate && solution.valid);
assert(solution.relativeEquationResidual < 1.0e-10);
assert(norm(solution.vector-exact)/norm(exact) < 1.0e-8);
end

function test_exploratory_forced_field_gate()
% A finite least-squares field may be retained for an explicitly labeled
% exploratory coefficient without changing its strict validity flag.
model.block = @(spec) struct('A',sparse(1,1),'Bslow',sparse(1,1)); %#ok<NASGU>
spec = wnl_spec(0,0,0,1,'exploratory_forced_test');
opts = struct('verbose',false,'forcedUseBlockGmres',false, ...
    'forcedUseRankAwareMinimumNorm',false, ...
    'forcedSolveResidualTolerance',0.5, ...
    'forcedExploratoryResidualTolerance',1.1, ...
    'stopOnUnconvergedForcedSolve',false, ...
    'forcedSolveMaxRestarts',0, ...
    'forcedUseAlgebraicCompletion',false);
solution = wnl_solve_forced(model,spec,1,{},opts);
assert(~solution.valid);
assert(solution.exploratoryUsable);
assert(solution.relativeEquationResidual == 1);
assert(solution.solveDiagnostics.totalIterations == 0);
assert(solution.solveDiagnostics.exploratoryEarlyStop);
end

function test_cylinder_temporal_schur_forced()
% Neighboring temporal blocks are coupled only through one active column.
% The reduced interface Schur solve must reproduce the complete global
% Floquet solution without reweighting or dropping any equation.
ndof = 5;
blocks.B0 = speye(ndof);
blocks.L0 = sparse(diag([1.1,1.7,2.3,3.2,4.1]));
blocks.L0(2,1) = 0.13;
blocks.L0(4,3) = -0.09;
blocks.Lplus = sparse([2,4],[5,5],[0.21,-0.08],ndof,ndof);
blocks.Lminus = sparse([1,3],[5,5],[-0.17,0.11],ndof,ndof);
omega = 1.9;
spec = wnl_spec(3,0.5,[-2,-1,0,1],ndof,'schur_forced_test');
spec.lambda = 0.37+0.04i;
config.omega = omega;
config.N = 1;
config.ndof = ndof;
config.mass = blocks.B0;
config.linearFourier = @(requestedSpec,k) ...
    schur_test_fourier(blocks,requestedSpec,k);
config.blockSolve = @(requestedSpec,forcing,options) ...
    vi_cylinder_wnl_forced_schur_solve( ...
    blocks,omega,requestedSpec,forcing,options);
model = wnl_fourier_model(config);
assembled = model.block(spec);
exact = (1:ndof*numel(spec.n)).' + ...
    1i*(ndof*numel(spec.n):-1:1).';
forcing = assembled.A*exact;
opts = struct('verbose',false,'forcedUseCylinderSchur',true, ...
    'forcedUseBlockGmres',false, ...
    'forcedUseRankAwareMinimumNorm',false, ...
    'forcedUseAlgebraicCompletion',false, ...
    'forcedSolveResidualTolerance',1.0e-11, ...
    'forcedSolveMaxRestarts',0);
solution = wnl_solve_forced(model,spec,forcing,{},opts);
assert(solution.modelSolveDiagnostics.attempted);
assert(solution.modelSolveDiagnostics.available);
assert(solution.modelSolveDiagnostics.passedPhysicalGate);
assert(solution.modelSolveDiagnostics.numberOfActiveColumns == 1);
assert(solution.modelSolveDiagnostics.reducedDimension == numel(spec.n));
assert(solution.relativeEquationResidual < 1.0e-11);
assert(norm(solution.vector-exact)/norm(exact) < 1.0e-10);
end

function value = schur_test_fourier(blocks,requestedSpec,k) %#ok<INUSD>
if k == 0
    value = blocks.L0;
elseif k == 1
    value = blocks.Lplus;
elseif k == -1
    value = blocks.Lminus;
else
    value = sparse(size(blocks.B0,1),size(blocks.B0,2));
end
end

function test_nonzero_rhs_algebraic_completion()
% Forced O(A^2) fields have A*q=f rather than A*q=0. Verify that the same
% descriptor partition removes algebraic inflation at a nonzero right-hand
% side without changing the dynamic variables.
A = sparse(5,5);
A(1,:) = [1,0,1,1,0];
A(2,:) = [0,1,0,1,1];
Bslow = sparse(diag([1,1,0,0,0]));
inflation = 1.0e12;
state = [3;4;inflation;1-inflation;inflation];
forcing = [4;5;0;0;0];
opts = struct('eigenpairAlgebraicCompletionMaxIterations',100, ...
    'eigenpairAlgebraicCompletionMaxRestarts',0, ...
    'eigenpairAlgebraicCompletionRegularization',1.0e-14, ...
    'eigenpairAlgebraicCompletionSolveTolerance',1.0e-12, ...
    'eigenpairAlgebraicCompletionEquationTolerance',1.0e-8);
[completed,information] = wnl_complete_algebraic_state( ...
    A,Bslow,state,[],opts,forcing);
assert(information.attempted && information.available && ...
    information.valid);
assert(norm(Bslow*(completed-state)) < 1.0e-12);
assert(norm(A*completed-forcing) < 1.0e-8);
assert(norm(completed) < 1.0e-6*norm(state));
end

function test_minimum_norm_algebraic_completion()
% Hold two dynamic variables fixed and remove a large algebraic null-space
% component without changing the descriptor state.
A = sparse(5,5);
A(1,:) = [1,0,1,1,0];
A(2,:) = [0,1,0,1,1];
Bslow = sparse(diag([1,1,0,0,0]));
inflation = 1.0e12;
state = [1;2;inflation;-1-inflation;-1+inflation];
opts = struct('eigenpairAlgebraicCompletionMaxIterations',100, ...
    'eigenpairAlgebraicCompletionMaxRestarts',0, ...
    'eigenpairAlgebraicCompletionRegularization',1.0e-14, ...
    'eigenpairAlgebraicCompletionSolveTolerance',1.0e-12, ...
    'eigenpairAlgebraicCompletionEquationTolerance',1.0e-8);
[completed,information] = wnl_complete_algebraic_state( ...
    A,Bslow,state,[],opts);
assert(information.attempted && information.available && ...
    information.valid);
assert(norm(Bslow*(completed-state)) < 1.0e-12);
assert(norm(A*completed) < 1.0e-8);
assert(norm(completed) < 1.0e-6*norm(state));
assert(information.fullNormReductionFactor < 1.0e-6);
end

function test_forced_completion_seed_selection()
% A completed field may satisfy the loose seed-growth bound while being a
% much worse refinement seed. Reject it unless it either improves the raw
% equation residual or substantially removes pressure/gauge inflation.
opts = struct('verbose',false,'forcedSolveMaxRestarts',0, ...
    'forcedAlgebraicCompletionMaxIterations',100, ...
    'forcedAlgebraicCompletionMaxRestarts',0, ...
    'forcedAlgebraicCompletionRegularization',1.0e6, ...
    'forcedAlgebraicCompletionSolveTolerance',1.0e-12, ...
    'forcedAlgebraicCompletionEquationTolerance',1.0e-8, ...
    'forcedAlgebraicCompletionMaximumFullNormRatio',0.25);

A = sparse([1,1;0,0]);
Bslow = sparse(diag([1,0]));
rawState = [1;1];
forcing = [2.001;0];
model.block = @(spec) forced_seed_test_block( ...
    spec,A,Bslow,rawState); %#ok<NASGU>
spec = wnl_spec(0,0,0,2,'reject_noncompacting_completion');
solution = wnl_solve_forced(model,spec,forcing,{},opts);
assert(solution.algebraicCompletion.attempted);
assert(solution.algebraicCompletion.available);
assert(~solution.algebraicCompletion.accepted);
assert(~solution.algebraicCompletion.acceptance.residualImproved);
assert(~solution.algebraicCompletion.acceptance. ...
    substantialCompaction);
assert(contains(solution.algebraicCompletion.acceptance.reason, ...
    'without substantial'));

A = sparse(5,5);
A(1,:) = [1,0,1,1,0];
A(2,:) = [0,1,0,1,1];
Bslow = sparse(diag([1,1,0,0,0]));
inflation = 1.0e12;
rawState = [3;4;inflation;1-inflation;inflation];
forcing = [4.001;5;0;0;0];
model.block = @(spec) forced_seed_test_block( ...
    spec,A,Bslow,rawState); %#ok<NASGU>
spec = wnl_spec(0,0,0,5,'accept_compacting_completion');
solution = wnl_solve_forced(model,spec,forcing,{},opts);
assert(solution.algebraicCompletion.attempted);
assert(solution.algebraicCompletion.available);
assert(solution.algebraicCompletion.accepted);
assert(solution.algebraicCompletion.acceptance. ...
    substantialCompaction);
end

function block = forced_seed_test_block(spec,A,Bslow,rawState) %#ok<INUSD>
block = struct('A',A,'Bslow',Bslow, ...
    'solve',@(forcing,blockSpec,options) rawState); %#ok<INUSD>
end

function test_temporal_block_gmres_forced()
% A coupled three-harmonic forced system should be solved by the scaled
% temporal-block path before the generic minimum-norm fallback is needed.
diagonalBlock = [4,1;1,3];
couplingBlock = -0.35*eye(2);
temporalCoupling = diag(ones(2,1),1)+diag(ones(2,1),-1);
A = sparse(kron(eye(3),diagonalBlock)+ ...
    kron(temporalCoupling,couplingBlock));
exact = (1:6).'+1i*(6:-1:1).';
forcing = A*exact;
model.block = @(spec) struct( ...
    'A',A,'Bslow',speye(6)); %#ok<NASGU>
spec = wnl_spec(0,0,1,2,'scaled_forced_block_gmres_test');
opts = struct('verbose',false,'forcedUseBlockGmres',true, ...
    'forcedSolveResidualTolerance',1.0e-10, ...
    'forcedSolveRefinementTolerance',1.0e-12, ...
    'forcedGmresRestart',10,'forcedGmresMaxCycles',10, ...
    'forcedSolveMaxRestarts',0);
solution = wnl_solve_forced(model,spec,forcing,{},opts);
assert(solution.gmresDiagnostics.attempted);
assert(solution.gmresDiagnostics.accepted);
assert(solution.gmresDiagnostics.improvedSeed);
assert(solution.gmresDiagnostics.passedPhysicalGate);
assert(~solution.gmresDiagnostics.minimumNormFallbackUsed);
assert(solution.gmresDiagnostics.relativeResidual < 1.0e-10);
assert(solution.valid);
assert(norm(solution.vector-exact)/norm(exact) < 1.0e-9);
assert(all(solution.gmresDiagnostics.blockRegularization > 0));
end

function test_descriptor_residual_rejects_algebraic_inflation()
% A pressure/gauge component in a zero-mass column must not make an
% unconverged physical eigenvector appear neutral.
A = sparse(diag([1,0]));
Bslow = sparse(diag([1,0]));
[directReference,referenceDetails] = wnl_descriptor_residual( ...
    A,Bslow,[1;0],'direct');
[directInflated,inflatedDetails] = wnl_descriptor_residual( ...
    A,Bslow,[1;1.0e12],'direct');
assert(abs(directInflated-directReference) < 10*eps);
assert(inflatedDetails.fullToPhysicalNormRatio > 1.0e11);
assert(referenceDetails.numberOfPhysicalEntries == 1);
legacyResidual = norm(A*[1;1.0e12]) / ...
    (max(1,sqrt(norm(A,1)*norm(A,inf)))*norm([1;1.0e12]));
assert(legacyResidual < 1.0e-10*directInflated);

[adjointReference,~] = wnl_descriptor_residual( ...
    A,Bslow,[1;0],'adjoint');
[adjointInflated,adjointDetails] = wnl_descriptor_residual( ...
    A,Bslow,[1;1.0e12],'adjoint');
assert(abs(adjointInflated-adjointReference) < 10*eps);
assert(adjointDetails.fullToPhysicalNormRatio > 1.0e11);
end

function test_output_save_path_resolution()
temporaryRoot = tempname;
[created,message] = mkdir(temporaryRoot);
assert(created,message);
cleanup = onCleanup(@() remove_test_directory(temporaryRoot)); %#ok<NASGU>
record.answer = 42;
[savedRecord,savedFile,information] = vi_save_output_record( ...
    record,fullfile('nested','wnl_test_record'),temporaryRoot);
expectedFile = fullfile(temporaryRoot,'nested','wnl_test_record.mat');
assert(strcmp(savedFile,expectedFile));
assert(isfile(savedFile));
assert(~information.usedFallback);
assert(strcmp(savedRecord.save.savedFile,savedFile));
loaded = load(savedFile,'output');
assert(loaded.output.answer == 42);
assert(strcmp(loaded.output.save.savedFile,savedFile));
end

function remove_test_directory(directory)
if isfolder(directory)
    rmdir(directory,'s');
end
end

function test_reference_eigenvalue_gate()
mode.spec = struct('label','reference_consistency_test', ...
    'lambda',0.102);
mode.tracking.eigenpairRefinement = struct('attempted',true, ...
    'initialLambda',0.1);
opts = wnl_options(struct( ...
    'maximumReducedFullEigenvalueRelativeMismatch',0.05, ...
    'maximumReducedFullEigenvalueAbsoluteMismatch',1.0e-3));
diagnostic = wnl_assert_reference_eigenvalue_consistent(mode,opts);
assert(diagnostic.consistent);
mode.spec.lambda = 0.12;
didReject = false;
try
    wnl_assert_reference_eigenvalue_consistent(mode,opts);
catch referenceError
    didReject = strcmp(referenceError.identifier, ...
        ['wnl_assert_reference_eigenvalue_consistent:', ...
         'MismatchTooLarge']);
end
assert(didReject);

cachedMode = mode;
cachedMode.tracking = struct('eigenpairRefinement', ...
    struct('attempted',false));
cachedMode.spec.reducedReferenceLambda = 0.1;
didRejectCached = false;
try
    wnl_assert_reference_eigenvalue_consistent(cachedMode,opts);
catch referenceError
    didRejectCached = strcmp(referenceError.identifier, ...
        ['wnl_assert_reference_eigenvalue_consistent:', ...
         'MismatchTooLarge']);
end
assert(didRejectCached);
end

function test_two_mode_shared_forced_fields()
config.omega = 1;
config.N = 0;
config.ndof = 1;
config.mass = 1;
config.linearFourier = @shared_field_test_linear;
config.quadraticLocal = @shared_field_test_quadratic;
config.cubicLocal = @shared_field_test_cubic;
model = wnl_fourier_model(config);

spec1 = model.makeSpec(1,0,'shared_mode_1');
spec1.direct = 1;
spec1.left = 1;
spec2 = model.makeSpec(3,0,'shared_mode_2');
spec2.direct = 1;
spec2.left = 1;
opts = struct('verbose',false);
result = wnl_analyze_mode_set(model,{spec1;spec2},opts);
assert(norm(result.g-[9,18;18,9],'fro') < 1.0e-12);
assert(result.optimization.twoModeSharedForcedFields);
assert(result.optimization.reusedForcedSolveCount == 4);
assert(result.optimization.actualForcedSolveCount == 6);
assert(result.optimization.conjugateDifferenceReused);
assert(result.cross{1,2}.reusedForcedFields.qBbarB);
assert(all(structfun(@(value) value, ...
    result.cross{2,1}.reusedForcedFields)));

baseline12 = wnl_cross_coefficient(model,result.modes{1}, ...
    result.modes{2},result.conjugateModes{2}, ...
    result.neutralModes,opts);
baseline21 = wnl_cross_coefficient(model,result.modes{2}, ...
    result.modes{1},result.conjugateModes{1}, ...
    result.neutralModes,opts);
assert(abs(result.g(1,2)-baseline12.g) < 1.0e-12);
assert(abs(result.g(2,1)-baseline21.g) < 1.0e-12);

duplicateOpts = struct('verbose',false, ...
    'reuseTwoModeForcedFields',false);
duplicate = wnl_analyze_mode_set(model,{spec1;spec2},duplicateOpts);
assert(norm(duplicate.g-result.g,'fro') < 1.0e-12);
assert(~duplicate.optimization.twoModeSharedForcedFields);
assert(duplicate.optimization.reusedForcedSolveCount == 0);
assert(duplicate.optimization.actualForcedSolveCount == 10);
end

function value = shared_field_test_linear(spec,k)
if k ~= 0
    value = sparse(1,1);
elseif ismember(abs(spec.m),[1,3])
    value = sparse(1,1);
else
    value = sparse(-1);
end
end

function value = shared_field_test_quadratic(a,b,~,~,~)
value = a.*b;
end

function value = shared_field_test_cubic(a,b,c,~,~,~,~)
value = a.*b.*c;
end

function test_cross_forced_fail_fast()
config.omega = 1;
config.N = 0;
config.ndof = 1;
config.mass = 1;
config.linearFourier = @fail_fast_test_linear;
config.quadraticLocal = @shared_field_test_quadratic;
config.cubicLocal = @shared_field_test_cubic;
model = wnl_fourier_model(config);
spec1 = model.makeSpec(1,0,'fail_fast_mode_1');
spec1.direct = 1;
spec1.left = 1;
spec2 = model.makeSpec(3,0,'fail_fast_mode_2');
spec2.direct = 1;
spec2.left = 1;
opts = struct('verbose',false,'forcedFailFast',true, ...
    'stopOnUnconvergedForcedSolve',true);
mode1 = wnl_compute_mode(model,spec1,opts);
mode2 = wnl_compute_mode(model,spec2,opts);
mode2Bar = wnl_conjugate_mode(model,mode2,opts);
neutralModes = wnl_unique_modes({mode1,mode2,mode2Bar});
result = wnl_cross_coefficient(model,mode1,mode2,mode2Bar, ...
    neutralModes,opts);
assert(~result.validCubicScaling && ~result.forcedSolvesValid);
assert(~isempty(result.qBbarB));
assert(isempty(result.qAbarB) && isempty(result.qAB));
assert(contains(result.message,'Fail-fast'));
end

function value = fail_fast_test_linear(spec,k)
if k ~= 0
    value = sparse(1,1);
elseif ismember(abs(spec.m),[1,3])
    value = sparse(1,1);
else
    % A=0 in every quadratic output block, so a nonzero forcing is
    % inconsistent and must fail before the next cross field is assembled.
    value = sparse(1,1);
end
end

function test_recovered_mode_cache()
temporaryRoot = tempname;
[created,message] = mkdir(temporaryRoot);
assert(created,message);
cleanup = onCleanup(@() remove_test_directory(temporaryRoot)); %#ok<NASGU>

input.numerics = struct('N',0,'Nr',3,'Ntheta',8);
input.boundary = struct('contactLine','free');
input.modes = struct('m',2,'radialIndex',1,'s',0.5, ...
    'label','cache_mode');
input.run.operatorFactory = 'cache_test_factory';
input.run.modeRecoveryOnly = false;
input.options.coefficientModeResidualTolerance = 1.0e-8;
parameters = struct('omegaStar',2,'R0',1,'C',0.1,'Bd',3, ...
    'At',0.2,'eta',0.3,'g_sgn',-1,'aAnalysis',0.7,'phase',0);
spec = struct('m',2,'s',0.5,'n',[-1,0],'ndof',2, ...
    'label','cache_mode','lambda',0.12,'radialIndex',1, ...
    'betaStar',2.5);
mode.spec = spec;
mode.spec.lambda = 0.121;
mode.vector = [1;2;3;4];
mode.left = [4;3;2;1];
mode.directResidual = 1.0e-11;
mode.leftResidual = 2.0e-11;
output.codeRelease = 'cache-test';
output.input = input;
output.parameters = parameters;
output.operatorMetadata.ndof = 2;
output.weaklyNonlinear.mode = mode;
cacheFile = fullfile(temporaryRoot,'cache.mat');
save(cacheFile,'output');

[updated,information] = vi_wnl_apply_recovery_cache( ...
    cacheFile,temporaryRoot,input,parameters,{spec},2);
assert(information.used && information.valid);
assert(isequal(updated{1}.direct,mode.vector));
assert(isequal(updated{1}.left,mode.left));
assert(abs(updated{1}.lambda-mode.spec.lambda) < eps);
assert(abs(updated{1}.reducedReferenceLambda-spec.lambda) < eps);

mode.directResidual = 2.0e-6;
output.weaklyNonlinear.mode = mode;
save(cacheFile,'output');
[resumed,warmInformation] = vi_wnl_apply_recovery_cache( ...
    cacheFile,temporaryRoot,input,parameters,{spec},2);
assert(~warmInformation.used && warmInformation.valid);
assert(warmInformation.warmStarted);
assert(isequal(resumed{1}.directLiftInitialVector,mode.vector));
assert(~isfield(resumed{1},'direct'));

changed = parameters;
changed.aAnalysis = 0.8;
[~,rejected] = vi_wnl_apply_recovery_cache( ...
    cacheFile,temporaryRoot,input,changed,{spec},2);
assert(~rejected.used && contains(rejected.reason,'aAnalysis'));
end

function test_runtime_profiles()
input.execution.profile = 'development';
input.execution.nonlinearTemporalOversampling = 2;
input.numerics.N = 11;
input.numerics.Nr = 12;
input.numerics.Ntheta = 32;
input.numerics.verticalGrid.type = 'multidomain';
input.numerics.verticalGrid.lowerBreaks = [];
input.numerics.verticalGrid.upperBreaks = [];
input.numerics.verticalGrid.pointsPerBoundaryLayer = 11;
input.numerics.verticalGrid.pointsInBulk = 15;
input.options.forcedUseBlockGmres = true;
input.options.directLiftMaxRestarts = 7;
input.options.modeTrackingMaxIterations = 20000;
input.options.eigenpairRefinementMaxSteps = 8;
input.options.eigenpairGmresMaxCycles = 80;
input.options.eigenpairLsqrFallbackMaxIterations = 3000;
input.options.forcedSolveMaxIterations = 10000;
input.options.forcedSolveMaxRestarts = 5;
[development,summary] = vi_wnl_apply_runtime_profile(input);
assert(development.numerics.N == 5);
assert(development.numerics.Nr == 9);
assert(development.numerics.Ntheta == 24);
assert(development.numerics.verticalGrid.pointsPerBoundaryLayer == 7);
assert(development.numerics.verticalGrid.pointsInBulk == 9);
assert(development.execution.nonlinearTemporalOversampling == 1);
assert(development.options.forcedUseBlockGmres);
assert(development.options.forcedGmresMaxCycles == 8);
assert(development.options.autoEscalateEigenpairRefinement);
assert(development.options.eigenpairGmresMaxCycles == 8);
assert(development.options.eigenpairLsqrFallbackMaxIterations == 0);
assert(development.options.eigenpairAlgebraicCompletionMaxIterations == 0);
assert(development.options.forcedAlgebraicCompletionMaxIterations == 0);
assert(strcmp(summary.profile,'development'));
assert(summary.autoEscalateEigenpairRefinement);

defaultProfileInput = rmfield(input,'execution');
[defaultProfile,defaultSummary] = ...
    vi_wnl_apply_runtime_profile(defaultProfileInput);
assert(strcmp(defaultProfile.execution.profile,'balanced'));
assert(strcmp(defaultSummary.profile,'balanced'));
assert(defaultProfile.options.eigenpairGmresMaxCycles == 15);
assert(defaultProfile.options.eigenpairLsqrFallbackMaxIterations == 750);
assert(defaultSummary.eigenpairCorrectionRegularization == 1.0e-12);
assert(defaultSummary.eigenpairCorrectionEquationTolerance == 1.0e-1);
assert(defaultSummary.eigenpairUseAlgebraicCompletion);
assert(defaultSummary.eigenpairAlgebraicCompletionMaxIterations == 500);
assert(defaultProfile.options.forcedAlgebraicCompletionMaxIterations == 500);
assert(defaultSummary.forcedUseAlgebraicCompletion);
assert(defaultSummary.forcedFailFast);
assert(defaultSummary.forcedUseCylinderSchur);
assert(defaultSummary.adaptiveDirectionalSteps);
assert(isequal(defaultSummary.quadraticStepMultipliers,[1,2,4,8]));
assert(isequal(defaultSummary.cubicStepMultipliers,[1,2,4]));
assert(defaultProfile.options.forcedUseBlockGmres);
assert(defaultProfile.options.forcedGmresMaxCycles == 30);

input.execution.profile = 'final';
[final,summary] = vi_wnl_apply_runtime_profile(input);
assert(final.numerics.N == 11);
assert(final.numerics.Nr == 12);
assert(final.numerics.Ntheta == 32);
assert(final.execution.nonlinearTemporalOversampling == 2);
assert(final.options.forcedUseBlockGmres);
assert(final.options.forcedGmresMaxCycles == 120);
assert(final.options.eigenpairAlgebraicCompletionMaxIterations == 1000);
assert(final.options.forcedAlgebraicCompletionMaxIterations == 1000);
assert(summary.isFinal);

badProfile = input;
badProfile.execution.profile = true;
profileTypeRejected = false;
try
    vi_wnl_apply_runtime_profile(badProfile);
catch profileError
    profileTypeRejected = strcmp(profileError.identifier, ...
        'vi_wnl_apply_runtime_profile:ProfileType');
end
assert(profileTypeRejected);
end

function test_shifted_floquet_block()
config.omega = 3.0;
config.N = 0;
config.ndof = 1;
config.mass = 2.0;
config.linearFourier = @(spec, k) shifted_scalar_operator(spec, k);
model = wnl_fourier_model(config);
spec = model.makeSpec(2, 0.25, 'shifted_block');
spec.lambda = 0.17-0.04i;
block = model.block(spec);
expected = 2*(spec.lambda+1i*config.omega*spec.s) - ...
    shifted_scalar_operator(spec, 0);
assert(abs(full(block.A)-expected) < 1.0e-13);
end

function value = shifted_scalar_operator(~, k)
if k == 0
    value = 0.31+0.08i;
else
    value = 0;
end
end

function test_shifted_spec_algebra()
config.omega = 1;
config.N = 1;
config.ndof = 1;
config.mass = 1;
config.linearFourier = @(~, ~) 0;
model = wnl_fourier_model(config);
a = model.makeSpec(2, 0.5, 'a');
a.lambda = 0.12+0.03i;
b = model.makeSpec(-1, 0, 'b');
b.lambda = -0.04+0.02i;
sumSpec = wnl_combine_spec(model, {a, b}, [1, 1], 'sum');
assert(sumSpec.m == 1 && abs(sumSpec.s-0.5) < eps);
assert(abs(sumSpec.lambda-(0.08+0.05i)) < 1.0e-14);

fieldA = wnl_make_field(a, ones(a.ndof, numel(a.n)));
barA = wnl_conjugate_field(model, fieldA, 'bar_a');
assert(abs(barA.spec.lambda-conj(a.lambda)) < 1.0e-14);
assert(~wnl_equivalent_spec(a, barA.spec));

aCopy = a;
assert(wnl_equivalent_spec(a, aCopy));
aCopy.lambda = a.lambda+1.0e-4;
assert(~wnl_equivalent_spec(a, aCopy));
end

function test_prescribed_lift_keeps_best_physical_residual()
% A restarted LSQR solve minimizes scaled, regularized equations, whereas
% WNL acceptance uses the unequilibrated full-operator residual. Verify that
% the returned lift is the best physical candidate over all attempts.
rng(19);
numberOfUnknowns = 12;
A = sparse(randn(numberOfUnknowns)+ ...
    diag(logspace(-5, 1, numberOfUnknowns)));
model = struct();
model.block = @(spec) struct('A', A, ...
    'Bslow', speye(numberOfUnknowns)); %#ok<NASGU>
spec = wnl_spec(2, 0, 0, numberOfUnknowns, ...
    'prescribed_lift_best_candidate_test');
spec.directSeed = complex(zeros(numberOfUnknowns, 1));
spec.directSeed(end) = 1;
spec.prescribedDofs = numberOfUnknowns;
opts = struct();
opts.verbose = false;
opts.left = ones(numberOfUnknowns, 1);
opts.modeTrackingRegularization = 1.0e-8;
opts.modeTrackingSolveTolerance = 1.0e-300;
opts.modeTrackingMaxIterations = 1;
opts.directLiftMaxRestarts = 2;
opts.modeResidualTolerance = 0;
mode = wnl_compute_mode(model, spec, opts);
diagnostics = mode.tracking.direct;
assert(diagnostics.numberOfAttempts == 3);
assert(diagnostics.restarts == 2);
assert(numel(diagnostics.attemptScaledResiduals) == 3);
[minimumResidual, minimumIndex] = ...
    min(diagnostics.attemptScaledResiduals);
assert(abs(diagnostics.scaledResidual-minimumResidual) <= ...
    20*eps*max(1, minimumResidual));
assert(diagnostics.selectedAttempt == minimumIndex);
assert(abs(mode.directResidual-minimumResidual) <= ...
    100*eps*max(1, minimumResidual));

warmSpec = spec;
warmSpec.directLiftInitialVector = mode.vector;
warmOpts = opts;
warmOpts.directLiftMaxRestarts = 0;
warmMode = wnl_compute_mode(model,warmSpec,warmOpts);
assert(warmMode.tracking.direct.warmStartUsed);
assert(isfinite(warmMode.tracking.direct.warmStartResidual));
assert(warmMode.directResidual <= ...
    warmMode.tracking.direct.warmStartResidual+100*eps);
end

function test_analysis_amplitude_resolution()
defaultOptions = wnl_options();
assert(defaultOptions.directLiftMaxRestarts == 3);
assert(defaultOptions.autoEscalateEigenpairRefinement);
assert(defaultOptions.eigenpairEscalationMaxSteps >= ...
    defaultOptions.eigenpairRefinementMaxSteps);
assert(defaultOptions.eigenpairEscalationMaximumResidualRatio > 1);
assert(defaultOptions.eigenpairCorrectionRegularization > 0);
assert(defaultOptions.eigenpairCorrectionEquationTolerance < 1);
assert(defaultOptions.eigenpairUseAlgebraicCompletion);
assert(defaultOptions.eigenpairAlgebraicCompletionMaxIterations > 0);
forcing.analysisAmplitude = 0.42;
linearSettings.accelerationOffsetFromCritical = 10;
[amplitude, details] = vi_resolve_analysis_amplitude( ...
    forcing, linearSettings, 0.3);
assert(abs(amplitude-0.42) < eps);
assert(abs(details.detuningFromCritical-0.12) < 10*eps);
assert(details.ignoredLegacyOffset && ~details.usedLegacyOffset);

forcing.analysisAmplitude = [];
[legacyAmplitude, legacyDetails] = ...
    vi_resolve_analysis_amplitude(forcing, ...
    struct('accelerationOffsetFromCritical', 0.05), 0.3);
assert(abs(legacyAmplitude-0.35) < eps);
assert(legacyDetails.usedLegacyOffset);

didReject = false;
try
    vi_resolve_analysis_amplitude( ...
        struct('analysisAmplitude', -0.1), struct(), 0.3);
catch amplitudeError
    didReject = strcmp(amplitudeError.identifier, ...
        'vi_resolve_analysis_amplitude:NegativeAmplitude');
end
assert(didReject);
end

function test_two_mode_initial_conditions_and_reconstruction()
modes(1).label = 'mode_one';
modes(2).label = 'mode_two';
settings.amplitudesOverH = [1e-3; 2e-4];
settings.phases = [0; pi/2];
initial = vi_mode_initial_conditions(modes, settings);
assert(initial.numberOfModes == 2);
assert(abs(initial.complexAmplitudesOverH(1)-1e-3) < 1e-15);
assert(abs(initial.complexAmplitudesOverH(2)-2e-4i) < 1e-15);

settings.amplitudesOverH(2) = 0;
zeroSeeded = vi_mode_initial_conditions(modes, settings);
assert(isequal(zeroSeeded.zeroSeededModes, 2));
badSettings = settings;
badSettings.phases = 0;
didReject = false;
try
    vi_mode_initial_conditions(modes, badSettings);
catch initialError
    didReject = strcmp(initialError.identifier, ...
        'vi_mode_initial_conditions:Size');
end
assert(didReject);

r = [0; 0.5; 1];
time = linspace(0, 2*pi, 17).';
spec1 = wnl_spec(1, 0, 0, numel(r), 'mode_one');
spec2 = wnl_spec(2, 0, 0, numel(r), 'mode_two');
mode1 = struct('spec', spec1, 'field', ...
    wnl_make_field(spec1, besselj(1, r)));
mode2 = struct('spec', spec2, 'field', ...
    wnl_make_field(spec2, 0.5*besselj(2, r)));
self1.qAA.field = wnl_make_field( ...
    wnl_spec(2, 0, 0, numel(r), 'self1_AA'), 0.2*ones(size(r)));
self1.qAbarA.field = wnl_make_field( ...
    wnl_spec(0, 0, 0, numel(r), 'self1_mean'), 0.1*ones(size(r)));
self2.qAA.field = wnl_make_field( ...
    wnl_spec(4, 0, 0, numel(r), 'self2_AA'), 0.15*ones(size(r)));
self2.qAbarA.field = wnl_make_field( ...
    wnl_spec(0, 0, 0, numel(r), 'self2_mean'), 0.05*ones(size(r)));
cross12.qAB.field = wnl_make_field( ...
    wnl_spec(3, 0, 0, numel(r), 'cross_sum'), 0.07*ones(size(r)));
cross12.qAbarB.field = wnl_make_field( ...
    wnl_spec(-1, 0, 0, numel(r), 'cross_difference'), ...
    0.03*ones(size(r)));
cross = cell(2, 2);
cross{1, 2} = cross12;
wnlResult = struct('modes', {{mode1; mode2}}, ...
    'self', {{self1; self2}}, 'cross', {cross});

linearResult.timeStar = time;
linearResult.forcingPeriods = time/(2*pi);
linearResult.displacementOverH = 1e-3*ones(size(time));
linearResult.floquetOscillation = ones(size(time));
linearResult.betaStar = 1;
metadata.layout.zeta = 1:numel(r);
metadata.discretization.r = r;
parameters.omegaStar = 1;
parameters.R0 = 1;
plotSettings.plotInterfaceDynamics = false;
plotSettings.includeSlavedHarmonics = true;
amplitudes = [1e-3*ones(size(time)), 2e-4*ones(size(time))];
result = vi_compare_interface_dynamics_modes(linearResult, ...
    amplitudes, wnlResult, metadata, parameters, plotSettings, 1);
assert(result.numberOfModes == 2);
assert(numel(result.secondOrderLabels) == 6);
assert(isequal(size(result.modeAmplitudesOverH), size(amplitudes)));
assert(norm(result.wnlTotalField-(result.wnlPrimaryField+ ...
    result.wnlSecondOrderField), 'fro') < 1e-14);
twoModeLinear = vi_compare_interface_dynamics_modes(linearResult, ...
    amplitudes, wnlResult, metadata, parameters, plotSettings, 1, ...
    amplitudes);
assert(twoModeLinear.linearReferenceIncludesBothModes);
assert(isequal(twoModeLinear.linearModeAmplitudesOverH,amplitudes));
assert(isequal(twoModeLinear.realModeAmplitudesLinear,real(amplitudes)));
assert(isequal(twoModeLinear.realModeAmplitudesWnl,real(amplitudes)));
assert(isequal(size(twoModeLinear.linearModeProbeSignals), ...
    size(amplitudes)));
assert(isequal(size(twoModeLinear.wnlPrimaryModeProbeSignals), ...
    size(amplitudes)));
assert(norm(sum(twoModeLinear.wnlPrimaryModeProbeSignals,2)- ...
    twoModeLinear.wnlPrimaryProbeSignal) < 1e-14);
assert(norm(twoModeLinear. ...
    realModeAmplitudeDifferenceLinearMinusWnl,'fro') == 0);
assert(isequal(size(twoModeLinear.linearField), ...
    size(twoModeLinear.wnlTotalField)));
end

function test_comparison_time_window()
periods = linspace(0, 20, 4001);
[indices, actualEnd] = vi_comparison_time_indices(periods, 4.0);
assert(abs(actualEnd-4.0) < 1.0e-14);
assert(indices(1) == 1 && indices(end) == 801);
[fullIndices, fullEnd] = vi_comparison_time_indices(periods, []);
assert(numel(fullIndices) == numel(periods));
assert(fullEnd == 20);

didReject = false;
try
    vi_comparison_time_indices(periods, 21);
catch timeError
    didReject = strcmp(timeError.identifier, ...
        'vi_comparison_time_indices:EndOutsideRecord');
end
assert(didReject);
end

function test_nonfinite_mode_gate()
mode = struct();
mode.spec = struct('label', 'nonfinite_gate_test');
mode.directResidual = 0;
mode.leftResidual = NaN;
opts = wnl_options(struct('stopOnUnconvergedMode', false));
didReject = false;
try
    wnl_assert_mode_converged(mode, opts);
catch gateError
    didReject = strcmp(gateError.identifier, ...
        'wnl_assert_mode_converged:NonfiniteResidual');
end
assert(didReject);

% A conjugate adjoint that has already passed the linear gate is not part
% of the strict nonresonant cubic projection. Verify that the component
% selector gates the direct field without silently accepting a bad direct
% residual.
componentOpts = wnl_options(struct( ...
    'modeResidualTolerance',1.0e-8, ...
    'checkAdjointModeResidual',false));
mode.directResidual = 1.0e-9;
wnl_assert_mode_converged(mode,componentOpts);
mode.directResidual = 1.0e-7;
didReject = false;
try
    wnl_assert_mode_converged(mode,componentOpts);
catch gateError
    didReject = strcmp(gateError.identifier, ...
        'wnl_assert_mode_converged:ResidualTooLarge');
end
assert(didReject);
end

function test_interface_dynamics_reconstruction()
r = [0; 0.5; 1.0];
time = linspace(0, 2*pi, 17).';
amplitude = 1.0e-2*ones(size(time));
shape = besselj(1, r);
shape = shape/max(abs(shape));

modeSpec = wnl_spec(1, 0, 0, numel(r), 'interface_test');
mode = struct();
mode.spec = modeSpec;
mode.field = wnl_make_field(modeSpec, shape);

qAASpec = wnl_spec(2, 0, 0, numel(r), 'interface_test_AA');
qMeanSpec = wnl_spec(0, 0, 0, numel(r), 'interface_test_mean');
self = struct();
self.qAA = struct('field', wnl_make_field( ...
    qAASpec, 0.2*ones(size(r))));
self.qAbarA = struct('field', wnl_make_field( ...
    qMeanSpec, 0.1*ones(size(r))));
wnlResult = struct('mode', mode, 'self', self);

linearResult = struct();
linearResult.timeStar = time;
linearResult.forcingPeriods = time/(2*pi);
linearResult.displacementOverH = amplitude;
linearResult.floquetOscillation = ones(size(time));
linearResult.betaStar = 1.0;
metadata = struct();
metadata.layout.zeta = 1:numel(r);
metadata.discretization.r = r;
parameters = struct('omegaStar', 1.0, 'R0', 1.0);

settings = struct();
settings.plotInterfaceDynamics = false;
settings.includeSlavedHarmonics = false;
primaryOnly = vi_compare_interface_dynamics(linearResult, amplitude, ...
    wnlResult, metadata, parameters, settings, 1.0);
assert(primaryOnly.fieldRmse < 1.0e-14);
assert(abs(primaryOnly.normalizedPrimaryCarrierOverlap-1) < 1.0e-14);

settings.includeSlavedHarmonics = true;
withSlaved = vi_compare_interface_dynamics(linearResult, amplitude, ...
    wnlResult, metadata, parameters, settings, 1.0);
assert(withSlaved.includedSecondAzimuthalHarmonic);
assert(withSlaved.includedMeanInterfaceCorrection);
assert(withSlaved.secondOrderRelativeL2 > 0);
assert(norm(withSlaved.wnlTotalField - ...
    (withSlaved.wnlPrimaryField + ...
    withSlaved.wnlSecondOrderField), 'fro') < 1.0e-14);
end

function test_temporal_block_gmres_adjoint()
% The two uncoupled temporal blocks are individually singular because of a
% constraint direction, while their coupled matrix has one simple physical
% null vector. This reproduces the flag-2 failure of an exact block inverse
% and exercises the regularized bordered block-Jacobi path.
temporalLaplacian = [1, -1; -1, 1];
spatialConstraint = diag([-1, 0]);
A = sparse(kron(temporalLaplacian, eye(2)) + ...
    kron(eye(2), spatialConstraint));
model = struct();
model.block = @(spec) struct('A', A, 'Bslow', speye(4)); %#ok<NASGU>
spec = wnl_spec(1, 0.5, 0, 2, 'block_gmres_adjoint_test');
spec.directSeed = [0;1;0;1];
opts = struct();
opts.verbose = false;
opts.modeTrackingRegularization = 0;
opts.adjointUseBlockGmres = true;
mode = wnl_compute_mode(model, spec, opts);
assert(mode.directResidual < 1.0e-12);
assert(mode.leftResidual < 1.0e-12);
assert(abs(mode.normalization-1) < 1.0e-12);
assert(mode.tracking.left.gmres.attempted);
assert(mode.tracking.left.gmres.succeeded);
assert(mode.tracking.left.candidateValid);
assert(all(mode.tracking.left.gmres.blockRegularization > 0));
assert(strcmp(mode.tracking.left.gmres.preconditioner, ...
    'regularized bordered block-Jacobi'));
assert(isfinite(mode.tracking.left.gmres.borderSchurComplement));
assert(abs(mode.tracking.left.gmres.borderSchurComplement) > 0);
assert(mode.tracking.left.gmres.constraintResidual < 1.0e-12);
end

function test_temporal_block_gmres_eigenpair()
% A three-harmonic affine pencil has its tracked central eigenvalue at 0.1
% but starts from lambda=0.2. The bordered direct correction should use the
% temporal-block preconditioner and recover the exact exponent in one step.
model = struct();
model.block = @affine_eigenpair_test_block;
spec = wnl_spec(1,0,1,1,'block_gmres_eigenpair_test');
spec.lambda = 0.2;
spec.directSeed = [0;1;0];
spec.prescribedDofs = 1;
opts = struct();
opts.verbose = false;
opts.modeTrackingRegularization = 0;
opts.refineOperatingPointEigenpair = true;
opts.eigenpairUseBlockGmres = true;
opts.eigenpairGmresAcceptanceTolerance = 1.0e-5;
opts.eigenpairRefinementMaxSteps = 3;
opts.eigenpairRefinementTolerance = 1.0e-12;
opts.eigenpairMaximumSeedResidualRatio = 1.0e12;
opts.stopEigenpairOnReferenceMismatch = false;
opts.autoEscalateEigenpairRefinement = false;
mode = wnl_compute_mode(model,spec,opts);
effectiveOpts = wnl_options(opts);
assert(abs(wnl_spec_lambda(mode.spec)-0.1) < 1.0e-10);
assert(mode.directResidual < 1.0e-10);
assert(mode.tracking.eigenpairRefinement.steps == 1);
assert(strcmp(mode.tracking.eigenpairRefinement.solverHistory{1}, ...
    'block-GMRES'));
assert(mode.tracking.eigenpairRefinement.linearSolveResidualHistory(1) < ...
    opts.eigenpairGmresAcceptanceTolerance);
assert(mode.tracking.eigenpairRefinement.minimumDescriptorRetention >= ...
    effectiveOpts.eigenpairMinimumDescriptorRetention);
assert(mode.tracking.eigenpairRefinement.minimumDescriptorOverlap >= ...
    effectiveOpts.eigenpairMinimumDescriptorOverlap);
assert(max(mode.tracking.eigenpairRefinement. ...
    phaseConstraintResidualHistory) <= ...
    effectiveOpts.eigenpairCorrectionConstraintTolerance);
assert(norm(mode.block.Bslow*mode.vector) > 0);
assert(strcmp(mode.tracking.eigenpairRefinement.phaseConstraintSource, ...
    'prescribed interface displacement'));
end

function test_physical_coordinate_lsqr_eigenpair()
% Exercise the pressure-safe LSQR correction path without block GMRES.
model = struct();
model.block = @affine_eigenpair_test_block;
spec = wnl_spec(1,0,1,1,'physical_lsqr_eigenpair_test');
spec.lambda = 0.2;
spec.directSeed = [0;1;0];
spec.prescribedDofs = 1;
opts = struct();
opts.verbose = false;
opts.modeTrackingRegularization = 0;
opts.refineOperatingPointEigenpair = true;
opts.eigenpairUseBlockGmres = false;
opts.eigenpairLsqrFallbackMaxIterations = 100;
opts.eigenpairLsqrFallbackMaxRestarts = 0;
opts.eigenpairRefinementMaxSteps = 2;
opts.eigenpairRefinementTolerance = 1.0e-12;
opts.eigenpairMaximumSeedResidualRatio = 1.0e12;
opts.stopEigenpairOnReferenceMismatch = false;
opts.autoEscalateEigenpairRefinement = false;
mode = wnl_compute_mode(model,spec,opts);
assert(abs(wnl_spec_lambda(mode.spec)-0.1) < 1.0e-10);
assert(mode.directResidual < 1.0e-10);
assert(mode.tracking.eigenpairRefinement.steps == 1);
assert(strcmp(mode.tracking.eigenpairRefinement.solverHistory{1}, ...
    'regularized LSQR fallback'));
assert(mode.tracking.eigenpairRefinement.linearSolveResidualHistory(1) <= ...
    1.0e-1);
assert(mode.directResidualDetails.fullToPhysicalNormRatio < 2);
end

function block = affine_eigenpair_test_block(spec)
lambda = wnl_spec_lambda(spec);
block.A = sparse(diag([lambda+1;lambda-0.1;lambda+2]));
% The descriptor is deliberately singular, as it is in the primitive
% velocity-pressure cylinder system. The tracked physical component is the
% central entry; algebraic first/third components have no slow derivative.
block.Bslow = sparse(diag([0,1,0]));
end

function test_inaccurate_seed_skips_eigenpair()
model = struct();
model.block = @affine_eigenpair_test_block;
spec = wnl_spec(1,0,1,1,'inaccurate_eigenpair_seed_test');
spec.lambda = 0.2;
spec.directSeed = [0;1;0];
opts = struct();
opts.verbose = false;
opts.modeTrackingRegularization = 0;
opts.refineOperatingPointEigenpair = true;
opts.modeResidualTolerance = 1.0e-12;
opts.eigenpairMaximumSeedResidualRatio = 1.0;
opts.autoEscalateEigenpairRefinement = false;
mode = wnl_compute_mode(model,spec,opts);
refinement = mode.tracking.eigenpairRefinement;
assert(refinement.skippedForInaccurateSeed);
assert(refinement.steps == 0);
assert(abs(wnl_spec_lambda(mode.spec)-0.2) < 1.0e-14);
assert(contains(refinement.stopReason,'safe physical Newton seed'));
end

function test_bordered_adjoint()
% A strongly nonnormal matrix with one simple physical null vector. The
% bordered adjoint solve must recover the left null vector and exact
% descriptor normalization without a penalty parameter.
A = sparse([1, 100, 0; 0, 1, 0; 0, 0, 0]);
model = struct();
model.block = @(spec) struct('A', A, 'Bslow', speye(3)); %#ok<NASGU>
spec = wnl_spec(2, 0, 0, 3, 'bordered_adjoint_test');
spec.directSeed = [0;0;1];
opts = struct();
opts.verbose = false;
opts.modeTrackingRegularization = 0;
mode = wnl_compute_mode(model, spec, opts);
assert(mode.directResidual < 1.0e-12);
assert(mode.leftResidual < 1.0e-12);
assert(abs(mode.normalization-1) < 1.0e-12);
assert(strcmp(mode.tracking.left.method, ...
    'exactly normalized bordered adjoint'));
end

function test_floquet_forcing_phase()
phase = 0.37;
B = vi_floquet_acceleration_matrix(5, phase);
assert(abs(B(2, 1)+exp(1i*phase)) < 1.0e-14);
assert(abs(B(1, 2)+exp(-1i*phase)) < 1.0e-14);
assert(norm(B-B', 'fro') < 1.0e-14);
end

function test_reduced_neutral_mode_consistency()
omegaStar = 8.0;
R0 = 35/22;
C = 9.82e-5;
Bd = 65.78;
At = 0.9976;
eta = 1.81e-2;
N = 3;
phase = 0.41;
[acceleration, thresholdVector, thresholdDiagnostics] = ...
    faradayFloquet_RT_boundary_GenEIG_cylindrical( ...
    omegaStar, R0, 2, 1, C, Bd, At, eta, N, 'SH', 1, phase);
[fixedVector, fixedDiagnostics] = ...
    vi_reduced_cylinder_mode_at_exponent(acceleration, omegaStar, ...
    R0, 2, 1, C, Bd, At, eta, N, 'SH', 1, ...
    0.5i*omegaStar, phase);
overlap = abs(thresholdVector'*fixedVector) / ...
    (norm(thresholdVector)*norm(fixedVector));
assert(thresholdDiagnostics.vectorResidual < 1.0e-10);
assert(fixedDiagnostics.vectorResidual < 1.0e-10);
assert(overlap > 1-1.0e-8);
end

function test_seeded_mode_tracking()
% Twenty globally smaller directions have zero interface projection. The
% constrained solve must still reconstruct component 21 and its adjoint.
model = struct();
model.block = @(spec) struct('A', ...
    sparse(diag([zeros(1, 20), 1.0e-6])), ...
    'Bslow', speye(21)); %#ok<NASGU>
spec = wnl_spec(1, 0, 0, 21, 'tracked_test');
spec.directSeed = [zeros(20, 1); 1];
opts = struct();
opts.verbose = false;
mode = wnl_compute_mode(model, spec, opts);
assert(abs(mode.vector(21)) > 1-1.0e-12);
assert(abs(mode.left(21)) > 1-1.0e-12);
assert(mode.tracking.used);
assert(strcmp(mode.tracking.direct.method, ...
    'interface-constrained full residual'));
end

function test_periodic_nonlinear_shift()
model = struct();
model.makeSpec = @(m, s, label) wnl_spec(m, s, 2, 1, label);
spec = model.makeSpec(1, 0, 'input');
coefficient = complex(zeros(1, numel(spec.n)));
coefficient(spec.n == 0) = 2;
field = wnl_make_field(spec, coefficient);
out = model.makeSpec(2, 0, 'shifted_square');
product = wnl_quadratic_convolution(field, field, out, ...
    @(a, b, sa, sb, so, k, nuA, nuB) ...
    shifted_product(a, b, k, nuA, nuB), [-1, 0, 1]); %#ok<NASGU>
assert(abs(product(out.n == -1)-4) < 1.0e-14);
assert(abs(product(out.n == 0)-8) < 1.0e-14);
assert(abs(product(out.n == 1)-12) < 1.0e-14);
end

function value = shifted_product(a, b, k, nuA, nuB) %#ok<INUSD>
value = (k+2)*a*b;
end

function test_subharmonic_conjugation()
model = struct();
model.makeSpec = @(m, s, label) wnl_spec(m, s, 2, 1, label);
spec = model.makeSpec(3, 0.5, 'subharmonic');
coeff = reshape(1:numel(spec.n), 1, []);
field = wnl_make_field(spec, coeff);
twice = wnl_conjugate_field(model, ...
    wnl_conjugate_field(model, field));
assert(wnl_equivalent_spec(field.spec, twice.spec));
assert(norm(field.coeff - twice.coeff) == 0);
end

function test_subharmonic_product()
model = struct();
model.makeSpec = @(m, s, label) wnl_spec(m, s, 2, 1, label);
spec = model.makeSpec(1, 0.5, 'subharmonic');
coeff = complex(zeros(1, numel(spec.n)));
coeff(spec.n == 0) = 2.0;
field = wnl_make_field(spec, coeff);
out = model.makeSpec(2, 0, 'square');
product = wnl_quadratic_convolution(field, field, out, ...
    @(a, b, sa, sb, so) a * b); %#ok<INUSD>
assert(abs(product(out.n == 1) - 4.0) < 1.0e-14);
assert(nnz(product) == 1);
end

function test_bessel_derivative_roots()
roots0 = bessel_derivative_root(0, 2);
roots1 = bessel_derivative_root(1, 2);
assert(max(abs(roots0 - [3.83170597020751, 7.01558666981562])) < 1.0e-9);
assert(max(abs(roots1 - [1.84118378134066, 5.33144277352503])) < 1.0e-9);
end

function test_multidomain_chebyshev_grid()
breaks = [-1, -0.82, -0.17, 0];
points = [6, 8, 6];
grid = vi_multidomain_chebyshev_grid(breaks, points);
assert(grid.numberOfElements == 3);
assert(grid.numberOfPoints == sum(points));
assert(isequal(size(grid.interfacePairs), [2, 2]));
assert(max(abs(grid.z(grid.interfacePairs(:, 1))- ...
    grid.z(grid.interfacePairs(:, 2)))) < 10*eps);
assert(max(abs(grid.z(grid.interfacePairs(:, 1))- ...
    breaks(2:end-1).')) < 10*eps);

% A degree-five polynomial is represented exactly in every element.
values = grid.z.^5-0.3*grid.z.^3+0.2*grid.z;
firstExact = 5*grid.z.^4-0.9*grid.z.^2+0.2;
secondExact = 20*grid.z.^3-1.8*grid.z;
assert(norm(grid.D*values-firstExact, inf) < 2.0e-10);
assert(norm(grid.D2*values-secondExact, inf) < 2.0e-9);

% Preserve the legacy one-domain path and exercise the automatic
% Stokes-layer partition independently of the full operator smoke test.
p.At = 0.9976;
p.C = 9.82e-5;
p.eta = 1.81e-2;
p.omegaStar = 8.0;
p.numerics.NzLower = 9;
p.numerics.NzUpper = 8;
p.numerics.verticalGrid.type = 'single';
singleLower = vi_cylinder_vertical_grid(p, 'lower');
assert(singleLower.numberOfElements == 1);
assert(singleLower.numberOfPoints == p.numerics.NzLower);
p.numerics.verticalGrid.type = 'multidomain';
automaticLower = vi_cylinder_vertical_grid(p, 'lower');
automaticUpper = vi_cylinder_vertical_grid(p, 'upper');
assert(automaticLower.numberOfElements == 3);
assert(automaticUpper.numberOfElements == 3);
assert(abs(automaticLower.z(1)+1) < 10*eps && ...
    abs(automaticLower.z(end)) < 10*eps);
assert(abs(automaticUpper.z(1)) < 10*eps && ...
    abs(automaticUpper.z(end)-1) < 10*eps);
assert(automaticLower.pointsInsideWallStokesLayer >= 4);
assert(automaticUpper.pointsInsideInterfaceStokesLayer >= 4);
end

function test_cylinder_operator_smoke()
p = struct();
p.omegaStar = 8.0;
p.R0 = 35/22;
p.C = 9.82e-5;
p.Bd = 65.78;
p.At = 0.9976;
p.eta = 1.81e-2;
p.g_sgn = -1;
p.aCritical = 0.2;
p.phase = 0;
p.numerics.Nr = 6;
p.numerics.NzUpper = 7;
p.numerics.NzLower = 7;
p.numerics.verticalGrid.type = 'multidomain';
p.numerics.verticalGrid.lowerBreaks = [-1, -0.8, -0.2, 0];
p.numerics.verticalGrid.lowerPoints = [4, 5, 4];
p.numerics.verticalGrid.upperBreaks = [0, 0.2, 0.8, 1];
p.numerics.verticalGrid.upperPoints = [4, 5, 4];
p.modes.m = 2;
p.modes.radialIndex = 1;
p.modes.betaStar = bessel_derivative_root(2, 1)/p.R0;
p.numerics.radialGrid.type = 'besselEnriched';
p.numerics.radialGrid.maximumProductOrder = 2;
p.numerics.radialGrid.fallbackToChebyshev = false;
p.numerics.Ntheta = 12;
p.numerics.quadraticStep = 2.0e-4;
p.numerics.cubicStep = 2.0e-3;
p.boundary.contactLine = 'free';
[operators, metadata] = cylinder_wnl_operators(p);
specA = wnl_spec(2, 0.5, 0, metadata.ndof, 'A');
specAA = wnl_spec(4, 0, 1, metadata.ndof, 'AA');
B0 = operators.B0(specA);
L0 = operators.Lhat(specA, 0);
assert(isequal(size(B0), [metadata.ndof, metadata.ndof]));
assert(isequal(size(L0), [metadata.ndof, metadata.ndof]));
assert(nnz(operators.Lhat(specA, 2)) == 0);
assert(nnz(operators.P(specA, 0, 0, 1)) == 0);
assert(strcmp(metadata.radialGrid.typeUsed, 'besselEnriched'));
assert(metadata.verticalGrid.lower.numberOfElements == 3);
assert(metadata.verticalGrid.upper.numberOfElements == 3);
matchingRows = [vertical_matching_rows(metadata.layout.d, ...
    metadata.verticalGrid.lower.interfacePairs, metadata.layout.nr); ...
    vertical_matching_rows(metadata.layout.l, ...
    metadata.verticalGrid.upper.interfacePairs, metadata.layout.nr)];
assert(nnz(B0(matchingRows, :)) == 0);
assert(all(full(sum(abs(L0(matchingRows, :)), 2)) > 0));

rng(7);
x = randn(metadata.ndof, 1)+1i*randn(metadata.ndof, 1);
y = randn(metadata.ndof, 1)+1i*randn(metadata.ndof, 1);
x = x/max(norm(x), 1);
y = y/max(norm(y), 1);
cxy = operators.C(x, y, specA, specA, specAA, 0, 0.5, 0.5);
cyx = operators.C(y, x, specA, specA, specAA, 0, 0.5, 0.5);
assert(all(isfinite(cxy)));
assert(norm(cxy-cyx) < 1.0e-7*max(1, norm(cxy)));

% Exercise the field-level path used by wnl_self_coefficient. This reaches
% the complete interface traction, including cylindrical r-theta shear.
fieldSpec = wnl_spec(2, 0, 0, metadata.ndof, 'field_A');
fieldOut = wnl_spec(4, 0, 0, metadata.ndof, 'field_AA');
fieldX = wnl_make_field(fieldSpec, x);
fieldY = wnl_make_field(fieldSpec, y);
fieldQuadratic = operators.CField(fieldX, fieldY, fieldOut);
fieldQuadraticReverse = operators.CField(fieldY, fieldX, fieldOut);
assert(isequal(size(fieldQuadratic), ...
    [metadata.ndof, numel(fieldOut.n)]));
assert(all(isfinite(fieldQuadratic(:))));
assert(norm(fieldQuadratic-fieldQuadraticReverse,'fro') < ...
    1.0e-7*max(1,norm(fieldQuadratic,'fro')));
assert(norm(fieldQuadratic(matchingRows, :), 'fro') < 1.0e-13);
interfaceVelocityRows = zeros(0,1);
for radialIndex = 2:metadata.layout.nr-1
    interfaceNode = radialIndex + ...
        (metadata.layout.nzD-1)*metadata.layout.nr;
    interfaceVelocityRows = [interfaceVelocityRows; ...
        metadata.layout.d.ur(interfaceNode); ...
        metadata.layout.d.ut(interfaceNode); ...
        metadata.layout.d.w(interfaceNode)]; %#ok<AGROW>
end
assert(norm(fieldQuadratic(interfaceVelocityRows,:),'fro') == 0);

% The axis regularity equations exchange their stored ur/ut row locations
% under m -> -m. Verify both the direct operator conjugacy and the distinct
% equation-row map required by an algebraic adjoint.
config.omega = p.omegaStar;
config.N = 1;
config.ndof = metadata.ndof;
model = vi_wnl_model(config,operators);
assert(isfield(model,'conjugateAdjointRows'));
sourceSpec = model.makeSpec(2,0.5,'conjugacy_source');
sourceSpec.lambda = 0.13+0.02i;
sourceBlock = model.block(sourceSpec);
rng(29);
sourceState = randn(size(sourceBlock.A,2),1) + ...
    1i*randn(size(sourceBlock.A,2),1);
sourceStateField = wnl_make_field(sourceSpec,sourceState);
conjugateStateField = wnl_conjugate_field( ...
    model,sourceStateField,'conjugacy_state_bar');
conjugateBlock = model.block(conjugateStateField.spec);
sourceResidualField = wnl_make_field(sourceSpec, ...
    sourceBlock.A*sourceState);
expectedConjugateResidual = wnl_conjugate_field( ...
    model,sourceResidualField,'conjugacy_residual_bar');
expectedConjugateResidual = model.conjugateAdjointRows( ...
    expectedConjugateResidual,sourceSpec);
directConjugacyError = norm(conjugateBlock.A* ...
    wnl_field_vector(conjugateStateField) - ...
    wnl_field_vector(expectedConjugateResidual));
directConjugacyScale = max(1,norm(sourceBlock.A*sourceState));
assert(directConjugacyError < 1.0e-11*directConjugacyScale);

sourceLeft = randn(size(sourceBlock.A,1),1) + ...
    1i*randn(size(sourceBlock.A,1),1);
sourceLeftField = wnl_make_field(sourceSpec,sourceLeft);
conjugateLeftField = wnl_conjugate_field( ...
    model,sourceLeftField,'conjugacy_left_bar');
conjugateLeftField = model.conjugateAdjointRows( ...
    conjugateLeftField,sourceSpec);
sourceAdjointResidualField = wnl_make_field(sourceSpec, ...
    sourceBlock.A'*sourceLeft);
expectedConjugateAdjointResidual = wnl_conjugate_field( ...
    model,sourceAdjointResidualField,'conjugacy_adjoint_residual_bar');
adjointConjugacyError = norm(conjugateBlock.A'* ...
    wnl_field_vector(conjugateLeftField) - ...
    wnl_field_vector(expectedConjugateAdjointResidual));
adjointConjugacyScale = max(1,norm(sourceBlock.A'*sourceLeft));
assert(adjointConjugacyError < 1.0e-11*adjointConjugacyScale);
end

function test_bessel_enriched_radial_grid()
p.R0 = 35/22;
p.numerics.Nr = 9;
p.numerics.radialGrid.type = 'besselEnriched';
p.numerics.radialGrid.maximumProductOrder = 2;
p.numerics.radialGrid.maximumConditionNumber = 1.0e8;
p.numerics.radialGrid.fallbackToChebyshev = false;
p.modes.m = 2;
p.modes.radialIndex = 6;
roots = bessel_derivative_root(2, 6);
p.modes.betaStar = roots(6)/p.R0;
grid = vi_cylinder_radial_grid(p);
assert(strcmp(grid.typeUsed, 'besselEnriched'));
assert(grid.conditionNumber < ...
    p.numerics.radialGrid.maximumConditionNumber);

beta = p.modes.betaStar;
for order = 1:3
    argument = beta*grid.r;
    value = besselj(order, argument);
    exactFirst = 0.5*beta*(besselj(order-1, argument) - ...
        besselj(order+1, argument));
    exactSecond = 0.25*beta^2*(besselj(order-2, argument) - ...
        2*besselj(order, argument) + besselj(order+2, argument));
    assert(norm(grid.D*value-exactFirst, inf) < 1.0e-9);
    assert(norm(grid.D2*value-exactSecond, inf) < 1.0e-8);
end

% A ninth-order polynomial collocation matrix cannot resolve the sixth
% radial Bessel branch comparably. This verifies that the reduction comes
% from the enriched basis rather than merely from changing node locations.
[chebyshevR, chebyshevD] = vi_chebyshev_lobatto( ...
    p.numerics.Nr, [0, p.R0]);
argument = beta*chebyshevR;
value = besselj(2, argument);
exactFirst = 0.5*beta*(besselj(1, argument)-besselj(3, argument));
chebyshevError = norm(chebyshevD*value-exactFirst, inf);
enrichedError = norm(grid.D*value-exactFirst, inf);
assert(chebyshevError > 1.0e3*max(enrichedError, eps));

p.numerics.radialGrid.type = 'chebyshev';
legacy = vi_cylinder_radial_grid(p);
assert(strcmp(legacy.typeUsed, 'chebyshev'));
assert(norm(legacy.r-chebyshevR, inf) < 10*eps);
assert(norm(legacy.D-chebyshevD, inf) < 1.0e-12);
end

function rows = vertical_matching_rows(layer, interfacePairs, nr)
rows = zeros(0, 1);
for interfaceIndex = 1:size(interfacePairs, 1)
    leftIz = interfacePairs(interfaceIndex, 1);
    rightIz = interfacePairs(interfaceIndex, 2);
    for ir = 2:nr-1
        leftNode = ir+(leftIz-1)*nr;
        rightNode = ir+(rightIz-1)*nr;
        rows = [rows; layer.ur(leftNode); layer.ut(leftNode); ...
            layer.w(leftNode); layer.p(leftNode); ...
            layer.ur(rightNode); layer.ut(rightNode); ...
            layer.w(rightNode)]; %#ok<AGROW>
    end
end
end
