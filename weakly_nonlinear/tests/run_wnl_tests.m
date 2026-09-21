function run_wnl_tests()
%RUN_WNL_TESTS Basic regression tests for the weakly nonlinear module.

thisFile = mfilename('fullpath');
moduleRoot = fileparts(fileparts(thisFile));
repositoryRoot = fileparts(moduleRoot);
addpath(moduleRoot);
addpath(fullfile(moduleRoot, 'examples'));
addpath(repositoryRoot);

wnl_demo_stuart_landau();
test_real_modal_coefficient_combinatorics();
test_complex_target_real_source_cross();
test_self_conjugate_field_projection();
test_self_conjugate_forced_solution_projection();
test_saved_real_coefficient_reprojection();
test_two_real_mode_shared_fields();
test_modal_amplitude_scaling();
test_radial_mode_correction();
test_radial_mode_prolongation();
test_phase_sensitive_cubic_terms();
test_phase_sensitive_coefficient_projection();
test_cubic_mode_set_completeness_gate();
test_chebyshev_quadrature();
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
test_model_solver_only_forced_failure();
test_cylinder_temporal_schur_forced();
test_cylinder_temporal_schur_pressure_nullspace();
test_cylinder_temporal_schur_pressure_quotient();
test_cylinder_temporal_schur_adjoint_pressure_nullspace();
test_temporal_block_gmres_forced();
test_bordered_adjoint();
test_model_specific_adjoint();
test_temporal_block_gmres_adjoint();
test_temporal_block_gmres_eigenpair();
test_physical_coordinate_lsqr_eigenpair();
test_prescribed_interface_refinement_constraint();
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
test_mixed_coordinate_cross_reconstruction();
test_analysis_amplitude_resolution();
test_comparison_time_window();
test_small_amplitude_cubic_transient();
test_runtime_profiles();
test_programmatic_input_merge();
test_two_mode_shared_forced_fields();
test_cross_forced_fail_fast();
test_recovered_mode_cache();
test_saved_coefficient_operator_version_gate();
test_output_save_path_resolution();
test_compact_output_record();
test_automatic_output_filename();
test_automatic_mode_labels();
fprintf('All WNL tests passed.\n');
end

function test_model_specific_adjoint()
model.block = @model_adjoint_test_block;
spec = wnl_spec(0,0,0,2,'model_adjoint_test');
spec.directSeed = [1;0];
spec.useDirectSeedAsRefinementState = true;
opts.verbose = false;
opts.refineOperatingPointEigenpair = false;
opts.modeResidualTolerance = 1.0e-9;
opts.coefficientModeResidualTolerance = 1.0e-9;
mode = wnl_compute_mode(model,spec,opts);
assert(mode.directResidual < 1.0e-13);
assert(mode.leftResidual < 1.0e-13);
assert(mode.tracking.left.modelSolver.accepted);
assert(strcmp(mode.tracking.left.method, ...
    'validated model-specific adjoint'));
end

function test_programmatic_input_merge()
defaults.alpha = 1;
defaults.nested.keep = 2;
defaults.nested.replace = 3;
defaults.array = [1,2];
override.nested.replace = 4;
override.nested.added = 5;
override.array = [7,8,9];
merged = vi_wnl_merge_input(defaults,override);
assert(merged.alpha == 1);
assert(merged.nested.keep == 2);
assert(merged.nested.replace == 4);
assert(merged.nested.added == 5);
assert(isequal(merged.array,[7,8,9]));
end

function block = model_adjoint_test_block(~)
block.A = sparse([0,0;0,2]);
block.Bslow = speye(2);
block.solveAdjoint = @(direct,spec,opts) ...
    model_adjoint_test_solve(direct,spec,opts); %#ok<INUSD>
end

function result = model_adjoint_test_solve(~,~,~)
result.vector = [1;0];
result.diagnostics.available = true;
result.diagnostics.method = 'test exact adjoint';
end

function test_real_modal_coefficient_combinatorics()
config.omega = 0;
config.N = 0;
config.ndof = 1;
config.mass = 1;
config.linearFourier = @(spec,k) real_mode_linear(spec,k); %#ok<NASGU>
config.quadraticLocal = @(a,b,varargin) a.*b;
config.cubicLocal = @(a,b,c,varargin) a.*b.*c;
model = wnl_fourier_model(config);
spec = model.makeSpec(0,0,'real_harmonic_mode');
spec.lambda = 1;
spec.direct = 1;
spec.left = 1;
result = wnl_analyze_single_mode(model,spec,struct('verbose',false));
assert(strcmp(result.coordinateType,'real'));
assert(abs(result.self.qAA.vector-1) < 1.0e-13);
assert(isempty(result.self.qAbarA));
% q_AA=1, so 2*C(phi,q_AA)+D(phi,phi,phi)=2+1=3.
assert(abs(result.g-3) < 1.0e-12);
assert(abs(result.gPhysicalPeak-3) < 1.0e-12);
assert(result.amplitudeScaleToPeakZetaOverH == 1);
assert(isreal(result.g));
end

function test_saved_real_coefficient_reprojection()
mode = struct('coordinateType','real', ...
    'spec',struct('m',0,'s',0,'lambda',1));
entry = @(raw) struct('g',NaN,'gUnprojected',raw, ...
    'validCubicScaling',false,'forcedSolvesValid',true, ...
    'forcedSolvesExploratoryUsable',true, ...
    'quadraticResonance',false,'message','withheld');
saved = struct('modes',{{mode;mode}}, ...
    'self',{{entry(2+3.0e-6i);entry(-4+5.0e-6i)}}, ...
    'cross',{{[],entry(6+7.0e-6i);entry(-8+9.0e-6i),[]}}, ...
    'g',complex(nan(2)),'forcedSolvesValid',true(2), ...
    'forcedSolvesExploratoryUsable',true(2), ...
    'quadraticResonances',[], ...
    'referenceType','operating-point Floquet reduction', ...
    'slowEnvelopeValid',false);
projected = wnl_reproject_real_mode_set(saved, ...
    struct('realCoefficientImaginaryTolerance',5.0e-6));
assert(projected.numericalCoefficientValidity);
assert(norm(projected.g-[2,6;-8,-4],'fro') < 1.0e-14);
assert(projected.smallAmplitudeTransientAvailable);
assert(~projected.quantitativelyValid);

saved.forcedSolvesValid(1,2) = false;
didReject = false;
try
    wnl_reproject_real_mode_set(saved, ...
        struct('realCoefficientImaginaryTolerance',5.0e-6));
catch errorResult
    didReject = strcmp(errorResult.identifier, ...
        'wnl_reproject_real_mode_set:UnconvergedForcedField');
end
assert(didReject);
end

function test_self_conjugate_field_projection()
config.omega = 1;
config.N = 1;
config.ndof = 1;
config.mass = 1;
config.linearFourier = @(~,~) sparse(1);
model = wnl_fourier_model(config);
spec = model.makeSpec(0,0,'real_field_projection');
field = wnl_make_field(spec,[1+2.0e-8i,2,1-1.0e-8i]);
[projected,diagnostics] = wnl_project_self_conjugate_field( ...
    model,field,struct('realFieldConjugacyTolerance',1.0e-6));
reflected = wnl_conjugate_field(model,projected,'projected_bar');
assert(diagnostics.accepted);
assert(diagnostics.relativeDefect > 0);
assert(norm(wnl_field_vector(projected)- ...
    wnl_field_vector(reflected)) < 1.0e-14);
[~,strictDiagnostics] = wnl_project_self_conjugate_field( ...
    model,field,struct('realFieldConjugacyTolerance',1.0e-12));
assert(~strictDiagnostics.accepted);
end

function test_self_conjugate_forced_solution_projection()
config.omega = 0;
config.N = 0;
config.ndof = 2;
config.mass = sparse(2,2);
config.linearFourier = @(~,k) projected_solution_linear(k);
model = wnl_fourier_model(config);
spec = model.makeSpec(0,0,'real_forced_projection');
solution.field = wnl_make_field(spec,[1;1i]);
solution.forcing = [1;0];
solution.valid = true;
solution.quadraticResonance = false;
solution.lambda = complex(zeros(0,1));
[projected,diagnostics] = ...
    wnl_project_self_conjugate_forced_solution(model,solution, ...
    struct('realFieldConjugacyTolerance',1.0e-12, ...
    'forcedSolveResidualTolerance',1.0e-10));
assert(~diagnostics.rawConjugacyAccepted);
assert(diagnostics.acceptedByProjectedEquation);
assert(diagnostics.accepted);
assert(diagnostics.projectedRelativeEquationResidual < 1.0e-14);
assert(norm(wnl_field_vector(projected)-[1;0]) < 1.0e-14);
end

function value = projected_solution_linear(k)
if k == 0
    value = sparse(diag([-1,0]));
else
    value = sparse(2,2);
end
end

function value = real_mode_linear(~,k)
if k == 0
    value = sparse(1);
else
    value = sparse(1,1);
end
end

function test_complex_target_real_source_cross()
% A real axisymmetric source does not make a nonaxisymmetric target's sum
% field self-conjugate. Its A*B^2 coefficient must remain in the target's
% complex coordinate without a reality projection.
config.omega = 0;
config.N = 0;
config.ndof = 1;
config.mass = 1;
config.linearFourier = @(~,~) sparse(1,1);
config.quadraticLocal = @(a,b,varargin) a.*b;
config.cubicLocal = @(a,b,c,varargin) a.*b.*c;
model = wnl_fourier_model(config);

specA = model.makeSpec(1,0,'complex_target');
specA.lambda = 1;
modeA = struct('spec',specA,'field',wnl_make_field(specA,1), ...
    'vector',1,'left',1,'normalization',1,'coordinateType','complex');
specB = model.makeSpec(0,0,'real_source');
specB.lambda = 2;
modeB = struct('spec',specB,'field',wnl_make_field(specB,1), ...
    'vector',1,'left',1,'normalization',1,'coordinateType','real');

result = wnl_real_source_cross_coefficient( ...
    model,modeA,modeB,{},struct('verbose',false),struct());
% q_BB=1/4 and q_AB=2/3, hence
% 2*q_BB + 2*q_AB + 3 = 29/6.
assert(result.validCubicScaling);
assert(abs(result.g-29/6) < 1.0e-12);
assert(~result.fieldReality.qAB.applied);
assert(~result.fieldReality.cubicForcing.applied);
end

function test_two_real_mode_shared_fields()
config.omega = 0;
config.N = 0;
config.ndof = 2;
config.mass = speye(2);
config.linearFourier = @(spec,k) two_real_mode_linear(spec,k); %#ok<NASGU>
config.quadraticLocal = @shared_field_test_quadratic;
config.cubicLocal = @shared_field_test_cubic;
model = wnl_fourier_model(config);
spec1 = model.makeSpec(0,0,'real_mode_1');
spec1.lambda = 1;
spec1.direct = [1;0];
spec1.left = [1;0];
spec2 = model.makeSpec(0,0.5,'real_mode_2');
spec2.lambda = 3;
spec2.direct = [0;1;0;1]/sqrt(2);
spec2.left = [0;1;0;1]/sqrt(2);
opts = struct('verbose',false);
shared = wnl_analyze_mode_set(model,{spec1;spec2},opts);
assert(all(strcmp(shared.coordinateTypes,'real')));
assert(all(isfinite(shared.g(:))) && isreal(shared.g));
assert(shared.optimization.baselineForcedSolveCount == 6);
assert(shared.optimization.actualForcedSolveCount == 3);
assert(shared.optimization.reusedForcedSolveCount == 3);
assert(isempty(shared.self{1}.qAbarA));
assert(isempty(shared.cross{1,2}.qAbarB));

targeted = wnl_analyze_mode_set(model,{spec1;spec2}, ...
    struct('verbose',false,'coefficientPairs',[2,1]));
assert(isequal(targeted.coefficientComputed,logical([0,0;1,0])));
assert(isnan(targeted.g(1,1)) && isnan(targeted.g(1,2)) && ...
    isnan(targeted.g(2,2)));
assert(abs(targeted.g(2,1)-shared.g(2,1)) < 1.0e-12);
assert(isempty(targeted.self{1}) && isempty(targeted.self{2}));
assert(isempty(targeted.cross{1,2}));
assert(~isempty(targeted.cross{2,1}));
assert(targeted.optimization.targetedCoefficientEvaluation);
assert(targeted.optimization.baselineForcedSolveCount == 2);
assert(targeted.optimization.actualForcedSolveCount == 2);

duplicate = wnl_analyze_mode_set(model,{spec1;spec2}, ...
    struct('verbose',false,'reuseTwoModeForcedFields',false));
assert(norm(duplicate.g-shared.g,'fro') < 1.0e-12);
assert(duplicate.optimization.actualForcedSolveCount == 6);
end

function value = two_real_mode_linear(~,k)
if k == 0
    value = sparse(diag([1,3]));
else
    value = sparse(2,2);
end
end

function test_modal_amplitude_scaling()
realSpec = wnl_spec(0,0,0,1,'real_mode');
realMode = struct('spec',realSpec,'coordinateType','real');
complexSpec = wnl_spec(2,0,0,1,'complex_mode');
complexMode = struct('spec',complexSpec,'coordinateType','complex');
scales = wnl_mode_amplitude_scale({realMode;complexMode});
assert(isequal(scales,[1;2]));
gInternal = [3,8;5,12];
gPhysical = wnl_physical_cubic_coefficients( ...
    gInternal,{realMode;complexMode});
assert(isequal(gPhysical,[3,2;5,3]));

modes(1) = struct('m',0,'s',0,'label','real_one');
modes(2) = struct('m',2,'s',0,'label','complex_two');
settings.amplitudesOverH = [1e-3;2e-3];
settings.phases = [pi;pi/3];
initial = vi_mode_initial_conditions(modes,settings);
assert(initial.complexAmplitudesOverH(1) == -1e-3);
assert(abs(initial.complexAmplitudesOverH(2)- ...
    2e-3*exp(1i*pi/3)) < 1.0e-15);
settings.phases(1) = pi/2;
didReject = false;
try
    vi_mode_initial_conditions(modes,settings);
catch phaseError
    didReject = strcmp(phaseError.identifier, ...
        'vi_mode_initial_conditions:RealModePhase');
end
assert(didReject);

end

function test_radial_mode_correction()
r = linspace(0,1,17).';
dr = r(2)-r(1);
quadrature = dr*ones(size(r));
quadrature([1,end]) = dr/2;
metadata.layout.zeta = 1:numel(r);
metadata.discretization.r = r;
metadata.discretization.radial.quadratureWeights = quadrature;
metadata.contactLine = 'free';
metadata.sidewallTangentialCondition = 'stressFree';
parameters.R0 = 1;
parameters.omegaStar = 1;

spec = wnl_spec(2,0,[-1,0],numel(r),'radial_correction_test');
spec.betaStar = 3.4;
besselShape = besselj(2,spec.betaStar*r);
besselShape = besselShape/max(abs(besselShape));
weight = quadrature.*r;
correctionShape = r.^2.*(1-r).^2;
correctionShape = correctionShape-besselShape* ...
    ((besselShape'*(weight.*correctionShape))/ ...
    (besselShape'*(weight.*besselShape)));
correctionShape = correctionShape/sqrt( ...
    correctionShape'*(weight.*correctionShape));
zeta = besselShape*[1,0.3i]+ ...
    correctionShape*[0.2,-0.05i];
mode.spec = spec;
mode.field = wnl_make_field(spec,zeta);

settings.radialCorrectionPolicy = 'auto';
settings.radialCorrectionRelativeL2Threshold = 1.0e-2;
settings.radialCorrectionTemporalSamples = 64;
analysis = vi_radial_mode_correction( ...
    mode,metadata,parameters,settings);
assert(analysis.includeInSurfacePattern);
assert(analysis.relativeCorrectionL2 > 1.0e-2);
assert(analysis.orthogonalityResidual < 1.0e-12);
assert(norm(analysis.fullCoefficients-( ...
    analysis.separableCoefficients+analysis.correctionCoefficients), ...
    'fro') < 1.0e-13);
assert(norm(analysis.selectedCoefficients-zeta,'fro') < 1.0e-13);

settings.radialCorrectionRelativeL2Threshold = 2;
omitted = vi_radial_mode_correction( ...
    mode,metadata,parameters,settings);
assert(~omitted.includeInSurfacePattern);
assert(norm(omitted.selectedCoefficients- ...
    omitted.separableCoefficients,'fro') < 1.0e-13);

time = linspace(0,1,21).';
linearResult.timeStar = time;
linearResult.forcingPeriods = time/(2*pi);
linearResult.displacementOverH = 1.0e-3*ones(size(time));
linearResult.floquetOscillation = ones(size(time));
linearResult.betaStar = spec.betaStar;
wnlResult.mode = mode;
wnlResult.self = struct();
amplitude = 1.0e-3*ones(size(time));
settings.plotInterfaceDynamics = false;
settings.includeSlavedHarmonics = false;
settings.radialCorrectionRelativeL2Threshold = 1.0e-2;
included = vi_compare_interface_dynamics(linearResult,amplitude, ...
    wnlResult,metadata,parameters,settings,2,amplitude);
assert(included.includedRadialModeCorrections);
assert(norm(included.linearField-included.wnlPrimaryField,'fro') ...
    < 1.0e-13);
settings.radialCorrectionPolicy = 'never';
separable = vi_compare_interface_dynamics(linearResult,amplitude, ...
    wnlResult,metadata,parameters,settings,2,amplitude);
assert(~separable.includedRadialModeCorrections);
assert(norm(included.wnlPrimaryField-separable.wnlPrimaryField,'fro') ...
    > 1.0e-6);

axisymmetricSpec = spec;
axisymmetricSpec.m = 0;
axisymmetricSpec.label = 'axisymmetric_radial_correction_test';
axisymmetricMode.spec = axisymmetricSpec;
axisymmetricMode.field = wnl_make_field(axisymmetricSpec,zeta);
settings.radialCorrectionPolicy = 'auto';
axisymmetric = vi_radial_mode_correction( ...
    axisymmetricMode,metadata,parameters,settings);
assert(axisymmetric.exactlySeparableBoundaryModel);
assert(~axisymmetric.includeInSurfacePattern);
assert(axisymmetric.relativeCorrectionL2 > 1.0e-2);
end

function test_phase_sensitive_cubic_terms()
subharmonic = wnl_spec(2,0.5,0,1,'subharmonic');
harmonic = wnl_spec(2,0,0,1,'harmonic');
differentM = wnl_spec(3,0,0,1,'different_m');
assert(wnl_phase_sensitive_coupling_allowed(subharmonic,harmonic));
assert(wnl_phase_sensitive_coupling_allowed(harmonic,subharmonic));
assert(~wnl_phase_sensitive_coupling_allowed( ...
    subharmonic,differentM));

A = [0.2+0.1i;-0.3+0.4i];
lambda = [0.5;-0.2];
g = [-1,2;3,-4];
h = [0,0.7-0.2i;-0.4,0];
computed = wnl_rhs_landau(0,A,lambda,g,h);
expected = lambda.*A+A.*(g*abs(A).^2)+conj(A).*(h*(A.^2));
assert(norm(computed-expected) < 1.0e-14);

time = linspace(0,0.2,41).';
limits = struct('maximumAmplitude',Inf, ...
    'maximumRelativeCorrection',Inf, ...
    'maximumNonlinearToLinearRateRatio',Inf, ...
    'linearRateFloor',1.0e-12);
corrected = vi_cubic_transient_correction( ...
    time,A,lambda,g,limits,h);
linear = exp(time*lambda.').*A.';
expectedPhase = complex(zeros(size(linear)));
for target = 1:2
    for source = 1:2
        if h(target,source) == 0
            continue;
        end
        exponent = conj(lambda(target))+2*lambda(source)-lambda(target);
        integral = expm1(exponent*time)/exponent;
        expectedPhase(:,target) = expectedPhase(:,target)+ ...
            exp(lambda(target)*time).*h(target,source)* ...
            conj(A(target))*A(source)^2.*integral;
    end
end
assert(norm(corrected.phaseSensitiveCorrections-expectedPhase,'fro') ...
    < 1.0e-13);
end

function test_phase_sensitive_coefficient_projection()
config.omega = 0;
config.N = 1;
config.ndof = 2;
config.mass = speye(2);
config.linearFourier = @phase_coefficient_linear;
config.quadraticLocal = @(a,b,varargin) complex(zeros(2,1));
config.cubicLocal = @phase_coefficient_cubic;
model = wnl_fourier_model(config);
target = model.makeSpec(2,0.25,'target');
target.lambda = 1;
target.direct = [1;0;1;0;0;0]/sqrt(2);
target.left = target.direct;
source = model.makeSpec(2,-0.25,'source');
source.lambda = 3;
source.direct = [0;0;0;1;0;0];
source.left = source.direct;
result = wnl_analyze_mode_set( ...
    model,{target;source},struct('verbose',false));
% D(e1,e2,e2)=e1. The phase-sensitive monomial has three ordered
% permutations, whereas A1*|A2|^2 has six.
assert(abs(result.phaseSensitiveG(1,2)-3) < 1.0e-12);
assert(abs(result.g(1,2)-6) < 1.0e-12);
assert(result.hasPhaseSensitiveCoupling);
end

function value = phase_coefficient_linear(~,k)
if k == 0
    value = sparse(diag([1,3]));
else
    value = sparse(2,2);
end
end

function value = phase_coefficient_cubic(a,b,c,varargin) %#ok<INUSD>
value = [a(1)*b(2)*c(2)+a(2)*b(1)*c(2)+ ...
    a(2)*b(2)*c(1);0];
end

function test_cubic_mode_set_completeness_gate()
harmonic = wnl_spec(2,0,0,1,'harmonic');
subharmonic = wnl_spec(2,0.5,0,1,'subharmonic');
differentM = wnl_spec(3,0,0,1,'different_m');
assert(wnl_validate_cubic_mode_set({harmonic;subharmonic}).supported);
assert(wnl_validate_cubic_mode_set({harmonic;differentM}).supported);
sameClass = wnl_spec(2,0,0,1,'same_class');
didReject = false;
try
    wnl_validate_cubic_mode_set({harmonic;sameClass});
catch modeSetError
    didReject = strcmp(modeSetError.identifier, ...
        'wnl_validate_cubic_mode_set:IncompleteCubicTensor');
end
assert(didReject);

azimuthalOne = wnl_spec(1,0,0,1,'azimuthal_one');
azimuthalThree = wnl_spec(3,0,0,1,'azimuthal_three');
didReject = false;
try
    wnl_validate_cubic_mode_set({azimuthalOne;azimuthalThree});
catch modeSetError
    didReject = strcmp(modeSetError.identifier, ...
        'wnl_validate_cubic_mode_set:IncompleteCubicTensor');
end
assert(didReject);

azimuthalThreeSubharmonic = ...
    wnl_spec(3,0.5,0,1,'azimuthal_three_subharmonic');
assert(wnl_validate_cubic_mode_set( ...
    {azimuthalOne;azimuthalThreeSubharmonic}).supported);
end

function test_chebyshev_quadrature()
[x,~,weights] = vi_chebyshev_lobatto(16,[0,1]);
for power = 0:12
    exact = 1/(power+1);
    assert(abs(weights'*(x.^power)-exact) < 5.0e-14);
end
roots = bessel_derivative_root(0,6);
shape = besselj(0,roots(6)*x);
spectralVolume = abs(weights'*(x.*shape));
trapezoidalVolume = abs(trapz(x,x.*shape));
assert(spectralVolume < trapezoidalVolume);
assert(spectralVolume < 1.0e-6);
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
input.boundary.sidewallTangentialCondition = 'stressFree';
stressFreeFilename = vi_wnl_output_filename(input);
assert(strcmp(stressFreeFilename, ...
    'vi_wnl_ag0-2_fHz-30_modes-m2l6-m2l2_wall-sf.mat'));
input.boundary.sidewallTangentialCondition = 'legacyFreeSlide';
assert(endsWith(vi_wnl_output_filename(input),'_wall-lfs.mat'));
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
    'forcedUseRankAwareMinimumNorm',true, ...
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
assert(~solution.rankAwareDiagnostics.attempted);
end

function test_model_solver_only_forced_failure()
% Large model-specific solves may deliberately fail fast without launching
% the global minimum-norm and iterative-refinement fallbacks.
model.block = @(spec) model_only_forced_test_block(); %#ok<NASGU>
spec = wnl_spec(0,0,0,1,'model_only_forced_test');
opts = struct('verbose',false,'forcedModelSolverOnly',true, ...
    'forcedUseBlockGmres',true, ...
    'forcedUseRankAwareMinimumNorm',true, ...
    'forcedUseAlgebraicCompletion',true, ...
    'forcedSolveResidualTolerance',1.0e-6);
solution = wnl_solve_forced(model,spec,1,{},opts);
assert(~solution.valid && solution.relativeEquationResidual == 1);
assert(solution.modelSolveDiagnostics.attempted);
assert(~solution.modelSolveDiagnostics.passedPhysicalGate);
assert(~solution.gmresDiagnostics.attempted);
assert(~solution.rankAwareDiagnostics.attempted);
assert(solution.solveDiagnostics.totalIterations == 0);
assert(~solution.algebraicCompletion.attempted);
end

function block = model_only_forced_test_block()
diagnostics = struct('attempted',true,'available',true, ...
    'passedPhysicalGate',false,'relativeResidual',1, ...
    'method','manufactured model solver');
block = struct('A',sparse(1,1),'Bslow',sparse(1,1), ...
    'solve',@(forcing,spec,opts) struct( ... %#ok<INUSD>
        'vector',complex(0),'diagnostics',diagnostics));
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
assert(any(strcmp(solution.modelSolveDiagnostics.blockFactorMethods, ...
    'equilibrated-lu')));
assert(solution.relativeEquationResidual < 1.0e-11);
assert(norm(solution.vector-exact)/norm(exact) < 1.0e-10);

opts.forcedCylinderSchurFactorMethod = 'columnQr';
opts.forcedModelSolverOnly = true;
lazyModel = model;
lazyModel.block = @(requestedSpec) error( ...
    'run_wnl_tests:UnexpectedBlockAssembly', ...
    'The model-only QR path assembled the complete block.'); %#ok<NASGU>
qrSolution = wnl_solve_forced(lazyModel,spec,forcing,{},opts);
assert(qrSolution.valid);
assert(any(strcmp(qrSolution.modelSolveDiagnostics.blockFactorMethods, ...
    'column-equilibrated-qr')));
assert(qrSolution.relativeEquationResidual < 1.0e-11);
assert(norm(qrSolution.vector-exact)/norm(exact) < 1.0e-10);
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

function test_cylinder_temporal_schur_pressure_nullspace()
% A pressure-only null coordinate must be fixed by a minimum-norm border,
% not by shifting all physical equations. The manufactured forcing is in
% the range of the singular complete Floquet operator.
ndof = 3;
blocks.B0 = sparse(diag([1,1,0]));
blocks.L0 = sparse(diag([-1.2,-2.1,0]));
blocks.Lplus = sparse(1,2,0.08,ndof,ndof);
blocks.Lminus = sparse(2,2,-0.05,ndof,ndof);
omega = 1.3;
spec = wnl_spec(0,0,[-1,0,1],ndof,'schur_pressure_null_test');
spec.lambda = 0.4;
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
exact = [1;2;0;3;4;0;5;6;0];
forcing = assembled.A*exact;
opts = struct('verbose',false,'forcedUseCylinderSchur',true, ...
    'forcedUseBlockGmres',false, ...
    'forcedUseRankAwareMinimumNorm',false, ...
    'forcedUseAlgebraicCompletion',false, ...
    'forcedCylinderSchurUseNullspaceBordering',true, ...
    'forcedCylinderSchurMaximumNullity',2, ...
    'forcedSolveResidualTolerance',1.0e-11, ...
    'forcedSolveMaxRestarts',0);
solution = wnl_solve_forced(model,spec,forcing,{},opts);
diagnostics = solution.modelSolveDiagnostics;
assert(solution.valid && diagnostics.passedPhysicalGate);
assert(all(diagnostics.blockNullityRange == 1));
assert(any(strcmp(diagnostics.blockFactorMethods, ...
    'fixed-range-bordered-lu')));
assert(norm(solution.vector-exact) < 1.0e-9*max(norm(exact),1));
end

function test_cylinder_temporal_schur_pressure_quotient()
% A pressure-only left-compatibility defect may be removed only after its
% right nullspace is certified algebraic. Use a tiny nonzero pressure pivot
% so the matrix is structurally full rank: the certified nullity hint must
% still take precedence over an unconstrained LU, which would manufacture a
% huge pressure state. The complete residual is retained diagnostically.
ndof = 2;
blocks.B0 = sparse(diag([1,0]));
blocks.L0 = sparse(diag([-1,-1.0e-14]));
blocks.Lplus = sparse(ndof,ndof);
blocks.Lminus = sparse(ndof,ndof);
omega = 1;
spec = wnl_spec(0,0,0,ndof,'schur_pressure_quotient_test');
spec.lambda = 0;
config.omega = omega;
config.N = 0;
config.ndof = ndof;
config.mass = blocks.B0;
config.linearFourier = @(requestedSpec,k) ...
    schur_test_fourier(blocks,requestedSpec,k);
config.blockSolve = @(requestedSpec,forcing,options) ...
    vi_cylinder_wnl_forced_schur_solve( ...
    blocks,omega,requestedSpec,forcing,options);
model = wnl_fourier_model(config);
opts = struct('verbose',false,'forcedUseCylinderSchur',true, ...
    'forcedCylinderSchurFactorMethod','columnQr', ...
    'forcedCylinderSchurUsePressureCompatibilityProjection',true, ...
    'forcedCylinderSchurPressureCompatibilityTolerance',1.0e-5, ...
    'forcedModelSolverOnly',true, ...
    'forcedSolveResidualTolerance',1.0e-8);
solution = wnl_solve_forced(model,spec,[2;1.0e-5],{},opts);
diagnostics = solution.modelSolveDiagnostics;
assert(solution.valid && diagnostics.passedPhysicalGate);
assert(diagnostics.pressureNullspaceCertified);
assert(sprank(model.block(spec).A) == ndof);
assert(all(diagnostics.blockNullityRange == 1));
assert(any(strcmp(diagnostics.blockFactorMethods, ...
    'fixed-range-bordered-lu')));
assert(diagnostics.pressureCompatibilityProjectionApplied);
assert(solution.relativeEquationResidual > ...
    opts.forcedSolveResidualTolerance);
assert(solution.acceptanceRelativeEquationResidual < 1.0e-12);
assert(abs(solution.vector(2)) < 1.0e-12);

opts.forcedCylinderSchurUsePressureCompatibilityProjection = false;
unprojected = wnl_solve_forced(model,spec,[2;1.0e-5],{},opts);
assert(~unprojected.valid);

opts.forcedCylinderSchurUsePressureCompatibilityProjection = true;
oversized = wnl_solve_forced(model,spec,[2;1.0e-3],{},opts);
assert(~oversized.valid);
assert(~oversized.modelSolveDiagnostics. ...
    pressureCompatibilityGatePassed);
end

function test_cylinder_temporal_schur_adjoint_pressure_nullspace()
% Three pressure-only null coordinates must not contaminate the unique
% physical adjoint of the coupled temporal system.
ndof = 2;
coupling = 0.1;
blocks.B0 = sparse(diag([1,0]));
blocks.L0 = sparse(diag([1,0]));
blocks.Lplus = sparse(1,1,coupling,ndof,ndof);
blocks.Lminus = blocks.Lplus;
spec = wnl_spec(0,0,[-1,0,1],ndof,'adjoint_pressure_null_test');
spec.lambda = 1+sqrt(2)*coupling;
directByHarmonic = [sqrt(2)/2,1,sqrt(2)/2;zeros(1,3)];
direct = directByHarmonic(:);
opts = struct('verbose',false, ...
    'adjointCylinderSchurNullTolerance',1.0e-10, ...
    'adjointCylinderSchurFactorResidualTolerance',1.0e-9);
solution = vi_cylinder_wnl_adjoint_schur_solve( ...
    blocks,0,spec,direct,opts);
config.omega = 0;
config.N = 1;
config.ndof = ndof;
config.mass = blocks.B0;
config.linearFourier = @(requestedSpec,k) ...
    schur_test_fourier(blocks,requestedSpec,k);
model = wnl_fourier_model(config);
assembled = model.block(spec);
relativeResidual = norm(assembled.A'*solution.vector) / ...
    max(norm(assembled.A,1)*norm(solution.vector),eps);
assert(solution.diagnostics.available);
assert(all(solution.diagnostics.blockNullityRange == 1));
assert(solution.diagnostics.maximumBorderMultiplierRatio < 1.0e-10);
assert(relativeResidual < 1.0e-11);
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

function test_compact_output_record()
spec = wnl_spec(0,0,0,2,'compact_storage_test');
spec.directSeed = [1;2];
spec.direct = [1;2];
spec.left = [0.5;0.25];
mode.spec = spec;
mode.block = struct('A',speye(2),'Bslow',speye(2),'spec',spec, ...
    'solve',[],'solveAdjoint',[]);
mode.vector = [1;2];
mode.field = wnl_make_field(spec,mode.vector);
mode.left = [0.5;0.25];
mode.leftField = wnl_make_field(spec,mode.left);
mode.normalization = mode.left'*mode.vector;
mode.directResidual = 1.0e-10;
mode.leftResidual = 2.0e-10;
mode.coordinateType = 'real';

forced.spec = spec;
forced.field = wnl_make_field(spec,[3;4]);
forced.vector = forced.field.coeff(:);
forced.forcing = [5;6];
forced.forcingNorm = norm(forced.forcing);
forced.relativeEquationResidual = 3.0e-10;
forced.valid = true;
self = struct('g',2,'gUnprojected',2,'qAA',forced, ...
    'qAbarA',[],'termSecondHarmonic',[7;8], ...
    'termDirectCubic',[9;10],'cubicForcing',[16;18], ...
    'forcedSolvesValid',true,'forcedSolvesExploratoryUsable',true, ...
    'validCubicScaling',true,'quadraticResonance',false);

wnl.modes = {mode};
wnl.conjugateModes = {mode};
wnl.neutralModes = {mode};
wnl.self = {self};
wnl.cross = cell(1);
phaseSensitive = struct('g',3,'gUnprojected',3, ...
    'qSourceSource',forced,'qTargetBarSource',forced, ...
    'termMixed',[1;2],'termSecondHarmonic',[3;4], ...
    'termDirectCubic',[5;6],'cubicForcing',[9;12], ...
    'forcedSolvesValid',true,'forcedSolvesExploratoryUsable',true, ...
    'validCubicScaling',true,'quadraticResonance',false);
wnl.phaseSensitive = {phaseSensitive};
wnl.g = 2;
wnl.gPhysicalPeak = 2;
output.weaklyNonlinear = wnl;
[compact,information] = vi_compact_output_record(output,'compact');

assert(strcmp(information.profile,'compact'));
assert(information.estimatedBytesAfter < ...
    information.estimatedBytesBefore);
assert(~isfield(compact.weaklyNonlinear,'conjugateModes'));
assert(~isfield(compact.weaklyNonlinear,'neutralModes'));
savedMode = compact.weaklyNonlinear.modes{1};
assert(~isfield(savedMode.block,'A'));
assert(isfield(savedMode.block,'Bslow'));
assert(~isfield(savedMode,'leftField'));
assert(~isfield(savedMode.spec,'directSeed'));
assert(isequal(savedMode.vector,mode.vector));
assert(isequal(savedMode.left,mode.left));
savedSelf = compact.weaklyNonlinear.self{1};
assert(~isfield(savedSelf,'termSecondHarmonic'));
assert(~isfield(savedSelf,'cubicForcing'));
assert(isfield(savedSelf,'storageProjectionDiagnostics'));
assert(~isfield(savedSelf.qAA,'vector'));
assert(~isfield(savedSelf.qAA,'forcing'));
assert(isequal(savedSelf.qAA.field.coeff,forced.field.coeff));
savedPhaseSensitive = compact.weaklyNonlinear.phaseSensitive{1};
assert(~isfield(savedPhaseSensitive,'termMixed'));
assert(~isfield(savedPhaseSensitive,'cubicForcing'));
assert(~isfield(savedPhaseSensitive,'qSourceSource'));
assert(~isfield(savedPhaseSensitive,'qTargetBarSource'));
assert(numel(savedPhaseSensitive.storageForcedFieldDiagnostics) == 2);
assert(compact.storage.capabilities.savedCoefficientPostprocessing);
assert(compact.storage.capabilities.modeRecoveryRestart);
assert(~compact.storage.capabilities.fullOperatorForensics);
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
spec2 = model.makeSpec(4,0,'shared_mode_2');
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
elseif ismember(abs(spec.m),[1,4])
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
    'betaStar',2.5,'reducedReferenceLambda',0.12);
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
output.operatorMetadata.linearOperatorVersion = ...
    'V2-pressure-compatible-sidewall';
output.operatorMetadata.adjointSolverVersion = ...
    'V2-pressure-nullspace-bordered-temporal-schur';
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

output.weaklyNonlinear.mode = mode;
output.weaklyNonlinear.mode.directResidual = 1.0e-11;
output.operatorMetadata = rmfield( ...
    output.operatorMetadata,'linearOperatorVersion');
save(cacheFile,'output');
[legacyWarmStart,legacyInformation] = vi_wnl_apply_recovery_cache( ...
    cacheFile,temporaryRoot,input,parameters,{spec},2);
assert(legacyInformation.warmStarted && ~legacyInformation.used);
assert(isequal(legacyWarmStart{1}.directLiftInitialVector,mode.vector));

changed = parameters;
changed.aAnalysis = 0.8;
[~,rejected] = vi_wnl_apply_recovery_cache( ...
    cacheFile,temporaryRoot,input,changed,{spec},2);
assert(~rejected.used && contains(rejected.reason,'aAnalysis'));

stressInput = input;
stressInput.boundary.sidewallTangentialCondition = 'stressFree';
stressSpec = rmfield(spec,'reducedReferenceLambda');
stressMode = mode;
stressMode.spec = stressSpec;
stressMode.directResidual = 1.0e-11;
output.input = stressInput;
output.parameters = parameters;
output.operatorMetadata.ndof = 2;
output.operatorMetadata.linearOperatorVersion = ...
    'V3-pressure-compatible-stress-free-sidewall';
output.weaklyNonlinear.mode = stressMode;
save(cacheFile,'output');
[stressUpdated,stressInformation] = vi_wnl_apply_recovery_cache( ...
    cacheFile,temporaryRoot,stressInput,parameters,{stressSpec},2);
assert(stressInformation.used);
assert(~isfield(stressUpdated{1},'reducedReferenceLambda'));

axisymmetricSavedInput = input;
axisymmetricSavedInput.modes.m = 0;
axisymmetricCurrentInput = axisymmetricSavedInput;
axisymmetricCurrentInput.boundary.sidewallTangentialCondition = ...
    'stressFree';
axisymmetricSpec = spec;
axisymmetricSpec.m = 0;
axisymmetricSpec.exactSeparableBesselInterface = true;
axisymmetricSpec.preservePrescribedDofsDuringRefinement = true;
axisymmetricMode = mode;
axisymmetricMode.spec = axisymmetricSpec;
axisymmetricMode.spec.lambda = 0.121;
output.input = axisymmetricSavedInput;
output.operatorMetadata.linearOperatorVersion = ...
    'V2-pressure-compatible-sidewall';
output.weaklyNonlinear.mode = axisymmetricMode;
save(cacheFile,'output');
[axisymmetricWarm,axisymmetricInformation] = ...
    vi_wnl_apply_recovery_cache(cacheFile,temporaryRoot, ...
    axisymmetricCurrentInput,parameters,{axisymmetricSpec},2);
assert(axisymmetricInformation.warmStarted);
assert(~axisymmetricInformation.used);
assert(contains(axisymmetricInformation.reason,'operator changed'));
assert(isequal(axisymmetricWarm{1}.directLiftInitialVector, ...
    axisymmetricMode.vector));
assert(abs(axisymmetricWarm{1}.lambda- ...
    axisymmetricSpec.reducedReferenceLambda) < eps);
end

function test_saved_coefficient_operator_version_gate()
temporaryRoot = tempname;
[created,message] = mkdir(temporaryRoot);
assert(created,message);
cleanup = onCleanup(@() remove_test_directory(temporaryRoot)); %#ok<NASGU>
cacheFile = fullfile(temporaryRoot,'legacy-cache.mat');
output = struct();
save(cacheFile,'output');
didRejectLegacy = false;
try
    vi_saved_small_amplitude_transient(cacheFile,struct());
catch cacheError
    didRejectLegacy = strcmp(cacheError.identifier, ...
        'vi_saved_small_amplitude_transient:LegacyLinearOperator');
end
assert(didRejectLegacy);

output.operatorMetadata.linearOperatorVersion = ...
    'V2-pressure-compatible-sidewall';
save(cacheFile,'output');
passedVersionGate = false;
try
    vi_saved_small_amplitude_transient(cacheFile,struct());
catch cacheError
    passedVersionGate = strcmp(cacheError.identifier, ...
        'vi_saved_small_amplitude_transient:MissingField');
end
assert(passedVersionGate);

output.operatorMetadata.linearOperatorVersion = ...
    'V3-pressure-compatible-stress-free-sidewall';
save(cacheFile,'output');
passedStressFreeVersionGate = false;
try
    vi_saved_small_amplitude_transient(cacheFile,struct());
catch cacheError
    passedStressFreeVersionGate = strcmp(cacheError.identifier, ...
        'vi_saved_small_amplitude_transient:MissingField');
end
assert(passedStressFreeVersionGate);
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
assert(defaultSummary.forcedCylinderSchurUseNullspaceBordering);
assert(defaultSummary.adaptiveDirectionalSteps);
assert(isequal(defaultSummary.quadraticStepMultipliers,[1,2,4,8]));
assert(isequal(defaultSummary.cubicStepMultipliers,[1,2,4]));
assert(defaultSummary.directionalStepEarlyStopTolerance == 1.0e-8);
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
spec1.betaStar = 1;
spec2.betaStar = 1;
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
    amplitudes, wnlResult, metadata, parameters, plotSettings, [2;2]);
assert(result.numberOfModes == 2);
assert(numel(result.secondOrderLabels) == 6);
assert(isequal(size(result.modeAmplitudesOverH), size(amplitudes)));
assert(norm(result.wnlTotalField-(result.wnlPrimaryField+ ...
    result.wnlSecondOrderField), 'fro') < 1e-14);
twoModeLinear = vi_compare_interface_dynamics_modes(linearResult, ...
    amplitudes, wnlResult, metadata, parameters, plotSettings, [2;2], ...
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

phaseMode1 = mode1;
phaseMode1.field = wnl_make_field(spec1,1i*besselj(1,r));
phaseWnlResult = struct('modes',{{phaseMode1;mode2}}, ...
    'self',{{struct();struct()}},'cross',{cell(2,2)});
phaseSettings = plotSettings;
phaseSettings.includeSlavedHarmonics = false;
phaseResult = vi_compare_interface_dynamics_modes(linearResult, ...
    amplitudes,phaseWnlResult,metadata,parameters,phaseSettings,[2;2], ...
    amplitudes);
assert(phaseResult.phaseAlignment(1) == 1);
assert(norm(phaseResult.linearModeProbeSignals(:,1)) < 1.0e-14);
assert(norm(phaseResult.wnlPrimaryModeProbeSignals(:,1)) < 1.0e-14);
end

function test_mixed_coordinate_cross_reconstruction()
r = [0;1];
time = [0;1];
spec1 = wnl_spec(0,0,0,numel(r),'real_mode');
spec2 = wnl_spec(2,0,0,numel(r),'complex_mode');
spec1.betaStar = 1;
spec2.betaStar = 1;
mode1 = struct('spec',spec1,'coordinateType','real', ...
    'field',wnl_make_field(spec1,zeros(size(r))));
mode2 = struct('spec',spec2,'coordinateType','complex', ...
    'field',wnl_make_field(spec2,zeros(size(r))));
self1 = struct();
self2 = struct();
cross12.qAB.field = wnl_make_field( ...
    wnl_spec(2,0,0,numel(r),'sum'),ones(size(r)));
cross12.qAbarB.field = wnl_make_field( ...
    wnl_spec(-2,0,0,numel(r),'difference'),ones(size(r)));
cross = cell(2,2);
cross{1,2} = cross12;
wnlResult = struct('modes',{{mode1;mode2}}, ...
    'self',{{self1;self2}},'cross',{cross});

linearResult.timeStar = time;
linearResult.forcingPeriods = time;
linearResult.displacementOverH = zeros(size(time));
linearResult.floquetOscillation = ones(size(time));
linearResult.betaStar = 1;
metadata.layout.zeta = 1:numel(r);
metadata.discretization.r = r;
parameters.omegaStar = 1;
parameters.R0 = 1;
settings.plotInterfaceDynamics = false;
settings.includeSlavedHarmonics = true;
amplitudes = [ones(size(time)),2*ones(size(time))];
result = vi_compare_interface_dynamics_modes(linearResult,amplitudes, ...
    wnlResult,metadata,parameters,settings,[1;2],amplitudes);
% Internal coordinates are a_1=a_2=1. q_12 and q_1bar2 are already the
% two conjugate members, so their total physical O(A1*A2) field is 2.
assert(norm(result.wnlSecondOrderField-2,'fro') < 1.0e-14);
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
modeSpec.betaStar = 1;
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
    wnlResult, metadata, parameters, settings, 2.0);
assert(primaryOnly.fieldRmse < 1.0e-14);
assert(abs(primaryOnly.normalizedPrimaryCarrierOverlap-1) < 1.0e-14);

settings.includeSlavedHarmonics = true;
withSlaved = vi_compare_interface_dynamics(linearResult, amplitude, ...
    wnlResult, metadata, parameters, settings, 2.0);
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

function test_prescribed_interface_refinement_constraint()
model = struct();
model.block = @separable_refinement_test_block;
spec = wnl_spec(0,0,0,2,'separable_refinement_test');
spec.lambda = 0.2;
spec.directSeed = [1;0];
spec.prescribedDofs = 1;
spec.preservePrescribedDofsDuringRefinement = true;
spec.useDirectSeedAsRefinementState = true;
opts = struct();
opts.verbose = false;
opts.modeTrackingRegularization = 0;
opts.refineOperatingPointEigenpair = true;
opts.eigenpairUseBlockGmres = true;
opts.eigenpairLsqrFallbackMaxIterations = 100;
opts.eigenpairLsqrFallbackMaxRestarts = 0;
opts.eigenpairCorrectionRegularization = 0;
opts.eigenpairRefinementMaxSteps = 2;
opts.eigenpairRefinementTolerance = 1.0e-12;
opts.eigenpairMaximumSeedResidualRatio = 1.0e12;
opts.stopEigenpairOnReferenceMismatch = false;
opts.autoEscalateEigenpairRefinement = false;
mode = wnl_compute_mode(model,spec,opts);
assert(abs(wnl_spec_lambda(mode.spec)-0.1) < 1.0e-10);
assert(mode.directResidual < 1.0e-10);
assert(abs(mode.vector(1)-1) < 1.0e-13);
assert(mode.tracking.eigenpairRefinement.steps == 1);
assert(strcmp(mode.tracking.eigenpairRefinement.solverHistory{1}, ...
    'regularized LSQR fallback'));
assert(mode.tracking.eigenpairRefinement. ...
    lastCorrectionInformation.preservedPrescribedDofs);
assert(mode.tracking.eigenpairRefinement. ...
    lastCorrectionInformation.numberOfPreservedDofs == 1);
end

function block = separable_refinement_test_block(spec)
lambda = wnl_spec_lambda(spec);
block.A = sparse([lambda-0.1,0;-2,1]);
block.Bslow = sparse([1,0;0,0]);
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
assert(sprank(0.37*B0-L0) == metadata.ndof);
assert(strcmp(metadata.linearOperatorVersion, ...
    'V2-pressure-compatible-sidewall'));
assert(strcmp(metadata.radialGrid.typeUsed, 'besselEnriched'));
assert(metadata.verticalGrid.lower.numberOfElements == 3);
assert(metadata.verticalGrid.upper.numberOfElements == 3);
defaultGaugeIr = max(2,min(metadata.layout.nr-1, ...
    ceil(metadata.layout.nr/2)));
assert(metadata.discretization.denseGaugeRadialIndex == defaultGaugeIr);
assert(metadata.pressureGauge.radialIndex == defaultGaugeIr);

gaugeParameters = p;
gaugeParameters.numerics.pressureGaugeRadialIndex = 3;
[gaugeOperators,gaugeMetadata] = cylinder_wnl_operators(gaugeParameters);
assert(gaugeMetadata.discretization.denseGaugeRadialIndex == 3);
assert(gaugeMetadata.pressureGauge.radialIndex == 3);
specGauge = wnl_spec(0,0,0,gaugeMetadata.ndof,'gauge_test');
[~,~,~,~,~,~,gaugeDiagnostics] = ...
    vi_cylinder_wnl_linear_operators( ...
    gaugeMetadata.discretization,specGauge.m);
gaugeNode = 3+(gaugeMetadata.discretization. ...
    denseGaugeVerticalIndex-1)*gaugeMetadata.layout.nr;
expectedGaugeRow = gaugeMetadata.layout.d.p(gaugeNode);
gaugeB0 = gaugeOperators.B0(specGauge);
assert(gaugeDiagnostics.gaugeRow == expectedGaugeRow);
assert(nnz(gaugeB0(expectedGaugeRow,:)) == 0);

% The legacy plus sign preserves the potential-flow Neumann Bessel seed.
% The true cylindrical stress-free minus sign must reject that separable
% m>0 velocity field and thereby permit a radial correction.
stressFreeParameters = p;
stressFreeParameters.boundary.sidewallTangentialCondition = 'stressFree';
[stressFreeOperators,stressFreeMetadata] = ...
    cylinder_wnl_operators(stressFreeParameters);
stressFreeB0 = stressFreeOperators.B0(specA);
stressFreeL0 = stressFreeOperators.Lhat(specA,0);
assert(strcmp(stressFreeMetadata.sidewallTangentialCondition, ...
    'stressFree'));
assert(strcmp(stressFreeMetadata.linearOperatorVersion, ...
    'V3-pressure-compatible-stress-free-sidewall'));
assert(sprank(0.37*stressFreeB0-stressFreeL0) == ...
    stressFreeMetadata.ndof);
azimuthalVelocity = zeros(metadata.layout.nr,1);
positiveRadius = metadata.discretization.r > 0;
azimuthalVelocity(positiveRadius) = p.modes.m* ...
    besselj(p.modes.m,p.modes.betaStar* ...
    metadata.discretization.r(positiveRadius))./ ...
    metadata.discretization.r(positiveRadius);
besselVelocityState = zeros(metadata.ndof,1);
for verticalIndex = 1:metadata.layout.nzD
    radialNodes = (1:metadata.layout.nr) + ...
        (verticalIndex-1)*metadata.layout.nr;
    besselVelocityState(metadata.layout.d.ut(radialNodes)) = ...
        azimuthalVelocity;
end
for verticalIndex = 1:metadata.layout.nzL
    radialNodes = (1:metadata.layout.nr) + ...
        (verticalIndex-1)*metadata.layout.nr;
    besselVelocityState(metadata.layout.l.ut(radialNodes)) = ...
        azimuthalVelocity;
end
sidewallUtRows = [metadata.layout.d.ut( ...
    metadata.layout.nr:metadata.layout.nr:metadata.layout.nD), ...
    metadata.layout.l.ut( ...
    metadata.layout.nr:metadata.layout.nr:metadata.layout.nL)];
legacyBesselWallResidual = norm( ...
    L0(sidewallUtRows,:)*besselVelocityState);
stressFreeBesselWallResidual = norm( ...
    stressFreeL0(sidewallUtRows,:)*besselVelocityState);
assert(legacyBesselWallResidual < ...
    1.0e-8*max(stressFreeBesselWallResidual,1));
assert(stressFreeBesselWallResidual > 1.0e-6);
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
fieldMinusSpec = wnl_spec(-2,0,0,gaugeMetadata.ndof, ...
    'gauge_source_minus');
fieldGaugeOut = wnl_spec(0,0,0,gaugeMetadata.ndof, ...
    'gauge_output');
fieldMinus = wnl_make_field(fieldMinusSpec,y);
gaugeQuadratic = gaugeOperators.CField( ...
    fieldX,fieldMinus,fieldGaugeOut);
assert(norm(gaugeQuadratic(expectedGaugeRow,:),'fro') == 0);
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

sideTangentialRows = zeros(0,1);
sideNormalRows = zeros(0,1);
sidePressureRows = zeros(0,1);
for verticalIndex = 1:metadata.layout.nzD
    sideNode = metadata.layout.nr + ...
        (verticalIndex-1)*metadata.layout.nr;
    sideTangentialRows = [sideTangentialRows; ... %#ok<AGROW>
        metadata.layout.d.ut(sideNode); ...
        metadata.layout.d.w(sideNode)];
    sideNormalRows(end+1,1) = ... %#ok<AGROW>
        metadata.layout.d.ur(sideNode);
    sidePressureRows(end+1,1) = ... %#ok<AGROW>
        metadata.layout.d.p(sideNode);
end
for verticalIndex = 1:metadata.layout.nzL
    sideNode = metadata.layout.nr + ...
        (verticalIndex-1)*metadata.layout.nr;
    sideTangentialRows = [sideTangentialRows; ... %#ok<AGROW>
        metadata.layout.l.ut(sideNode); ...
        metadata.layout.l.w(sideNode)];
    sideNormalRows(end+1,1) = ... %#ok<AGROW>
        metadata.layout.l.ur(sideNode);
    sidePressureRows(end+1,1) = ... %#ok<AGROW>
        metadata.layout.l.p(sideNode);
end
assert(norm(fieldQuadratic(sideTangentialRows,:),'fro') == 0);
assert(nnz(B0(sidePressureRows,:)) == 0);
assert(all(full(sum(abs(L0(sidePressureRows,:)),2)) > 0));
assert(norm(fieldQuadratic(sidePressureRows,:),'fro') > 1.0e-12);

% A free contact line has zeta_r=0, so the mapped free-slide sidewall rows
% are linear. A pinned line fixes zeta but permits zeta_r~=0 and must retain
% the quadratic ALE radial-derivative correction.
p.boundary.contactLine = 'pinned';
[pinnedOperators,pinnedMetadata] = cylinder_wnl_operators(p);
assert(pinnedMetadata.ndof == metadata.ndof);
xPinned = x;
yPinned = y;
xPinned(metadata.layout.zeta(end)) = 0;
yPinned(metadata.layout.zeta(end)) = 0;
pinnedX = wnl_make_field(fieldSpec,xPinned);
pinnedY = wnl_make_field(fieldSpec,yPinned);
pinnedQuadratic = pinnedOperators.CField(pinnedX,pinnedY,fieldOut);
assert(norm(pinnedQuadratic(sideTangentialRows,:),'fro') > 1.0e-12);
assert(norm(pinnedQuadratic(sideNormalRows,:),'fro') == 0);

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

% A fixed reference library must produce exactly the same radial operator
% when the dynamically retained mode changes.
referenceModes = repmat(p.modes,2,1);
referenceModes(1) = p.modes;
referenceModes(2).m = 0;
referenceModes(2).radialIndex = 2;
referenceRoots = bessel_derivative_root(0,2);
referenceModes(2).betaStar = referenceRoots(2)/p.R0;
p.numerics.radialGrid.referenceModes = referenceModes;
p.modes = referenceModes(1);
fixedA = vi_cylinder_radial_grid(p);
p.modes = referenceModes(2);
fixedB = vi_cylinder_radial_grid(p);
assert(fixedA.modeSetIndependent && fixedB.modeSetIndependent);
assert(strcmp(fixedA.basisModeSource,'referenceModes'));
assert(isequal(fixedA.selectedBasisLabels, ...
    fixedB.selectedBasisLabels));
assert(isequal(fixedA.D,fixedB.D));
assert(isequal(fixedA.D2,fixedB.D2));
assert(numel(fixedA.basisModes) == 2);

% A refined grid can preserve every function selected on a coarser grid,
% making the enriched spaces nested instead of repivoting from scratch.
p.modes = referenceModes(1);
p.numerics.Nr = 11;
p.numerics.radialGrid.preferredBasisLabels = ...
    fixedA.selectedBasisLabels;
nested = vi_cylinder_radial_grid(p);
assert(all(ismember(fixedA.selectedBasisLabels, ...
    nested.selectedBasisLabels)));
assert(numel(nested.selectedBasisLabels) == 11);

p.numerics.Nr = 9;
p.numerics.radialGrid.type = 'chebyshev';
legacy = vi_cylinder_radial_grid(p);
assert(strcmp(legacy.typeUsed, 'chebyshev'));
assert(legacy.modeSetIndependent);
assert(norm(legacy.r-chebyshevR, inf) < 10*eps);
assert(norm(legacy.D-chebyshevD, inf) < 1.0e-12);
end

function test_radial_mode_prolongation()
sourceNr = 6;
targetNr = 8;
nzD = 2;
nzL = 3;
[sourceR,~,~] = vi_chebyshev_lobatto(sourceNr,[0,1]);
[targetR,~,~] = vi_chebyshev_lobatto(targetNr,[0,1]);
sourceLayout = test_cylinder_layout(sourceNr,nzD,nzL);
targetLayout = test_cylinder_layout(targetNr,nzD,nzL);
numberOfHarmonics = 2;
coeff = complex(zeros(sourceLayout.ndof,numberOfHarmonics));
expected = complex(zeros(targetLayout.ndof,numberOfHarmonics));
variables = {'ur','ut','w','p'};
for variableIndex = 1:numel(variables)
    name = variables{variableIndex};
    [coeff(sourceLayout.d.(name),:), ...
        expected(targetLayout.d.(name),:)] = polynomial_block( ...
        sourceR,targetR,nzD,numberOfHarmonics,variableIndex);
    [coeff(sourceLayout.l.(name),:), ...
        expected(targetLayout.l.(name),:)] = polynomial_block( ...
        sourceR,targetR,nzL,numberOfHarmonics,variableIndex+4);
end
[coeff(sourceLayout.zeta,:),expected(targetLayout.zeta,:)] = ...
    polynomial_block(sourceR,targetR,1,numberOfHarmonics,9);
mode.field.coeff = coeff;
metadata.layout = sourceLayout;
metadata.radialGrid.r = sourceR;
[seed,diagnostics] = vi_prolong_cylinder_mode_radially( ...
    mode,metadata,targetR);
actual = reshape(seed,targetLayout.ndof,numberOfHarmonics);
assert(norm(actual-expected,'fro') < 1.0e-10);
assert(diagnostics.finite);
assert(diagnostics.interpolationConstantResidual < 1.0e-12);
end

function [sourceValues,targetValues] = polynomial_block( ...
        sourceR,targetR,nz,nh,offset)
sourceValues = complex(zeros(numel(sourceR)*nz,nh));
targetValues = complex(zeros(numel(targetR)*nz,nh));
for harmonicIndex = 1:nh
    for verticalIndex = 1:nz
        coefficient = offset+verticalIndex+1i*harmonicIndex;
        source = coefficient+sourceR+0.3*sourceR.^2-0.2*sourceR.^3;
        target = coefficient+targetR+0.3*targetR.^2-0.2*targetR.^3;
        sourceRows = (verticalIndex-1)*numel(sourceR)+(1:numel(sourceR));
        targetRows = (verticalIndex-1)*numel(targetR)+(1:numel(targetR));
        sourceValues(sourceRows,harmonicIndex) = source;
        targetValues(targetRows,harmonicIndex) = target;
    end
end
end

function layout = test_cylinder_layout(nr,nzD,nzL)
nD = nr*nzD;
nL = nr*nzL;
next = 0;
layout.d.ur = next+(1:nD); next = next+nD;
layout.d.ut = next+(1:nD); next = next+nD;
layout.d.w = next+(1:nD); next = next+nD;
layout.d.p = next+(1:nD); next = next+nD;
layout.l.ur = next+(1:nL); next = next+nL;
layout.l.ut = next+(1:nL); next = next+nL;
layout.l.w = next+(1:nL); next = next+nL;
layout.l.p = next+(1:nL); next = next+nL;
layout.zeta = next+(1:nr); next = next+nr;
layout.ndof = next;
layout.nr = nr;
layout.nzD = nzD;
layout.nzL = nzL;
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
