function outputs = vi_run_30hz_common_operator_cases(stage,cases,nr)
%VI_RUN_30HZ_COMMON_OPERATOR_CASES Shared-operator 30 Hz coefficient audit.
%
% outputs = vi_run_30hz_common_operator_cases('recovery')
% outputs = vi_run_30hz_common_operator_cases('self')
% outputs = vi_run_30hz_common_operator_cases('all-self')
% report  = vi_run_30hz_common_operator_cases('verify')
%
% Every retained pair uses one Bessel-enriched radial space built from the
% union m0l6, m0l2, m2l6, and m2l2. The temporal, vertical, azimuthal, wall,
% and pressure-gauge discretizations are also fixed. The self stage requests
% only diagonal cubic coefficients so duplicated modes are independently
% recomputed without paying for cross terms.

if nargin < 1 || isempty(stage)
    stage = 'verify';
end
stage = validatestring(stage, ...
    {'recovery','self','full','all-self','all-full','verify'});
definitions = case_definitions();
if nargin < 2 || isempty(cases)
    cases = fieldnames(definitions);
elseif ischar(cases) || isstring(cases)
    cases = cellstr(cases);
end
if nargin < 3 || isempty(nr)
    nr = 16;
end
validateattributes(nr,{'numeric'}, ...
    {'scalar','integer','>=',11,'finite'});
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
referenceModes = master_reference_modes();

if strcmp(stage,'verify')
    files = cell(numel(cases),1);
    for caseIndex = 1:numel(cases)
        name = validatestring(cases{caseIndex},fieldnames(definitions));
        files{caseIndex} = output_file( ...
            repositoryRoot,definitions.(name),nr);
    end
    outputs = vi_compare_shared_mode_invariance(files);
    reportFile = fullfile(repositoryRoot,sprintf( ...
        'vi_wnl_commonradial_30Hz_Nr%d_invariance_report.mat',nr));
    report = outputs; %#ok<NASGU>
    save(reportFile,'report','-v7');
    fprintf('Saved %s\n',reportFile);
    return;
end

outputs = cell(numel(cases),1);
for caseIndex = 1:numel(cases)
    name = validatestring(cases{caseIndex},fieldnames(definitions));
    definition = definitions.(name);
    if startsWith(stage,'all-')
        vi_wnl_run_with_input(make_override(definition,'recovery', ...
            nr,referenceModes,repositoryRoot));
        coefficientStage = char(extractAfter(stage,'all-'));
        outputs{caseIndex} = vi_wnl_run_with_input(make_override( ...
            definition,coefficientStage,nr,referenceModes, ...
            repositoryRoot));
    else
        outputs{caseIndex} = vi_wnl_run_with_input(make_override( ...
            definition,stage,nr,referenceModes,repositoryRoot));
    end
end
end

function definitions = case_definitions()
definitions.m0l6_m0l2 = struct('m',[0,0],'l',[6,2],'s',[0.5,0]);
definitions.m0l6_m2l2 = struct('m',[0,2],'l',[6,2],'s',[0.5,0]);
definitions.m2l6_m2l2 = struct('m',[2,2],'l',[6,2],'s',[0.5,0]);
end

function modes = master_reference_modes()
modes = repmat(struct('m',0,'radialIndex',1),4,1);
modes(1).m = 0; modes(1).radialIndex = 6;
modes(2).m = 0; modes(2).radialIndex = 2;
modes(3).m = 2; modes(3).radialIndex = 6;
modes(4).m = 2; modes(4).radialIndex = 2;
end

function input = make_override(definition,stage,nr,referenceModes,root)
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.dimensional.frequencyHz = 30;
input.forcing.analysisAmplitude = 3;
input.numberOfModes = 2;
input.modes = repmat(struct('m',0,'radialIndex',1,'s',0),2,1);
for modeIndex = 1:2
    input.modes(modeIndex).m = definition.m(modeIndex);
    input.modes(modeIndex).radialIndex = definition.l(modeIndex);
    input.modes(modeIndex).s = definition.s(modeIndex);
end
input.numerics.N = 11;
input.numerics.Nr = nr;
input.numerics.Ntheta = 32;
input.numerics.radialGrid.type = 'besselEnriched';
input.numerics.radialGrid.referenceModes = referenceModes;
input.numerics.radialGrid.maximumProductOrder = 2;
input.numerics.radialGrid.maximumConditionNumber = 1.0e10;
input.numerics.radialGrid.independenceTolerance = 1.0e-10;
input.numerics.radialGrid.fallbackToChebyshev = false;
input.numerics.verticalGrid.type = 'multidomain';
input.numerics.verticalGrid.pointsPerBoundaryLayer = 11;
input.numerics.verticalGrid.pointsInBulk = 15;
input.boundary.contactLine = 'free';
input.boundary.sidewallTangentialCondition = 'stressFree';
% This mode audit gate remains ten times tighter than the forced-field gate.
% Radial convergence is assessed separately by repeating the common-grid
% calculation at a second Nr, rather than changing the gate between cases.
input.options.coefficientModeResidualTolerance = 1.0e-7;
input.options.eigenpairRefinementTolerance = 1.0e-7;
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
input.options.stopOnUnconvergedMode = true;
input.options.stopOnUnconvergedForcedSolve = true;
input.options.forcedFailFast = true;
input.options.realCoefficientImaginaryTolerance = 5.0e-6;
input.options.coefficientPairs = [1,1;2,2];
input.weaklyNonlinear.reference = 'analysisAmplitude';
input.weaklyNonlinear.transientModel = 'cubicEnvelope';
input.weaklyNonlinear.autoSelectSlowModes = false;
input.weaklyNonlinear.requireSlowModesForComparison = false;
input.initialConditions.amplitudesOverH = [1.0e-3;1.0e-3];
input.initialConditions.phases = [0;0];
input.comparison.runGrowthRateSweep = false;
input.comparison.plotInterfaceDynamics = false;
input.comparison.plotRadialModeCorrections = false;
input.comparison.plotSecondModeInitialAmplitudeSweep = false;
input.run.outputFile = output_file(root,definition,nr);
input.run.recoveredModeFile = input.run.outputFile;
input.run.saveResults = true;
input.run.saveProfile = 'compact';
input.run.weaklyNonlinear = true;
input.run.linearPreview = true;
input.run.plotLinearTransient = false;
input.run.plotHarmonicNorms = false;
input.run.compareLinearAndWNL = false;
switch stage
    case 'recovery'
        input.run.postprocessSavedCoefficientsOnly = false;
        input.run.modeRecoveryOnly = true;
        input.run.reuseRecoveredModes = false;
        input.run.requireRecoveredModes = false;
        input.run.allowModeRecoveryDuringCoefficientRun = false;
    case {'self','full'}
        input.run.postprocessSavedCoefficientsOnly = false;
        input.run.modeRecoveryOnly = false;
        input.run.reuseRecoveredModes = true;
        input.run.requireRecoveredModes = true;
        input.run.allowModeRecoveryDuringCoefficientRun = false;
        if strcmp(stage,'full')
            input.options.coefficientPairs = [];
        end
end
end

function filename = output_file(root,definition,nr)
modeTokens = cell(1,2);
for modeIndex = 1:2
    modeTokens{modeIndex} = sprintf('m%dl%d', ...
        definition.m(modeIndex),definition.l(modeIndex));
end
filename = fullfile(root,sprintf( ...
    'vi_wnl_commonradial_ag0-3_fHz-30_modes-%s_Nr%d.mat', ...
    strjoin(modeTokens,'-'),nr));
end
