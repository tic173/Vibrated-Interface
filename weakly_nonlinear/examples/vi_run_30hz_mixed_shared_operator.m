function output = vi_run_30hz_mixed_shared_operator(stage,nr)
%VI_RUN_30HZ_MIXED_SHARED_OPERATOR Mixed m0l6-m2l2 shared radial audit.
%
% vi_run_30hz_mixed_shared_operator('recovery')
% output = vi_run_30hz_mixed_shared_operator('full')
%
% The radial space includes the retained m0l6 and m2l2 branches plus the
% m2l6 branch needed to preserve the validated nonaxisymmetric family
% operator. Its previously selected completion columns are preferred.

if nargin < 1 || isempty(stage)
    stage = 'recovery';
end
stage = validatestring(stage,{'recovery','full'});
if nargin < 2 || isempty(nr)
    nr = 16;
end
validateattributes(nr,{'numeric'},{'scalar','integer','>=',11});

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
templateFile = fullfile(root,'weakly_nonlinear','data', ...
    'vi_wnl_ag0-3_fHz-30_modes-m0l6-m2l2.mat');
outputFile = fullfile(root,'weakly_nonlinear','data',sprintf( ...
    'vi_wnl_mixedshared_ag0-3_fHz-30_modes-m0l6-m2l2_Nr%d.mat',nr));
canonicalM2File = fullfile(root,'weakly_nonlinear','data', ...
    'vi_wnl_canonicalself_ag0-3_fHz-30_modes-m2l6-m2l2_Nr16.mat');
if strcmp(stage,'recovery')
    sourceFile = templateFile;
else
    sourceFile = outputFile;
end
if ~isfile(sourceFile)
    error('vi_run_30hz_mixed_shared_operator:MissingSource', ...
        'Required source file is missing: %s',sourceFile);
end
if ~isfile(canonicalM2File)
    error('vi_run_30hz_mixed_shared_operator:MissingCanonicalBasis', ...
        'Canonical m=2 radial basis is missing: %s',canonicalM2File);
end

saved = load(vi_wnl_resolve_data_file(sourceFile),'output');
canonicalM2 = load(vi_wnl_resolve_data_file(canonicalM2File),'output');
preferredBasisLabels = canonicalM2.output.operatorMetadata.radialGrid. ...
    selectedBasisLabels;
input = saved.output.input;
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.numerics.Nr = nr;
input.numerics.radialGrid.referenceModes = shared_reference_modes();
input.numerics.radialGrid.preferredBasisLabels = preferredBasisLabels;
if strcmp(stage,'recovery') && nr ~= ...
        canonicalM2.output.operatorMetadata.layout.nr
    [targetR,~,~] = vi_chebyshev_lobatto( ...
        nr,[0,canonicalM2.output.parameters.R0]);
    [input.modes(2).directSeed,seedDiagnostics] = ...
        vi_prolong_cylinder_mode_radially( ...
        canonicalM2.output.weaklyNonlinear.modes{2}, ...
        canonicalM2.output.operatorMetadata,targetR);
    input.modes(2).directSeedDiagnostics = seedDiagnostics;
end
input.options.coefficientModeResidualTolerance = 1.0e-7;
input.options.eigenpairRefinementTolerance = 1.0e-7;
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
input.options.stopOnUnconvergedMode = true;
input.options.stopOnUnconvergedForcedSolve = true;
input.options.forcedFailFast = true;
input.options.realCoefficientImaginaryTolerance = 5.0e-6;
input.options.coefficientPairs = [];
input.weaklyNonlinear.reference = 'analysisAmplitude';
input.weaklyNonlinear.transientModel = 'cubicEnvelope';
input.weaklyNonlinear.autoSelectSlowModes = false;
input.weaklyNonlinear.requireSlowModesForComparison = false;
input.comparison.runGrowthRateSweep = false;
input.comparison.plotInterfaceDynamics = false;
input.comparison.plotRadialModeCorrections = false;
input.comparison.plotSecondModeInitialAmplitudeSweep = false;
input.run.outputFile = outputFile;
input.run.recoveredModeFile = outputFile;
input.run.saveResults = true;
input.run.saveProfile = 'compact';
input.run.weaklyNonlinear = true;
input.run.linearPreview = true;
input.run.plotLinearTransient = false;
input.run.plotHarmonicNorms = false;
input.run.compareLinearAndWNL = false;
input.run.postprocessSavedCoefficientsOnly = false;
input.run.allowModeRecoveryDuringCoefficientRun = false;
if strcmp(stage,'recovery')
    input.run.modeRecoveryOnly = true;
    input.run.reuseRecoveredModes = false;
    input.run.requireRecoveredModes = false;
else
    input.run.modeRecoveryOnly = false;
    input.run.reuseRecoveredModes = true;
    input.run.requireRecoveredModes = true;
end
output = vi_wnl_run_with_input(input);
end

function modes = shared_reference_modes()
modes = repmat(struct('m',0,'radialIndex',1),3,1);
modes(1).m = 0; modes(1).radialIndex = 6;
modes(2).m = 2; modes(2).radialIndex = 6;
modes(3).m = 2; modes(3).radialIndex = 2;
end
