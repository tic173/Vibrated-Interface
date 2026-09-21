function outputs = vi_run_30hz_complete_cubic_cases(cases)
%VI_RUN_30HZ_COMPLETE_CUBIC_CASES Recompute complete 30 Hz G and H matrices.
%
% outputs = vi_run_30hz_complete_cubic_cases()
% output = vi_run_30hz_complete_cubic_cases('m0l6_m2l2')
%
% Each case reuses only validated direct/adjoint Floquet modes. All slaved
% fields, ordinary cubic coefficients G, and phase-sensitive coefficients H
% are recomputed with the current solver release.

definitions = case_definitions();
if nargin < 1 || isempty(cases)
    cases = fieldnames(definitions);
elseif ischar(cases) || isstring(cases)
    cases = cellstr(cases);
end

outputs = cell(numel(cases),1);
for caseIndex = 1:numel(cases)
    caseName = validatestring(cases{caseIndex},fieldnames(definitions));
    definition = definitions.(caseName);
    if ~isfile(definition.sourceFile)
        error('vi_run_30hz_complete_cubic_cases:MissingRecovery', ...
            'Validated recovery cache is missing: %s', ...
            definition.sourceFile);
    end
    saved = load(definition.sourceFile,'output');
    input = complete_coefficient_input(saved.output.input,definition);
    outputs{caseIndex} = vi_wnl_run_with_input(input);
end

if numel(outputs) == 1
    outputs = outputs{1};
end
end

function definitions = case_definitions()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
definitions.m0l6_m0l2 = definition(root,'m0l6-m0l2', ...
    'vi_wnl_canonicalself_ag0-3_fHz-30_modes-m0l6-m0l2_Nr16.mat');
definitions.m0l6_m2l2 = definition(root,'m0l6-m2l2', ...
    'vi_wnl_ag0-3_fHz-30_modes-m0l6-m2l2.mat');
definitions.m2l6_m2l2 = definition(root,'m2l6-m2l2', ...
    'vi_wnl_canonicalself_ag0-3_fHz-30_modes-m2l6-m2l2_Nr16.mat');
end

function value = definition(root,modeToken,sourceName)
value.sourceFile = fullfile(root,sourceName);
value.outputFile = fullfile(root,sprintf( ...
    'vi_wnl_complete_ag0-3_fHz-30_modes-%s_Nr16.mat',modeToken));
end

function input = complete_coefficient_input(input,definition)
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.options.coefficientPairs = [];
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
input.options.stopOnUnconvergedMode = true;
input.options.stopOnUnconvergedForcedSolve = true;
input.options.forcedFailFast = true;
input.options.realCoefficientImaginaryTolerance = 5.0e-6;
input.weaklyNonlinear.reference = 'analysisAmplitude';
input.weaklyNonlinear.transientModel = 'cubicEnvelope';
input.weaklyNonlinear.autoSelectSlowModes = false;
input.weaklyNonlinear.requireSlowModesForComparison = false;
input.comparison.runGrowthRateSweep = false;
input.comparison.plotInterfaceDynamics = false;
input.comparison.plotRadialModeCorrections = false;
input.comparison.plotSecondModeInitialAmplitudeSweep = false;
input.run.outputFile = definition.outputFile;
input.run.recoveredModeFile = definition.sourceFile;
input.run.saveResults = true;
input.run.saveProfile = 'compact';
input.run.weaklyNonlinear = true;
input.run.linearPreview = true;
input.run.plotLinearTransient = false;
input.run.plotHarmonicNorms = false;
input.run.compareLinearAndWNL = false;
input.run.postprocessSavedCoefficientsOnly = false;
input.run.modeRecoveryOnly = false;
input.run.reuseRecoveredModes = true;
input.run.requireRecoveredModes = true;
input.run.allowModeRecoveryDuringCoefficientRun = false;
end
