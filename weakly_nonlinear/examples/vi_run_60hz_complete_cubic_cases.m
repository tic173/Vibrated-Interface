function outputs = vi_run_60hz_complete_cubic_cases(cases)
%VI_RUN_60HZ_COMPLETE_CUBIC_CASES Recompute complete 60 Hz G and H matrices.
%
% outputs = vi_run_60hz_complete_cubic_cases()
% output = vi_run_60hz_complete_cubic_cases('m0l9_m0l2')
%
% The two cases run sequentially while sharing one two-worker process pool.
% Validated direct/adjoint modes are reused, but every slaved field and
% every ordinary or phase-sensitive cubic coefficient is recomputed with
% the current solver release.

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
        error('vi_run_60hz_complete_cubic_cases:MissingRecovery', ...
            'Validated recovery cache is missing: %s', ...
            definition.sourceFile);
    end
    saved = load(vi_wnl_resolve_data_file(definition.sourceFile),'output');
    input = complete_coefficient_input(saved.output.input,definition);
    fprintf('\nStarting complete 60 Hz coefficient case %d/%d: %s\n', ...
        caseIndex,numel(cases),caseName);
    outputs{caseIndex} = vi_wnl_run_with_input(input);
end

if numel(outputs) == 1
    outputs = outputs{1};
end
end

function definitions = case_definitions()
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
definitions.m0l9_m0l2 = definition(root, ...
    'vi_wnl_quantitative_ag0-4_fHz-60_modes-m0l9-m0l2.mat', ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m0l9-m0l2_Nr16.mat',true);
definitions.m2l9_m2l2 = definition(root, ...
    ['vi_wnl_complete_ag0-4_fHz-60_modes-m2l9-m2l2_', ...
     'chebyshev_Nr36.mat'], ...
    ['vi_wnl_complete_ag0-4_fHz-60_modes-m2l9-m2l2_', ...
     'chebyshev_Nr36.mat'],false);
end

function value = definition(root,sourceName,outputName,useNullspaceBordering)
value.sourceFile = fullfile(root,'weakly_nonlinear','data',sourceName);
value.outputFile = fullfile(root,'weakly_nonlinear','data',outputName);
value.useNullspaceBordering = useNullspaceBordering;
end

function input = complete_coefficient_input(input,definition)
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.options.coefficientPairs = [];
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
% For the m=2 pair, its mixed m=0 half-period block has a left nullspace
% that changes under the Floquet mass shift. Column QR resolves each
% temporal diagonal in its own range instead of reusing a fixed border.
input.options.forcedCylinderSchurUseNullspaceBordering = ...
    definition.useNullspaceBordering;
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
