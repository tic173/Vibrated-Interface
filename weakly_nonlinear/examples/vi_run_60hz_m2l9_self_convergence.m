function output = vi_run_60hz_m2l9_self_convergence(stage,nr)
%VI_RUN_60HZ_M2L9_SELF_CONVERGENCE Chebyshev audit of the 60 Hz m2l9 self term.
%
% vi_run_60hz_m2l9_self_convergence('recovery',40)
% output = vi_run_60hz_m2l9_self_convergence('coefficient',40)

if nargin < 1 || isempty(stage)
    stage = 'recovery';
end
stage = validatestring(stage,{'recovery','coefficient'});
if nargin < 2 || isempty(nr)
    nr = 40;
end
validateattributes(nr,{'numeric'},{'scalar','integer','>=',28});

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sourceFile = nearest_two_mode_cache(root,nr);
outputFile = fullfile(root,'weakly_nonlinear','data',sprintf( ...
    'vi_wnl_selfaudit_ag0-4_fHz-60_mode-m2l9_chebyshev_Nr%d.mat',nr));
if strcmp(stage,'coefficient')
    sourceFile = outputFile;
end
if ~isfile(sourceFile)
    error('vi_run_60hz_m2l9_self_convergence:MissingSource', ...
        'Required recovery cache is missing: %s',sourceFile);
end

saved = load(vi_wnl_resolve_data_file(sourceFile),'output');
input = saved.output.input;
input.modes = input.modes(1);
input.numberOfModes = 1;
input.initialConditions.amplitudesOverH = ...
    input.initialConditions.amplitudesOverH(1);
input.initialConditions.phases = input.initialConditions.phases(1);
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.numerics.Nr = nr;
input.numerics.radialGrid.type = 'chebyshev';
input.numerics.radialGrid.referenceModes = [];
input.numerics.radialGrid.preferredBasisLabels = {};
if strcmp(stage,'recovery')
    [targetR,~,~] = vi_chebyshev_lobatto( ...
        nr,[0,saved.output.parameters.R0]);
    [input.modes.directSeed,seedDiagnostics] = ...
        vi_prolong_cylinder_mode_radially( ...
        saved.output.weaklyNonlinear.modes{1}, ...
        saved.output.operatorMetadata,targetR);
    input.modes.directSeedDiagnostics = seedDiagnostics;
end
input.options.coefficientPairs = [];
input.options.coefficientModeResidualTolerance = 1.0e-7;
input.options.eigenpairRefinementTolerance = 1.0e-7;
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
input.options.forcedCylinderSchurUseNullspaceBordering = false;
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
input.run.modeRecoveryOnly = strcmp(stage,'recovery');
input.run.reuseRecoveredModes = strcmp(stage,'coefficient');
input.run.requireRecoveredModes = strcmp(stage,'coefficient');
output = vi_wnl_run_with_input(input);
end

function filename = nearest_two_mode_cache(root,nr)
candidates = dir(fullfile(root,'weakly_nonlinear','data', ...
    ['vi_wnl_complete_ag0-4_fHz-60_modes-m2l9-m2l2_', ...
     'chebyshev_Nr*.mat']));
bestNr = -Inf;
filename = '';
for candidateIndex = 1:numel(candidates)
    token = regexp(candidates(candidateIndex).name, ...
        '_Nr(\d+)\.mat$','tokens','once');
    if isempty(token)
        continue;
    end
    candidateNr = str2double(token{1});
    if candidateNr < nr && candidateNr > bestNr
        bestNr = candidateNr;
        filename = fullfile(candidates(candidateIndex).folder, ...
            candidates(candidateIndex).name);
    end
end
end
