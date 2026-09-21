function output = vi_run_60hz_m2_radial_refinement(stage,nr,gridType)
%VI_RUN_60HZ_M2_RADIAL_REFINEMENT Refine the 60 Hz m2l9-m2l2 audit.
%
% vi_run_60hz_m2_radial_refinement('recovery',18)
% output = vi_run_60hz_m2_radial_refinement('self',18)
% output = vi_run_60hz_m2_radial_refinement('full',18)
% output = vi_run_60hz_m2_radial_refinement('self',28,'chebyshev')
%
% Both saved Nr=16 full-state modes are interpolated to the requested radial
% space. Recovery is certified before any cubic fields are solved.

if nargin < 1 || isempty(stage)
    stage = 'recovery';
end
stage = validatestring(stage,{'recovery','self','self1','self2','full'});
if nargin < 2 || isempty(nr)
    nr = 18;
end
validateattributes(nr,{'numeric'},{'scalar','integer','>=',17});
if nargin < 3 || isempty(gridType)
    gridType = 'besselEnriched';
end
gridType = validatestring(gridType,{'besselEnriched','chebyshev'});

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
baseSourceFile = fullfile(root, ...
    'vi_wnl_ag0-4_fHz-60_modes-m2l9-m2l2.mat');
if strcmp(gridType,'chebyshev')
    gridTag = '_chebyshev';
else
    gridTag = '';
end
outputFile = fullfile(root,sprintf( ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m2l9-m2l2%s_Nr%d.mat', ...
    gridTag,nr));
if ~isfile(baseSourceFile)
    error('vi_run_60hz_m2_radial_refinement:MissingSource', ...
        'Validated Nr=16 recovery cache is missing: %s',baseSourceFile);
end
if ~strcmp(stage,'recovery') && ~isfile(outputFile)
    error('vi_run_60hz_m2_radial_refinement:MissingRecovery', ...
        ['Run the recovery stage first; its certified cache is ', ...
         'missing: %s'],outputFile);
end

if strcmp(stage,'recovery')
    seedFile = nearest_lower_grid_cache(root,gridTag,nr,baseSourceFile);
    source = load(seedFile,'output');
    saved = source;
else
    saved = load(outputFile,'output');
    source = saved;
end
input = saved.output.input;
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.numerics.Nr = nr;
input.numerics.radialGrid.type = gridType;
if strcmp(gridType,'besselEnriched')
    input.numerics.radialGrid.referenceModes = reference_modes();
    input.numerics.radialGrid.preferredBasisLabels = ...
        source.output.operatorMetadata.radialGrid.selectedBasisLabels;
else
    input.numerics.radialGrid.referenceModes = [];
    input.numerics.radialGrid.preferredBasisLabels = {};
end
if strcmp(stage,'recovery')
    [targetR,~,~] = vi_chebyshev_lobatto( ...
        nr,[0,source.output.parameters.R0]);
    for modeIndex = 1:2
        [input.modes(modeIndex).directSeed,seedDiagnostics] = ...
            vi_prolong_cylinder_mode_radially( ...
            source.output.weaklyNonlinear.modes{modeIndex}, ...
            source.output.operatorMetadata,targetR);
        input.modes(modeIndex).directSeedDiagnostics = seedDiagnostics;
    end
end
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
if strcmp(stage,'self')
    input.options.coefficientPairs = [1,1;2,2];
elseif strcmp(stage,'self1')
    input.options.coefficientPairs = [1,1];
elseif strcmp(stage,'self2')
    input.options.coefficientPairs = [2,2];
else
    input.options.coefficientPairs = [];
end
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

function filename = nearest_lower_grid_cache(root,gridTag,nr,fallback)
pattern = sprintf( ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m2l9-m2l2%s_Nr*.mat', ...
    gridTag);
candidates = dir(fullfile(root,pattern));
bestNr = -Inf;
filename = fallback;
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

function modes = reference_modes()
modes = repmat(struct('m',2,'radialIndex',1),2,1);
modes(1).radialIndex = 9;
modes(2).radialIndex = 2;
end
