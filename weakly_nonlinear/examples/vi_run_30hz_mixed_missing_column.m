function output = vi_run_30hz_mixed_missing_column()
%VI_RUN_30HZ_MIXED_MISSING_COLUMN Recompute the m0-source mixed column.
%
% Retain the validated pair-specific operator and direct/adjoint modes, but
% allow global physical-residual refinement when the cylinder Schur seed
% narrowly misses the strict forced-field gate. Only g11 and g21 are
% requested; g12 and g22 are already strict in the complete V70 run.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
sourceFile = fullfile(root, ...
    'vi_wnl_complete_ag0-3_fHz-30_modes-m0l6-m2l2_Nr16.mat');
outputFile = fullfile(root, ...
    'vi_wnl_targeted_ag0-3_fHz-30_modes-m0l6-m2l2_source1_Nr16.mat');
if ~isfile(sourceFile)
    error('vi_run_30hz_mixed_missing_column:MissingSource', ...
        'Validated mixed recovery is missing: %s',sourceFile);
end

saved = load(sourceFile,'output');
input = saved.output.input;
input.execution.profile = 'final';
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.options.coefficientPairs = [1,1;2,1];
input.options.coefficientModeResidualTolerance = 1.0e-7;
input.options.eigenpairRefinementTolerance = 1.0e-7;
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
input.options.forcedUseBlockGmres = false;
input.options.forcedUseRankAwareMinimumNorm = false;
input.options.forcedCylinderSchurUseNullspaceBordering = false;
input.options.forcedCylinderSchurFactorMethod = 'columnQr';
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
input.run.recoveredModeFile = sourceFile;
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
output = vi_wnl_run_with_input(input);
end
