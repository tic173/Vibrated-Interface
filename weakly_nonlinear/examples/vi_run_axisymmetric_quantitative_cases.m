function outputs = vi_run_axisymmetric_quantitative_cases(stage,cases)
%VI_RUN_AXISYMMETRIC_QUANTITATIVE_CASES Run the two strict m=0 cases.
%
% outputs = vi_run_axisymmetric_quantitative_cases('recovery')
% outputs = vi_run_axisymmetric_quantitative_cases('coefficients')
% outputs = vi_run_axisymmetric_quantitative_cases('reproject')
% outputs = vi_run_axisymmetric_quantitative_cases('all')
%
% The recovery stage writes coefficient-ready direct/adjoint modes. The
% coefficient stage requires that cache and accepts cubic coefficients only
% when every complete nonlinear forced-field residual passes 1e-6.

if nargin < 1 || isempty(stage)
    stage = 'recovery';
end
stage = validatestring(stage, ...
    {'recovery','coefficients','reproject','trajectory','all'});
if nargin < 2 || isempty(cases)
    cases = {'ag3_30Hz_m0l6_m0l2','ag4_60Hz_m0l9_m0l2'};
elseif ischar(cases) || isstring(cases)
    cases = cellstr(cases);
end
repositoryRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
definitions = quantitative_case_definitions(repositoryRoot);
outputs = cell(numel(cases),1);
for caseIndex = 1:numel(cases)
    caseName = validatestring(cases{caseIndex},fieldnames(definitions));
    definition = definitions.(caseName);
    switch stage
        case 'recovery'
            outputs{caseIndex} = vi_wnl_run_with_input( ...
                make_override(definition,'recovery'));
        case 'coefficients'
            outputs{caseIndex} = vi_wnl_run_with_input( ...
                make_override(definition,'coefficients'));
        case 'trajectory'
            outputs{caseIndex} = vi_wnl_run_with_input( ...
                make_override(definition,'trajectory'));
        case 'reproject'
            outputs{caseIndex} = ...
                vi_reproject_saved_real_coefficients( ...
                definition.outputFile,5.0e-6,true);
        case 'all'
            vi_wnl_run_with_input(make_override(definition,'recovery'));
            outputs{caseIndex} = vi_wnl_run_with_input( ...
                make_override(definition,'coefficients'));
    end
end
end

function definitions = quantitative_case_definitions(repositoryRoot)
definitions.ag3_30Hz_m0l6_m0l2 = struct( ...
    'frequencyHz',30,'acceleration',3,'radialIndices',[6,2], ...
    'verticalPoints',35,'outputFile',fullfile(repositoryRoot,'weakly_nonlinear','data', ...
    'vi_wnl_quantitative_ag0-3_fHz-30_modes-m0l6-m0l2.mat'));
definitions.ag4_60Hz_m0l9_m0l2 = struct( ...
    'frequencyHz',60,'acceleration',4,'radialIndices',[9,2], ...
    'verticalPoints',31,'outputFile',fullfile(repositoryRoot,'weakly_nonlinear','data', ...
    'vi_wnl_quantitative_ag0-4_fHz-60_modes-m0l9-m0l2.mat'));
end

function input = make_override(definition,stage)
input.execution.profile = 'final';
input.dimensional.frequencyHz = definition.frequencyHz;
input.forcing.analysisAmplitude = definition.acceleration;
input.numberOfModes = 2;
input.modes(1).m = 0;
input.modes(1).radialIndex = definition.radialIndices(1);
input.modes(1).s = 0.5;
input.modes(2).m = 0;
input.modes(2).radialIndex = definition.radialIndices(2);
input.modes(2).s = 0;
input.numerics.N = 11;
input.numerics.Nr = 16;
input.numerics.Ntheta = 32;
input.numerics.verticalGrid.type = 'multidomain';
input.numerics.verticalGrid.pointsPerBoundaryLayer = ...
    definition.verticalPoints;
input.numerics.verticalGrid.pointsInBulk = ...
    definition.verticalPoints+4;
input.execution.startParallelPool = true;
input.execution.parallelWorkers = 2;
input.boundary.contactLine = 'free';
input.boundary.sidewallTangentialCondition = 'stressFree';
input.options.coefficientModeResidualTolerance = 5.0e-7;
input.options.eigenpairRefinementTolerance = 5.0e-7;
input.options.adjointUseModelSolver = true;
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedCylinderSchurUseEquilibratedBlocks = true;
input.options.forcedCylinderSchurFactorMethod = 'columnQr';
input.options.forcedCylinderSchurUsePressureCompatibilityProjection = true;
input.options.forcedCylinderSchurPressureCompatibilityTolerance = 1.0e-5;
input.options.forcedModelSolverOnly = true;
% This gate should not be tighter than the O(1e-6) mode/forced-field
% accuracy used by these refined production cases. The exact m=0
% coefficients are real; this permits at most five parts per million of
% numerical imaginary leakage before projecting onto that real coordinate.
input.options.realCoefficientImaginaryTolerance = 5.0e-6;
input.options.stopOnUnconvergedMode = true;
input.options.stopOnUnconvergedForcedSolve = true;
input.options.forcedFailFast = true;
input.weaklyNonlinear.reference = 'analysisAmplitude';
input.weaklyNonlinear.transientModel = 'cubicEnvelope';
input.weaklyNonlinear.autoSelectSlowModes = false;
input.weaklyNonlinear.requireSlowModesForComparison = false;
input.weaklyNonlinear.allowExploratoryTrajectoryWithUnconvergedForcedFields = ...
    false;
input.initialConditions.amplitudesOverH = [1.0e-3;1.0e-3];
input.initialConditions.phases = [0;0];
input.comparison.endForcingPeriod = 6;
input.comparison.snapshotForcingPeriod = 6;
input.comparison.runGrowthRateSweep = false;
input.run.outputFile = definition.outputFile;
input.run.recoveredModeFile = definition.outputFile;
input.run.saveResults = true;
input.run.saveProfile = 'compact';
input.run.weaklyNonlinear = true;
input.run.compareLinearAndWNL = true;
switch stage
    case 'recovery'
        input.run.postprocessSavedCoefficientsOnly = false;
        input.run.modeRecoveryOnly = true;
        input.run.reuseRecoveredModes = false;
        input.run.requireRecoveredModes = false;
        input.run.allowModeRecoveryDuringCoefficientRun = false;
        input.run.plotLinearTransient = false;
    case 'coefficients'
        input.run.postprocessSavedCoefficientsOnly = false;
        input.run.modeRecoveryOnly = false;
        input.run.reuseRecoveredModes = true;
        input.run.requireRecoveredModes = true;
        input.run.allowModeRecoveryDuringCoefficientRun = false;
    case 'trajectory'
        input.run.postprocessSavedCoefficientsOnly = true;
end
end
