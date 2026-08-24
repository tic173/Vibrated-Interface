function result = vi_saved_small_amplitude_transient(source, settings)
%VI_SAVED_SMALL_AMPLITUDE_TRANSIENT Postprocess saved operating-point data.
% Reuses stored lambda_j and g(j,k); it does not rerun mode or forced solves.
% Supports the fully coupled cubic envelope and the first finite-time cubic
% correction selected by input.weaklyNonlinear.transientModel.
% SOURCE may be a MAT filename or the same input structure used by the main
% driver. For an input structure, the cache path and postprocessing settings
% are resolved automatically.

if nargin < 2 || isempty(settings)
    settings = struct();
end
[filename,settings,sourceInput] = resolve_source(source,settings);
saved = load(filename, 'output');
if ~isfield(saved, 'output')
    error('vi_saved_small_amplitude_transient:MissingOutput', ...
        'The MAT file must contain the driver output structure.');
end
output = saved.output;
if ~isempty(sourceInput)
    assert_matching_case(sourceInput,output.input);
end
required = {'input','initialConditions','parameters','linear', ...
    'weaklyNonlinear'};
for fieldIndex = 1:numel(required)
    if ~isfield(output,required{fieldIndex})
        error('vi_saved_small_amplitude_transient:MissingField', ...
            'output.%s is required.',required{fieldIndex});
    end
end
wnl = output.weaklyNonlinear;
if isempty(wnl)
    error('vi_saved_small_amplitude_transient:MissingCoefficients', ...
        'The saved result does not contain weakly nonlinear coefficients.');
end
if isfield(wnl,'linearCoefficients')
    lambda = wnl.linearCoefficients(:);
else
    lambda = wnl.linearCoefficient(:);
end
g = wnl.g;
if any(~isfinite(real(g(:)))) || any(~isfinite(imag(g(:))))
    error('vi_saved_small_amplitude_transient:NonfiniteCoefficients', ...
        'The saved cubic coefficient matrix contains nonfinite values.');
end

amplitudeScale = 1;
if isfield(output,'operatorMetadata') && ...
        isfield(output.operatorMetadata,'zetaOverHPerUnitAmplitude')
    amplitudeScale = ...
        output.operatorMetadata.zetaOverHPerUnitAmplitude;
elseif isfield(output.input,'comparison') && ...
        isfield(output.input.comparison,'zetaOverHPerUnitAmplitude')
    amplitudeScale = ...
        output.input.comparison.zetaOverHPerUnitAmplitude;
end
gInZetaUnits = g/abs(amplitudeScale)^2;
savedA0 = output.initialConditions.complexAmplitudesOverH(:);
A0 = resolve_initial_amplitudes(settings,savedA0,numel(lambda));

defaultEndPeriod = get_nested(output.input, ...
    {'comparison','endForcingPeriod'}, ...
    output.linear.forcingPeriods(end));
endPeriod = get_setting(settings,'endForcingPeriod',defaultEndPeriod);
[indices,~] = vi_comparison_time_indices( ...
    output.linear.forcingPeriods,endPeriod);
timeRequested = output.linear.timeStar(indices);
periodsRequested = output.linear.forcingPeriods(indices);

limits = struct();
limits.maximumAmplitude = get_setting(settings,'maximumAmplitude', ...
    get_nested(output.input,{'comparison','maximumWnlAmplitudeOverH'},0.5));
limits.maximumRelativeCorrection = get_setting(settings, ...
    'maximumRelativeCorrection',get_nested(output.input, ...
    {'comparison','maximumRelativeCubicCorrection'},0.10));
limits.maximumNonlinearToLinearRateRatio = get_setting(settings, ...
    'maximumNonlinearToLinearRateRatio',get_nested(output.input, ...
    {'comparison','maximumNonlinearToLinearRateRatio'},0.10));
limits.linearRateFloor = 1.0e-12*max(1,output.parameters.omegaStar);

savedTransientModel = get_nested(output.input, ...
    {'weaklyNonlinear','transientModel'},'smallAmplitudeCorrection');
transientModel = validatestring(get_setting(settings, ...
    'transientModel',savedTransientModel), ...
    {'smallAmplitudeCorrection','cubicEnvelope'});
transient = build_saved_transient(transientModel,timeRequested,A0, ...
    lambda,gInZetaUnits,limits);
numberOfTimes = numel(transient.time);
periods = periodsRequested(1:numberOfTimes);
modeLabels = saved_mode_labels(wnl,numel(A0));

strictCoefficients = false;
if isfield(wnl,'forcedSolvesValid')
    strictCoefficients = all(wnl.forcedSolvesValid(:));
end
exploratoryCoefficients = strictCoefficients;
if isfield(wnl,'forcedSolvesExploratoryUsable')
    exploratoryCoefficients = ...
        all(wnl.forcedSolvesExploratoryUsable(:));
end
if ~strictCoefficients
    if exploratoryCoefficients
        warning('vi_saved_small_amplitude_transient:ExploratoryCoefficients', ...
            ['The saved forced fields passed only the exploratory residual ', ...
             'gate. The transient uses coefficients that are not strictly ', ...
             'converged.']);
    else
        warning('vi_saved_small_amplitude_transient:UnconvergedCoefficients', ...
            ['At least one saved forced field missed its residual gate. ', ...
             'Treat the resulting nonlinear trajectory as an unconverged ', ...
             'diagnostic.']);
    end
end

result = transient;
result.sourceFile = filename;
result.forcingPeriods = periods;
result.modeLabels = modeLabels;
result.lambda = lambda;
result.gInZetaUnits = gInZetaUnits;
result.initialAmplitudesOverH = A0;
result.savedInitialAmplitudesOverH = savedA0;
result.initialConditionsOverridden = any(A0 ~= savedA0);
result.strictCoefficients = strictCoefficients;
result.exploratoryCoefficients = exploratoryCoefficients;
plotEnabled = get_setting(settings,'plot',true);

result.interfaceDynamics = [];
result.linearProbeDisplacementOverH = [];
result.nonlinearProbeDisplacementOverH = [];
result.probeDisplacementDifferenceLinearMinusNonlinear = [];
% Backward-compatible aliases retained for existing postprocessing scripts.
result.cubicCorrectedProbeDisplacementOverH = [];
result.probeDisplacementDifferenceLinearMinusCubic = [];
result.secondModeInitialAmplitudeSweep = [];
if isfield(output,'operatorMetadata') && ...
        isstruct(output.operatorMetadata)
    usedIndices = indices(1:numberOfTimes);
    linearInterfaceResult = truncate_linear_result( ...
        output.linear,usedIndices);
    interfaceSettings = output.input.comparison;
    if ~isempty(sourceInput) && isfield(sourceInput,'comparison')
        interfaceSettings = overlay_struct( ...
            interfaceSettings,sourceInput.comparison);
    end
    interfaceSettings.plotInterfaceDynamics = plotEnabled && ...
        get_nested(interfaceSettings,{'plotInterfaceDynamics'},true);
    interfaceSettings.nonlinearModelLabel = transient.modelLabel;
    interfaceDynamics = vi_compare_interface_dynamics_modes( ...
        linearInterfaceResult,transient.nonlinearAmplitudes,wnl, ...
        output.operatorMetadata,output.parameters,interfaceSettings, ...
        amplitudeScale,transient.linearAmplitudes);
    result.interfaceDynamics = interfaceDynamics;
    result.linearProbeDisplacementOverH = ...
        interfaceDynamics.linearProbeSignal;
    result.nonlinearProbeDisplacementOverH = ...
        interfaceDynamics.wnlTotalProbeSignal;
    result.probeDisplacementDifferenceLinearMinusNonlinear = ...
        interfaceDynamics.probeDifferenceLinearMinusWnl;
    result.cubicCorrectedProbeDisplacementOverH = ...
        result.nonlinearProbeDisplacementOverH;
    result.probeDisplacementDifferenceLinearMinusCubic = ...
        result.probeDisplacementDifferenceLinearMinusNonlinear;
    sweepAmplitudes = get_setting(settings, ...
        'secondModeInitialAmplitudeSweepOverH',[]);
    if numel(A0) >= 2 && ~isempty(sweepAmplitudes)
        result.secondModeInitialAmplitudeSweep = ...
            build_second_mode_probe_sweep( ...
            sweepAmplitudes,timeRequested,periodsRequested,indices,A0, ...
            lambda,gInZetaUnits,limits,output,wnl,amplitudeScale, ...
            interfaceSettings,modeLabels,transientModel);
    end
else
    warning('vi_saved_small_amplitude_transient:NoInterfaceMetadata', ...
        ['The saved result has no operator metadata, so modal amplitudes ', ...
         'were compared without reconstructing probe displacement.']);
end

fprintf('\nSaved-coefficient transient postprocessor\n');
fprintf('  source                    = %s\n',filename);
fprintf('  transient model           = %s\n',transient.modelLabel);
fprintf('  returned end period       = %.6g T_f\n',periods(end));
fprintf('  completed requested window= %d\n', ...
    transient.completedRequestedWindow);
for modeIndex = 1:numel(A0)
    fprintf(['  mode %d max trajectory difference/rate ratio ', ...
        '= %.6g / %.6g\n'],modeIndex, ...
        max(transient.relativeTrajectoryDifferences(:,modeIndex)), ...
        max(transient.nonlinearToLinearRateRatios(:,modeIndex)));
end
if ~transient.completedRequestedWindow
    fprintf('  stop reason               = %s\n', ...
        strjoin(transient.stopReasons,', '));
end
if ~isempty(result.interfaceDynamics)
    fprintf('  probe r/R, theta          = %.6g, %.6g rad\n', ...
        result.interfaceDynamics.probeRadiusOverR, ...
        result.interfaceDynamics.probeTheta);
    fprintf('  probe displacement RMSE   = %.6g [zeta/h]\n', ...
        result.interfaceDynamics.probeRmse);
end

if plotEnabled
    plot_transient(result);
    plot_landau_terms(result);
    if get_setting(settings,'plotSecondModeInitialAmplitudeSweep',true) && ...
            ~isempty(result.secondModeInitialAmplitudeSweep)
        plot_second_mode_probe_sweep( ...
            result.secondModeInitialAmplitudeSweep);
    end
end
end

function transient = build_saved_transient(transientModel,timeRequested,A0, ...
        lambda,g,limits)
switch transientModel
    case 'smallAmplitudeCorrection'
        transient = vi_cubic_transient_correction( ...
            timeRequested,A0,lambda,g,limits);
        transient.modelLabel = 'small-amplitude cubic correction';
        transient.nonlinearAmplitudes = transient.correctedAmplitudes;
        transient.trajectoryDifferenceLinearMinusNonlinear = ...
            transient.linearAmplitudes-transient.nonlinearAmplitudes;
        transient.relativeTrajectoryDifferences = ...
            transient.relativeCorrections;
        termState = transient.linearAmplitudes;
        transient.termEvaluationState = 'exact linear trajectory';
    case 'cubicEnvelope'
        integration = vi_integrate_landau_limited( ...
            timeRequested,A0,lambda,g,limits.maximumAmplitude);
        transient = integration;
        elapsedTime = transient.time-transient.time(1);
        transient.linearAmplitudes = ...
            exp(elapsedTime*lambda.') .* A0.';
        transient.nonlinearAmplitudes = transient.amplitudes;
        transient.correctedAmplitudes = transient.nonlinearAmplitudes;
        transient.trajectoryDifferenceLinearMinusNonlinear = ...
            transient.linearAmplitudes-transient.nonlinearAmplitudes;
        transient.relativeTrajectoryDifferences = abs( ...
            transient.trajectoryDifferenceLinearMinusNonlinear)./ ...
            max(abs(transient.linearAmplitudes),1.0e-14);
        transient.model = 'fully coupled cubic-envelope equation';
        transient.modelLabel = 'coupled cubic envelope';
        transient.stopReasons = {};
        if transient.stoppedAtAmplitudeLimit
            transient.stopReasons = {'amplitude limit'};
        end
        termState = transient.nonlinearAmplitudes;
        transient.termEvaluationState = 'coupled nonlinear trajectory';
end
transient.transientModel = transientModel;
[linearTerms,nonlinearTerms,crossTerms,rateCoefficients,rateRatios] = ...
    landau_term_diagnostics(termState,lambda,g,limits.linearRateFloor);
transient.linearRhsTerms = linearTerms;
transient.nonlinearRhsTerms = nonlinearTerms;
transient.totalRhsTerms = linearTerms+nonlinearTerms;
transient.crossCubicRhsTerms = crossTerms;
transient.nonlinearRateCoefficients = rateCoefficients;
transient.nonlinearToLinearRateRatios = rateRatios;
end

function [linearTerms,nonlinearTerms,crossTerms,rateCoefficients,ratios] = ...
        landau_term_diagnostics(amplitudes,lambda,g,linearRateFloor)
numberOfModes = numel(lambda);
linearTerms = amplitudes.*lambda.';
rateCoefficients = abs(amplitudes).^2*g.';
nonlinearTerms = amplitudes.*rateCoefficients;
crossTerms = complex(zeros(size(amplitudes,1),numberOfModes,numberOfModes));
for targetMode = 1:numberOfModes
    for sourceMode = 1:numberOfModes
        crossTerms(:,targetMode,sourceMode) = amplitudes(:,targetMode).* ...
            g(targetMode,sourceMode).*abs(amplitudes(:,sourceMode)).^2;
    end
end
ratios = abs(rateCoefficients)./max(abs(lambda.'),linearRateFloor);
end

function plot_transient(result)
numberOfModes = numel(result.lambda);
colors = lines(numberOfModes);
figure('Name','Saved-coefficient linear and nonlinear trajectories', ...
    'Position',[100,100,1100,720]);
tiledlayout(2,2,'TileSpacing','compact','Padding','loose');

nexttile;
hold on;
handles = gobjects(2*numberOfModes,1);
labels = cell(2*numberOfModes,1);
for modeIndex = 1:numberOfModes
    handles(2*modeIndex-1) = plot(result.forcingPeriods, ...
        real(result.linearAmplitudes(:,modeIndex)),'-', ...
        'Color',colors(modeIndex,:),'LineWidth',1.3);
    handles(2*modeIndex) = plot(result.forcingPeriods, ...
        real(result.nonlinearAmplitudes(:,modeIndex)),'--', ...
        'Color',colors(modeIndex,:),'LineWidth',1.5);
    labels{2*modeIndex-1} = sprintf('linear: %s', ...
        result.modeLabels{modeIndex});
    labels{2*modeIndex} = sprintf('%s: %s',result.modelLabel, ...
        result.modeLabels{modeIndex});
end
yline(0,':','HandleVisibility','off'); grid on;
xlabel('forcing periods, t/T_f'); ylabel('Re(A_j), \zeta/h');
legend(handles,labels,'Interpreter','none','Location','best');
title('corresponding signed modal trajectories');

nexttile;
plot(result.forcingPeriods,real( ...
    result.trajectoryDifferenceLinearMinusNonlinear),'LineWidth',1.3);
yline(0,':','HandleVisibility','off'); grid on;
xlabel('forcing periods, t/T_f'); ylabel('Re(A_{L,j}-A_{NL,j}), \zeta/h');
legend(result.modeLabels,'Interpreter','none','Location','best');
title('linear minus nonlinear trajectory');

nexttile;
plot(result.forcingPeriods,result.relativeTrajectoryDifferences, ...
    'LineWidth',1.3); grid on;
xlabel('forcing periods, t/T_f');
ylabel('|A_{L,j}-A_{NL,j}|/|A_{L,j}|');
legend(result.modeLabels,'Interpreter','none','Location','best');
title('relative trajectory difference');

nexttile;
plot(result.forcingPeriods,result.nonlinearToLinearRateRatios, ...
    'LineWidth',1.3); grid on;
xlabel('forcing periods, t/T_f');
ylabel('|N_j(A)|/|\lambda_j A_j|');
legend(result.modeLabels,'Interpreter','none','Location','best');
title('nonlinear-to-linear RHS ratio');

if ~isempty(result.interfaceDynamics)
    figure('Name','Probe displacement comparison', ...
        'Position',[140,140,900,700]);
    tiledlayout(2,1,'TileSpacing','compact','Padding','loose');
    nexttile;
    plot(result.forcingPeriods, ...
        result.linearProbeDisplacementOverH,'LineWidth',1.3);
    hold on;
    plot(result.forcingPeriods, ...
        result.nonlinearProbeDisplacementOverH,'--','LineWidth',1.5);
    yline(0,':','HandleVisibility','off'); grid on;
    xlabel('forcing periods, t/T_f'); ylabel('\zeta/h at probe');
    legend('corresponding linear',result.modelLabel,'Location','best');
    title(sprintf('probe r/R=%.4g, \\theta=%.4g', ...
        result.interfaceDynamics.probeRadiusOverR, ...
        result.interfaceDynamics.probeTheta));

    nexttile;
    plot(result.forcingPeriods, ...
        result.probeDisplacementDifferenceLinearMinusNonlinear, ...
        'LineWidth',1.3);
    hold on; yline(0,':','HandleVisibility','off'); grid on;
    xlabel('forcing periods, t/T_f');
    ylabel('(linear-nonlinear) \zeta/h');
    title(sprintf('corresponding probe difference, RMSE=%.4g', ...
        result.interfaceDynamics.probeRmse));
end
end

function plot_landau_terms(result)
numberOfModes = numel(result.lambda);
figure('Name','Corresponding linear and nonlinear modal terms', ...
    'Position',[180,120,1100,max(420,340*numberOfModes)]);
layout = tiledlayout(numberOfModes,2, ...
    'TileSpacing','compact','Padding','loose');
for modeIndex = 1:numberOfModes
    nexttile(2*modeIndex-1);
    plot(result.forcingPeriods,real(result.linearRhsTerms(:,modeIndex)), ...
        'LineWidth',1.3);
    hold on;
    plot(result.forcingPeriods,real(result.nonlinearRhsTerms(:,modeIndex)), ...
        '--','LineWidth',1.4);
    yline(0,':','HandleVisibility','off'); grid on;
    xlabel('forcing periods, t/T_f'); ylabel('Re(dA_j/dt^*)');
    legend('\lambda_j A_j','N_j(A)','Location','best');
    title(sprintf('signed RHS terms: %s',result.modeLabels{modeIndex}), ...
        'Interpreter','none');

    nexttile(2*modeIndex);
    plot(result.forcingPeriods,abs(result.linearRhsTerms(:,modeIndex)), ...
        'LineWidth',1.3);
    hold on;
    plot(result.forcingPeriods,abs(result.nonlinearRhsTerms(:,modeIndex)), ...
        '--','LineWidth',1.4);
    grid on;
    xlabel('forcing periods, t/T_f'); ylabel('|dA_j/dt^*|');
    legend('|\lambda_j A_j|','|N_j(A)|','Location','best');
    title(sprintf('RHS magnitudes: mode %d',modeIndex));
end
sgtitle(layout,sprintf('terms evaluated on the %s', ...
    result.termEvaluationState));
end

function sweep = build_second_mode_probe_sweep(amplitudes,timeRequested, ...
        periodsRequested,indices,A0,lambda,gInZetaUnits,limits,output,wnl, ...
        amplitudeScale,interfaceSettings,modeLabels,transientModel)
validateattributes(amplitudes,{'numeric'}, ...
    {'vector','real','nonnegative','finite','nonempty'});
amplitudes = unique(amplitudes(:).','stable');
numberOfCases = numel(amplitudes);
numberOfTimes = numel(timeRequested);
numberOfModes = numel(A0);
linearSignals = NaN(numberOfTimes,numberOfModes,numberOfCases);
nonlinearSignals = NaN(numberOfTimes,numberOfModes,numberOfCases);
controlledEndPeriods = NaN(numberOfCases,1);
completedRequestedWindow = false(numberOfCases,1);
stopReasons = cell(numberOfCases,1);
phase2 = angle(A0(2));
sweepInterfaceSettings = interfaceSettings;
sweepInterfaceSettings.plotInterfaceDynamics = false;
sweepInterfaceSettings.verbose = false;

for caseIndex = 1:numberOfCases
    caseA0 = A0;
    caseA0(2) = amplitudes(caseIndex)*exp(1i*phase2);
    caseTransient = build_saved_transient( ...
        transientModel,timeRequested,caseA0,lambda,gInZetaUnits,limits);
    caseTimeCount = numel(caseTransient.time);
    caseLinearResult = truncate_linear_result( ...
        output.linear,indices(1:caseTimeCount));
    caseDynamics = vi_compare_interface_dynamics_modes( ...
        caseLinearResult,caseTransient.nonlinearAmplitudes,wnl, ...
        output.operatorMetadata,output.parameters,sweepInterfaceSettings, ...
        amplitudeScale,caseTransient.linearAmplitudes);
    linearSignals(1:caseTimeCount,:,caseIndex) = ...
        caseDynamics.linearModeProbeSignals;
    nonlinearSignals(1:caseTimeCount,:,caseIndex) = ...
        caseDynamics.wnlPrimaryModeProbeSignals;
    controlledEndPeriods(caseIndex) = ...
        periodsRequested(caseTimeCount);
    completedRequestedWindow(caseIndex) = ...
        caseTransient.completedRequestedWindow;
    stopReasons{caseIndex} = caseTransient.stopReasons;
end

sweep = struct();
sweep.secondModeInitialAmplitudesOverH = amplitudes(:);
sweep.fixedFirstModeInitialAmplitudeOverH = A0(1);
sweep.secondModeInitialPhase = phase2;
sweep.forcingPeriods = periodsRequested(:);
sweep.modeLabels = modeLabels;
sweep.transientModel = transientModel;
sweep.nonlinearModelLabel = caseTransient.modelLabel;
sweep.linearPrimaryModeProbeSignals = linearSignals;
sweep.nonlinearPrimaryModeProbeSignals = nonlinearSignals;
sweep.linearMinusNonlinearPrimaryModeProbeSignals = ...
    linearSignals-nonlinearSignals;
% Backward-compatible aliases.
sweep.cubicCorrectedPrimaryModeProbeSignals = nonlinearSignals;
sweep.linearMinusCubicPrimaryModeProbeSignals = ...
    sweep.linearMinusNonlinearPrimaryModeProbeSignals;
sweep.controlledEndForcingPeriods = controlledEndPeriods;
sweep.completedRequestedWindow = completedRequestedWindow;
sweep.stopReasons = stopReasons;
sweep.probeRadiusOverR = caseDynamics.probeRadiusOverR;
sweep.probeTheta = caseDynamics.probeTheta;
end

function plot_second_mode_probe_sweep(sweep)
amplitudes = sweep.secondModeInitialAmplitudesOverH;
numberOfCases = numel(amplitudes);
colors = lines(numberOfCases);
labels = arrayfun(@(value) sprintf('A_2(0) = %.0e',value), ...
    amplitudes,'UniformOutput',false);
periods = sweep.forcingPeriods;

figure('Name','Second-mode initial-amplitude probe sweep', ...
    'Position',[100,100,1200,760]);
layout = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');
for modeIndex = 1:2
    nexttile(modeIndex);
    hold on;
    amplitudeHandles = gobjects(numberOfCases,1);
    for caseIndex = 1:numberOfCases
        plot(periods,sweep.linearPrimaryModeProbeSignals( ...
            :,modeIndex,caseIndex),'--','Color',colors(caseIndex,:), ...
            'LineWidth',1.15,'HandleVisibility','off');
        amplitudeHandles(caseIndex) = plot(periods, ...
            sweep.nonlinearPrimaryModeProbeSignals( ...
            :,modeIndex,caseIndex),'Color',colors(caseIndex,:), ...
            'LineWidth',1.35);
    end
    yline(0,':','HandleVisibility','off'); grid on;
    currentAxes = gca;
    currentAxes.YAxis.Exponent = 0;
    xlabel('forcing periods, t/T_f'); ylabel('\zeta_j/h at probe');
    title(sprintf('mode %d primary: %s',modeIndex, ...
        sweep.modeLabels{modeIndex}),'Interpreter','none');
    if modeIndex == 1
        nonlinearStyle = plot(NaN,NaN,'k-','LineWidth',1.35);
        linearStyle = plot(NaN,NaN,'k--','LineWidth',1.15);
        legend([amplitudeHandles;nonlinearStyle;linearStyle], ...
            [labels(:);{sweep.nonlinearModelLabel;'pure linear'}], ...
            'Location','best');
    end

    nexttile(modeIndex+2);
    hold on;
    for caseIndex = 1:numberOfCases
        plot(periods,sweep.linearMinusNonlinearPrimaryModeProbeSignals( ...
            :,modeIndex,caseIndex),'Color',colors(caseIndex,:), ...
            'LineWidth',1.35);
    end
    yline(0,':','HandleVisibility','off'); grid on;
    currentAxes = gca;
    currentAxes.YAxis.Exponent = 0;
    xlabel('forcing periods, t/T_f');
    ylabel('(linear-nonlinear) \zeta_j/h');
    title(sprintf('linear minus nonlinear: mode %d',modeIndex));
end
sgtitle(layout,sprintf(['probe r/R=%.4g, \\theta=%.4g; ', ...
    'fixed A_1(0)=%.3g; %s'],sweep.probeRadiusOverR,sweep.probeTheta, ...
    abs(sweep.fixedFirstModeInitialAmplitudeOverH), ...
    sweep.nonlinearModelLabel));
end

function labels = saved_mode_labels(wnl,numberOfModes)
labels = cell(numberOfModes,1);
if isfield(wnl,'modes')
    for modeIndex = 1:numberOfModes
        labels{modeIndex} = wnl.modes{modeIndex}.spec.label;
    end
elseif isfield(wnl,'mode')
    labels{1} = wnl.mode.spec.label;
end
for modeIndex = 1:numberOfModes
    if isempty(labels{modeIndex})
        labels{modeIndex} = sprintf('mode %d',modeIndex);
    end
end
end

function value = get_setting(settings,name,defaultValue)
if isfield(settings,name) && ~isempty(settings.(name))
    value = settings.(name);
else
    value = defaultValue;
end
end

function A0 = resolve_initial_amplitudes(settings,savedA0,numberOfModes)
if isfield(settings,'complexAmplitudesOverH') && ...
        ~isempty(settings.complexAmplitudesOverH)
    A0 = settings.complexAmplitudesOverH(:);
else
    amplitudes = abs(savedA0);
    phases = angle(savedA0);
    if isfield(settings,'amplitudesOverH') && ...
            ~isempty(settings.amplitudesOverH)
        amplitudes = settings.amplitudesOverH(:);
    end
    if isfield(settings,'phases') && ~isempty(settings.phases)
        phases = settings.phases(:);
    end
    validateattributes(amplitudes,{'numeric'}, ...
        {'vector','real','nonnegative','finite','numel',numberOfModes});
    validateattributes(phases,{'numeric'}, ...
        {'vector','real','finite','numel',numberOfModes});
    A0 = amplitudes.*exp(1i*phases);
end
if numel(A0) ~= numberOfModes || any(~isfinite(real(A0))) || ...
        any(~isfinite(imag(A0)))
    error('vi_saved_small_amplitude_transient:InitialConditions', ...
        'The initial amplitudes must contain one finite value per mode.');
end
end

function value = get_nested(source,path,defaultValue)
value = source;
for fieldIndex = 1:numel(path)
    if ~isstruct(value) || ~isfield(value,path{fieldIndex})
        value = defaultValue;
        return;
    end
    value = value.(path{fieldIndex});
end
end

function [filename,settings,sourceInput] = resolve_source(source,settings)
sourceInput = [];
if isstruct(source)
    sourceInput = source;
    settings = inherit_input_settings(settings,sourceInput);
    requestedFile = '';
    if isfield(sourceInput,'run') && ...
            isfield(sourceInput.run,'outputFile') && ...
            ~isempty(sourceInput.run.outputFile) && ...
            ~strcmpi(strtrim(string(sourceInput.run.outputFile)),'auto')
        requestedFile = char(sourceInput.run.outputFile);
    end
    if isempty(requestedFile)
        requestedFile = vi_wnl_output_filename(sourceInput);
    end
    filename = resolve_cache_path(requestedFile);
else
    validateattributes(source, {'char','string'}, {'scalartext'});
    filename = resolve_cache_path(char(source));
end
end

function settings = inherit_input_settings(settings,input)
if isfield(input,'initialConditions')
    if ~isfield(settings,'amplitudesOverH') && ...
            isfield(input.initialConditions,'amplitudesOverH')
        settings.amplitudesOverH = ...
            input.initialConditions.amplitudesOverH;
    end
    if ~isfield(settings,'phases') && ...
            isfield(input.initialConditions,'phases')
        settings.phases = input.initialConditions.phases;
    end
end
if ~isfield(settings,'transientModel') && ...
        isfield(input,'weaklyNonlinear') && ...
        isfield(input.weaklyNonlinear,'transientModel')
    settings.transientModel = input.weaklyNonlinear.transientModel;
end
if ~isfield(input,'comparison')
    return;
end
mapping = { ...
    'endForcingPeriod','endForcingPeriod'; ...
    'maximumWnlAmplitudeOverH','maximumAmplitude'; ...
    'maximumRelativeCubicCorrection','maximumRelativeCorrection'; ...
    'maximumNonlinearToLinearRateRatio', ...
        'maximumNonlinearToLinearRateRatio'; ...
    'secondModeInitialAmplitudeSweepOverH', ...
        'secondModeInitialAmplitudeSweepOverH'; ...
    'plotSecondModeInitialAmplitudeSweep', ...
        'plotSecondModeInitialAmplitudeSweep'};
for mappingIndex = 1:size(mapping,1)
    inputName = mapping{mappingIndex,1};
    settingName = mapping{mappingIndex,2};
    if ~isfield(settings,settingName) && ...
            isfield(input.comparison,inputName)
        settings.(settingName) = input.comparison.(inputName);
    end
end
end

function filename = resolve_cache_path(requestedFile)
moduleRoot = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(moduleRoot);
requestedFile = char(requestedFile);
[~,name,extension] = fileparts(requestedFile);
if isempty(extension)
    extension = '.mat';
    requestedFile = [requestedFile,extension];
end
outputName = [name,extension];
if is_absolute_path(requestedFile)
    % A cache may retain an absolute path from before the project moved.
    % Preserve that path first, then relocate the same file by basename.
    candidates = {requestedFile, ...
        fullfile(repositoryRoot,outputName), ...
        fullfile(repositoryRoot,'weakly_nonlinear','results',outputName), ...
        fullfile(tempdir,outputName)};
else
    candidates = {fullfile(repositoryRoot,requestedFile), ...
        fullfile(repositoryRoot,'weakly_nonlinear','results',outputName), ...
        fullfile(tempdir,outputName)};
end
candidates = unique(candidates,'stable');
for candidateIndex = 1:numel(candidates)
    if isfile(candidates{candidateIndex})
        filename = candidates{candidateIndex};
        return;
    end
end
filename = candidates{1};
error('vi_saved_small_amplitude_transient:CacheNotFound', ...
    'Could not find the automatically selected coefficient cache: %s', ...
    filename);
end

function assert_matching_case(input,savedInput)
required = {'dimensional','forcing','numberOfModes','modes'};
for fieldIndex = 1:numel(required)
    if ~isfield(input,required{fieldIndex}) || ...
            ~isfield(savedInput,required{fieldIndex})
        error('vi_saved_small_amplitude_transient:CaseIdentity', ...
            'Both current and saved inputs must contain %s.', ...
            required{fieldIndex});
    end
end
assert_same_scalar(input.dimensional,savedInput.dimensional, ...
    'frequencyHz','vibration frequency');
assert_same_scalar(input.forcing,savedInput.forcing, ...
    'analysisAmplitude','analysis acceleration');
if input.numberOfModes ~= savedInput.numberOfModes
    error('vi_saved_small_amplitude_transient:CaseMismatch', ...
        'The retained mode count does not match the saved cache.');
end
for modeIndex = 1:input.numberOfModes
    names = {'m','radialIndex','s'};
    for nameIndex = 1:numel(names)
        name = names{nameIndex};
        assert_same_scalar(input.modes(modeIndex), ...
            savedInput.modes(modeIndex),name, ...
            sprintf('mode %d %s',modeIndex,name));
    end
end
end

function assert_same_scalar(current,saved,name,description)
if ~isfield(current,name) || ~isfield(saved,name) || ...
        ~isequaln(current.(name),saved.(name))
    error('vi_saved_small_amplitude_transient:CaseMismatch', ...
        'The current %s does not match the saved cache.',description);
end
end

function tf = is_absolute_path(pathValue)
tf = startsWith(pathValue,'/') || startsWith(pathValue,'\') || ...
    ~isempty(regexp(pathValue,'^[A-Za-z]:[\\/]', 'once'));
end

function linearResult = truncate_linear_result(linearResult,indices)
timeFields = {'timeStar','timeSeconds','forcingPeriods','periodicPart', ...
    'floquetOscillation','displacementOverH','growthEnvelope'};
for fieldIndex = 1:numel(timeFields)
    name = timeFields{fieldIndex};
    if isfield(linearResult,name) && ...
            numel(linearResult.(name)) >= max(indices)
        value = linearResult.(name);
        linearResult.(name) = value(indices);
    end
end
end

function destination = overlay_struct(destination,source)
names = fieldnames(source);
for fieldIndex = 1:numel(names)
    destination.(names{fieldIndex}) = source.(names{fieldIndex});
end
end
