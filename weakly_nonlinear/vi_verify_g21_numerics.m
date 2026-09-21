function report = vi_verify_g21_numerics(cacheFile,settings)
%VI_VERIFY_G21_NUMERICS Targeted numerical checks for one saved g21 value.
%
% report = vi_verify_g21_numerics()
% report = vi_verify_g21_numerics(cacheFile)
% report = vi_verify_g21_numerics(cacheFile,settings)
%
% The default 30 Hz, a/g=3 cache uses mode 2 as the target and mode 1 as
% the source, so the computed entry multiplies A_2*A_1^2. Only this cross
% coefficient is recomputed. The default suite changes the pressure-gauge
% row and the nonlinear directional-difference steps without changing the
% physical parameters or modal normalization.
%
% Useful settings are:
%   sweep                       'gauge' (default), 'step', or 'all'
%   gaugeRadialIndices          interior indices; default [4,8,12] for Nr=16
%   stepScales                  default [0.5,1,2]
%   targetIndex, sourceIndex    default 2,1
%   plot                        default true
%   verbose                     default true
%   useParallelNonlinearActions default false
%   coefficientModeResidualTolerance
%                               optional verification-only override
%   reuseSavedBaseline           default true
%
% No coefficient cache or full forced field is saved by this function.

moduleRoot = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(moduleRoot);
if nargin < 1 || isempty(cacheFile)
    cacheFile = fullfile(repositoryRoot, ...
        'vi_wnl_quantitative_ag0-3_fHz-30_modes-m0l6-m0l2.mat');
end
if nargin < 2 || isempty(settings)
    settings = struct();
end
validateattributes(settings,{'struct'},{'scalar'});

cacheFile = resolve_cache_file(cacheFile,repositoryRoot);
loaded = load(cacheFile,'output');
if ~isfield(loaded,'output') || ~isstruct(loaded.output)
    error('vi_verify_g21_numerics:MissingOutput', ...
        'The cache must contain a scalar output structure.');
end
saved = loaded.output;
validate_saved_output(saved);

defaults = default_settings(saved);
settings = vi_wnl_merge_input(defaults,settings);
settings.sweep = validatestring(settings.sweep,{'gauge','step','all'});
validateattributes(settings.targetIndex,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(settings.sourceIndex,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(settings.stepScales,{'numeric'}, ...
    {'vector','real','positive','finite','nonempty'});
settings.stepScales = unique(settings.stepScales(:).','stable');
validateattributes(settings.plot,{'logical'},{'scalar'});
validateattributes(settings.verbose,{'logical'},{'scalar'});
validateattributes(settings.useParallelNonlinearActions, ...
    {'logical'},{'scalar'});
validateattributes(settings.reuseSavedBaseline,{'logical'},{'scalar'});
validateattributes(settings.gaugeRelativeTolerance,{'numeric'}, ...
    {'scalar','real','nonnegative','finite'});
validateattributes(settings.stepRelativeTolerance,{'numeric'}, ...
    {'scalar','real','nonnegative','finite'});
if ~isempty(settings.coefficientModeResidualTolerance)
    validateattributes(settings.coefficientModeResidualTolerance, ...
        {'numeric'},{'scalar','real','positive','finite'});
end

modes = saved.weaklyNonlinear.modes;
numberOfModes = numel(modes);
if max(settings.targetIndex,settings.sourceIndex) > numberOfModes
    error('vi_verify_g21_numerics:ModeIndex', ...
        'The requested target/source index exceeds the saved mode count.');
end
target = modes{settings.targetIndex};
source = modes{settings.sourceIndex};
if target.spec.m ~= 0 || source.spec.m ~= 0
    error('vi_verify_g21_numerics:AxisymmetricOnly', ...
        ['The pressure-gauge verification currently requires two ', ...
         'axisymmetric modes.']);
end

nr = saved.input.numerics.Nr;
if isempty(settings.gaugeRadialIndices)
    settings.gaugeRadialIndices = default_gauge_indices(nr);
end
validateattributes(settings.gaugeRadialIndices,{'numeric'}, ...
    {'vector','integer','>=',2,'<=',nr-1,'finite','nonempty'});
settings.gaugeRadialIndices = unique( ...
    settings.gaugeRadialIndices(:).','stable');

savedGauge = saved_pressure_gauge(saved,nr);
baseQuadraticStep = numeric_field(saved.input.numerics, ...
    'quadraticStep',2.0e-4);
baseCubicStep = numeric_field(saved.input.numerics, ...
    'cubicStep',2.0e-3);
configurations = make_configurations(settings,savedGauge, ...
    baseQuadraticStep,baseCubicStep);

report = struct();
report.cacheFile = cacheFile;
report.targetIndex = settings.targetIndex;
report.sourceIndex = settings.sourceIndex;
report.targetLabel = target.spec.label;
report.sourceLabel = source.spec.label;
report.equationTerm = sprintf('A_%d*A_%d^2', ...
    settings.targetIndex,settings.sourceIndex);
report.savedG21 = saved.weaklyNonlinear.g( ...
    settings.targetIndex,settings.sourceIndex);
report.savedG21Raw = saved_cross_raw(saved,settings);
report.savedGaugeRadialIndex = savedGauge;
report.settings = settings;
report.configurations = repmat(empty_result(),numel(configurations),1);

fprintf('\nTargeted g%d%d verification: %s <- |%s|^2\n', ...
    settings.targetIndex,settings.sourceIndex, ...
    report.targetLabel,report.sourceLabel);
fprintf('  saved coefficient = %.12g %+.12gi\n', ...
    real(report.savedG21),imag(report.savedG21));
fprintf(['  %d unique configuration(s); only the requested cross ', ...
    'coefficient will be evaluated.\n'],numel(configurations));

for configurationIndex = 1:numel(configurations)
    configuration = configurations(configurationIndex);
    fprintf(['\n[%d/%d] gauge ir=%d, quadratic step %.3e, ', ...
        'cubic step %.3e\n'],configurationIndex,numel(configurations), ...
        configuration.gaugeRadialIndex,configuration.quadraticStep, ...
        configuration.cubicStep);
    try
        if settings.reuseSavedBaseline && is_saved_baseline( ...
                configuration,savedGauge,baseQuadraticStep,baseCubicStep)
            value = saved_baseline_result(saved,settings);
            fprintf(['  reused saved baseline g%d%d = %.12g %+.12gi ', ...
                '(no nonlinear solve)\n'],settings.targetIndex, ...
                settings.sourceIndex,real(value.g21Raw), ...
                imag(value.g21Raw));
        else
            value = evaluate_configuration(saved,configuration, ...
                savedGauge,settings);
        end
    catch evaluationError
        value = empty_result();
        value.errorIdentifier = evaluationError.identifier;
        value.errorMessage = evaluationError.message;
        value.wallClockSeconds = NaN;
        warning('vi_verify_g21_numerics:ConfigurationFailed', ...
            'Configuration %d failed: %s',configurationIndex, ...
            evaluationError.message);
    end
    value.gaugeSweepMember = configuration.gaugeSweepMember;
    value.stepSweepMember = configuration.stepSweepMember;
    value.gaugeRadialIndex = configuration.gaugeRadialIndex;
    value.quadraticStep = configuration.quadraticStep;
    value.cubicStep = configuration.cubicStep;
    value.stepScale = configuration.stepScale;
    report.configurations(configurationIndex) = value;
end

report.gaugeTest = assess_group(report.configurations, ...
    'gaugeSweepMember',settings.gaugeRelativeTolerance, ...
    report.savedG21Raw);
report.stepTest = assess_group(report.configurations, ...
    'stepSweepMember',settings.stepRelativeTolerance, ...
    report.savedG21Raw);
print_assessment('pressure-gauge',report.gaugeTest);
print_assessment('directional-step',report.stepTest);
if settings.plot
    report.figureNumber = plot_report(report);
else
    report.figureNumber = [];
end
end

function settings = default_settings(saved)
settings = struct();
settings.sweep = 'gauge';
settings.gaugeRadialIndices = [];
settings.stepScales = [0.5,1,2];
settings.targetIndex = 2;
settings.sourceIndex = 1;
settings.plot = true;
settings.verbose = true;
settings.useParallelNonlinearActions = false;
settings.coefficientModeResidualTolerance = [];
settings.reuseSavedBaseline = true;
settings.gaugeRelativeTolerance = 0.05;
settings.stepRelativeTolerance = 0.05;
if isfield(saved.input,'execution') && ...
        isfield(saved.input.execution,'reportTiming')
    settings.verbose = logical(saved.input.execution.reportTiming);
end
end

function configurations = make_configurations(settings,savedGauge, ...
        quadraticStep,cubicStep)
template = struct('gaugeRadialIndex',savedGauge, ...
    'quadraticStep',quadraticStep,'cubicStep',cubicStep, ...
    'stepScale',1,'gaugeSweepMember',false,'stepSweepMember',false);
configurations = repmat(template,0,1);
if any(strcmp(settings.sweep,{'gauge','all'}))
    for gauge = settings.gaugeRadialIndices
        value = template;
        value.gaugeRadialIndex = gauge;
        value.gaugeSweepMember = true;
        configurations = append_unique(configurations,value);
    end
end
if any(strcmp(settings.sweep,{'step','all'}))
    for scale = settings.stepScales
        value = template;
        value.quadraticStep = scale*quadraticStep;
        value.cubicStep = scale*cubicStep;
        value.stepScale = scale;
        value.stepSweepMember = true;
        configurations = append_unique(configurations,value);
    end
end
end

function configurations = append_unique(configurations,value)
match = find(arrayfun(@(existing) ...
    existing.gaugeRadialIndex == value.gaugeRadialIndex && ...
    existing.quadraticStep == value.quadraticStep && ...
    existing.cubicStep == value.cubicStep,configurations),1);
if isempty(match)
    configurations(end+1,1) = value;
else
    configurations(match).gaugeSweepMember = ...
        configurations(match).gaugeSweepMember || value.gaugeSweepMember;
    configurations(match).stepSweepMember = ...
        configurations(match).stepSweepMember || value.stepSweepMember;
end
end

function value = evaluate_configuration(saved,configuration, ...
        savedGauge,settings)
wallClock = tic;
parameters = saved.parameters;
parameters.numerics = saved.input.numerics;
parameters.numerics.pressureGaugeRadialIndex = ...
    configuration.gaugeRadialIndex;
parameters.numerics.quadraticStep = configuration.quadraticStep;
parameters.numerics.cubicStep = configuration.cubicStep;
if isfield(saved.input,'execution')
    parameters.execution = saved.input.execution;
else
    parameters.execution = struct();
end
parameters.execution.useParallelNonlinearActions = ...
    settings.useParallelNonlinearActions;
parameters.execution.reportTiming = settings.verbose;
if isfield(saved.input,'weaklyNonlinear') && ...
        isfield(saved.input.weaklyNonlinear,'reference') && ...
        strcmpi(saved.input.weaklyNonlinear.reference,'analysisAmplitude')
    parameters.aCritical = parameters.aAnalysis;
    parameters.operatingAmplitude = parameters.aAnalysis;
end

factory = str2func(saved.input.run.operatorFactory);
[operators,metadata] = factory(parameters);
config = struct('omega',parameters.omegaStar, ...
    'N',saved.input.numerics.N,'ndof',metadata.ndof);
model = vi_wnl_model(config,operators);
model.residualLayout = metadata;

specs = cell(numel(saved.weaklyNonlinear.modes),1);
sameGauge = configuration.gaugeRadialIndex == savedGauge;
for modeIndex = 1:numel(specs)
    cachedMode = saved.weaklyNonlinear.modes{modeIndex};
    spec = remove_restart_fields(cachedMode.spec);
    if sameGauge
        spec.direct = cachedMode.vector(:);
        spec.left = cachedMode.left(:);
    else
        spec.directSeed = shift_pressure_gauge(cachedMode.vector, ...
            spec,metadata,configuration.gaugeRadialIndex);
        spec.useDirectSeedAsRefinementState = true;
        spec.directSeedSource = 'saved mode shifted to verification gauge';
    end
    specs{modeIndex} = spec;
end

opts = saved.input.options;
opts.coefficientPairs = [settings.targetIndex,settings.sourceIndex];
opts.reuseTwoModeForcedFields = false;
opts.refineOperatingPointEigenpair = false;
opts.detuning = [];
opts.normalizeDirect = metadata.normalizeDirect;
opts.verbose = settings.verbose;
strictModeTolerance = opts.coefficientModeResidualTolerance;
if ~isempty(settings.coefficientModeResidualTolerance)
    opts.coefficientModeResidualTolerance = ...
        settings.coefficientModeResidualTolerance;
    opts.modeResidualTolerance = max(opts.modeResidualTolerance, ...
        opts.coefficientModeResidualTolerance);
end
result = wnl_analyze_mode_set(model,specs,opts);
entry = result.cross{settings.targetIndex,settings.sourceIndex};
target = result.modes{settings.targetIndex};

value = empty_result();
value.success = entry.validCubicScaling && ...
    isfinite(real(entry.gUnprojected)) && ...
    isfinite(imag(entry.gUnprojected));
value.coefficientAccepted = entry.validCubicScaling;
value.g21 = entry.g;
value.g21Raw = entry.gUnprojected;
value.termProjections = projected_terms(entry,target);
value.termNames = {'mean','mixed','direct cubic'};
value.modeDirectResiduals = cellfun(@(mode) ...
    mode.directResidual,result.modes);
value.modeAdjointResiduals = cellfun(@(mode) ...
    mode.leftResidual,result.modes);
value.modeGateToleranceUsed = opts.coefficientModeResidualTolerance;
value.strictModeGatePassed = max([value.modeDirectResiduals; ...
    value.modeAdjointResiduals]) <= strictModeTolerance;
value.targetConditionIndicator = norm(target.left)* ...
    norm(target.block.Bslow*target.vector)/ ...
    max(abs(target.normalization),eps);
value.forcedRelativeResiduals = [ ...
    forced_relative_residual(entry.qBB), ...
    forced_relative_residual(entry.qAB)];
value.wallClockSeconds = toc(wallClock);
fprintf(['  g%d%d = %.12g %+.12gi; accepted=%d; ', ...
    'mode residual max %.3e; strict mode gate=%d; ', ...
    'forced residual max %.3e; %.1f s\n'], ...
    settings.targetIndex,settings.sourceIndex, ...
    real(value.g21Raw),imag(value.g21Raw), ...
    value.coefficientAccepted, ...
    max([value.modeDirectResiduals;value.modeAdjointResiduals]), ...
    value.strictModeGatePassed, ...
    max(value.forcedRelativeResiduals),value.wallClockSeconds);
end

function spec = remove_restart_fields(spec)
names = {'direct','left','directSeed','directLiftInitialVector', ...
    'recoveryCacheFile'};
present = intersect(fieldnames(spec),names);
if ~isempty(present)
    spec = rmfield(spec,present);
end
end

function vector = shift_pressure_gauge(vector,spec,metadata,gaugeIr)
layout = metadata.layout;
coefficients = reshape(vector,metadata.ndof,numel(spec.n));
gaugeIz = metadata.discretization.denseGaugeVerticalIndex;
gaugeNode = gaugeIr+(gaugeIz-1)*layout.nr;
for harmonicIndex = 1:size(coefficients,2)
    offset = coefficients(layout.d.p(gaugeNode),harmonicIndex);
    coefficients(layout.d.p,harmonicIndex) = ...
        coefficients(layout.d.p,harmonicIndex)-offset;
    coefficients(layout.l.p,harmonicIndex) = ...
        coefficients(layout.l.p,harmonicIndex)-offset;
end
vector = coefficients(:);
end

function projections = projected_terms(entry,target)
names = {'termMean','termMixed','termDirectCubic'};
projections = complex(nan(1,numel(names)));
for index = 1:numel(names)
    if isfield(entry,names{index}) && ...
            ~isempty(entry.(names{index}))
        projections(index) = target.left'*entry.(names{index})(:)/ ...
            target.normalization;
    end
end
end

function value = forced_relative_residual(solution)
value = NaN;
for name = {'relativeEquationResidual','forcingRelativeResidual', ...
        'relativeResidual'}
    if isstruct(solution) && isfield(solution,name{1}) && ...
            isscalar(solution.(name{1}))
        value = solution.(name{1});
        return;
    end
end
end

function value = empty_result()
value = struct('success',false,'coefficientAccepted',false, ...
    'g21',NaN,'g21Raw',NaN,'termNames',{{}}, ...
    'termProjections',complex(nan(1,3)), ...
    'modeDirectResiduals',[],'modeAdjointResiduals',[], ...
    'forcedRelativeResiduals',[],'targetConditionIndicator',NaN, ...
    'modeGateToleranceUsed',NaN,'strictModeGatePassed',false, ...
    'wallClockSeconds',NaN,'errorIdentifier','','errorMessage','', ...
    'gaugeSweepMember',false,'stepSweepMember',false, ...
    'gaugeRadialIndex',NaN,'quadraticStep',NaN,'cubicStep',NaN, ...
    'stepScale',NaN);
end

function assessment = assess_group(configurations,memberName,tolerance, ...
        referenceValue)
selected = [configurations.(memberName)] & ...
    [configurations.success];
values = [configurations(selected).g21Raw];
assessment = struct('numberSuccessful',nnz(selected), ...
    'numberCompared',0,'values',values,'relativeSpread',NaN, ...
    'maximumImaginaryLeakage',NaN, ...
    'maximumRelativeDeviationFromSaved',NaN, ...
    'tolerance',tolerance,'passed',false,'conclusion','not evaluated');
if isempty(values)
    return;
end
comparisonValues = values;
sameAsReference = abs(values-referenceValue) <= ...
    100*eps(max(1,abs(referenceValue)));
if ~any(sameAsReference)
    comparisonValues(end+1) = referenceValue;
end
assessment.numberCompared = numel(comparisonValues);
assessment.values = comparisonValues;
realValues = real(comparisonValues);
assessment.relativeSpread = (max(realValues)-min(realValues))/ ...
    max(abs(median(realValues)),eps);
assessment.maximumRelativeDeviationFromSaved = max( ...
    abs(comparisonValues-referenceValue))/max(abs(referenceValue),eps);
assessment.maximumImaginaryLeakage = max(abs(imag(comparisonValues))./ ...
    max(abs(realValues),eps));
assessment.passed = numel(comparisonValues) >= 2 && ...
    assessment.maximumRelativeDeviationFromSaved <= tolerance;
if numel(comparisonValues) < 2
    assessment.conclusion = 'fewer than two successful configurations';
elseif assessment.passed
    assessment.conclusion = 'stable within the requested tolerance';
else
    assessment.conclusion = 'numerically sensitive';
end
end

function print_assessment(name,value)
fprintf('\n%s test: %s\n',name,value.conclusion);
if value.numberSuccessful > 0
    fprintf(['  successful configurations=%d, compared values=%d, ', ...
        'relative spread=%.3g, ', ...
        'maximum change from saved=%.3g, ', ...
        'maximum relative imaginary leakage=%.3g\n'], ...
        value.numberSuccessful,value.numberCompared,value.relativeSpread, ...
        value.maximumRelativeDeviationFromSaved, ...
        value.maximumImaginaryLeakage);
end
end

function figureNumber = plot_report(report)
figureHandle = figure('Name','g21 numerical verification');
layout = tiledlayout(1,2,'TileSpacing','compact', ...
    'Padding','compact');
title(layout,sprintf('%s <- |%s|^2', ...
    report.targetLabel,report.sourceLabel),'Interpreter','none');

nexttile;
selected = [report.configurations.gaugeSweepMember] & ...
    [report.configurations.success];
plot([report.configurations(selected).gaugeRadialIndex], ...
    real([report.configurations(selected).g21Raw]),'o-', ...
    'LineWidth',1.5,'MarkerSize',7);
yline(real(report.savedG21Raw),'--','saved g_{21}');
grid on;
xlabel('pressure-gauge radial index');
ylabel('Re(g_{21})');
title('gauge invariance');

nexttile;
selected = [report.configurations.stepSweepMember] & ...
    [report.configurations.success];
semilogx([report.configurations(selected).stepScale], ...
    real([report.configurations(selected).g21Raw]),'o-', ...
    'LineWidth',1.5,'MarkerSize',7);
yline(real(report.savedG21Raw),'--','saved g_{21}');
grid on;
xlabel('quadratic and cubic base-step scale');
ylabel('Re(g_{21})');
title('directional-step convergence');
figureNumber = figureHandle.Number;
end

function indices = default_gauge_indices(nr)
indices = unique(max(2,min(nr-1,[ceil(nr/2),round(nr/4), ...
    round(3*nr/4)])),'stable');
end

function gauge = saved_pressure_gauge(saved,nr)
gauge = max(2,min(nr-1,ceil(nr/2)));
if isfield(saved,'operatorMetadata') && ...
        isfield(saved.operatorMetadata,'pressureGauge') && ...
        isfield(saved.operatorMetadata.pressureGauge,'radialIndex')
    gauge = saved.operatorMetadata.pressureGauge.radialIndex;
elseif isfield(saved,'operatorMetadata') && ...
        isfield(saved.operatorMetadata,'discretization') && ...
        isfield(saved.operatorMetadata.discretization, ...
        'denseGaugeRadialIndex')
    gauge = saved.operatorMetadata.discretization. ...
        denseGaugeRadialIndex;
end
end

function value = numeric_field(source,name,defaultValue)
value = defaultValue;
if isfield(source,name) && isnumeric(source.(name)) && ...
        isscalar(source.(name))
    value = source.(name);
end
end

function tf = is_saved_baseline(configuration,savedGauge, ...
        quadraticStep,cubicStep)
tf = configuration.gaugeRadialIndex == savedGauge && ...
    configuration.quadraticStep == quadraticStep && ...
    configuration.cubicStep == cubicStep;
end

function value = saved_baseline_result(saved,settings)
entry = saved.weaklyNonlinear.cross{ ...
    settings.targetIndex,settings.sourceIndex};
target = saved.weaklyNonlinear.modes{settings.targetIndex};
value = empty_result();
value.g21 = entry.g;
value.g21Raw = entry.gUnprojected;
value.coefficientAccepted = entry.validCubicScaling;
value.success = value.coefficientAccepted && ...
    isfinite(real(value.g21Raw)) && isfinite(imag(value.g21Raw));
value.termNames = {'mean','mixed','direct cubic'};
value.termProjections = saved_term_projections(entry,target);
value.modeDirectResiduals = cellfun(@(mode) ...
    mode.directResidual,saved.weaklyNonlinear.modes);
value.modeAdjointResiduals = cellfun(@(mode) ...
    mode.leftResidual,saved.weaklyNonlinear.modes);
value.modeGateToleranceUsed = ...
    saved.input.options.coefficientModeResidualTolerance;
value.strictModeGatePassed = max([value.modeDirectResiduals; ...
    value.modeAdjointResiduals]) <= value.modeGateToleranceUsed;
if isfield(target,'storageDiagnostics')
    value.targetConditionIndicator = ...
        target.storageDiagnostics.eigenpairConditionIndicator;
elseif isfield(target,'block') && isfield(target.block,'Bslow')
    value.targetConditionIndicator = norm(target.left)* ...
        norm(target.block.Bslow*target.vector)/ ...
        max(abs(target.normalization),eps);
end
qBB = saved_source_mean_field(saved,settings.sourceIndex);
qAB = entry.qAB;
value.forcedRelativeResiduals = [ ...
    forced_relative_residual(qBB),forced_relative_residual(qAB)];
value.wallClockSeconds = 0;
end

function projections = saved_term_projections(entry,target)
names = {'termMean','termMixed','termDirectCubic'};
if all(cellfun(@(name) isfield(entry,name) && ...
        ~isempty(entry.(name)),names))
    projections = projected_terms(entry,target);
    return;
end
projections = complex(nan(1,numel(names)));
if ~isfield(entry,'storageProjectionDiagnostics')
    return;
end
stored = entry.storageProjectionDiagnostics;
for index = 1:numel(names)
    position = find(strcmp({stored.name},names{index}),1);
    if ~isempty(position)
        projections(index) = stored(position).projection;
    end
end
end

function solution = saved_source_mean_field(saved,sourceIndex)
solution = struct();
if isfield(saved.weaklyNonlinear,'self') && ...
        numel(saved.weaklyNonlinear.self) >= sourceIndex && ...
        isstruct(saved.weaklyNonlinear.self{sourceIndex})
    entry = saved.weaklyNonlinear.self{sourceIndex};
    if isfield(entry,'qAA')
        solution = entry.qAA;
    elseif isfield(entry,'qAbarA')
        solution = entry.qAbarA;
    end
end
end

function value = saved_cross_raw(saved,settings)
entry = saved.weaklyNonlinear.cross{ ...
    settings.targetIndex,settings.sourceIndex};
value = saved.weaklyNonlinear.g( ...
    settings.targetIndex,settings.sourceIndex);
if isstruct(entry) && isfield(entry,'gUnprojected') && ...
        isscalar(entry.gUnprojected)
    value = entry.gUnprojected;
end
end

function filename = resolve_cache_file(filename,repositoryRoot)
filename = char(filename);
if ~isfile(filename)
    candidate = fullfile(repositoryRoot,filename);
    if isfile(candidate)
        filename = candidate;
    end
end
if ~isfile(filename)
    error('vi_verify_g21_numerics:MissingCache', ...
        'Could not find coefficient cache: %s',filename);
end
filename = char(java.io.File(filename).getCanonicalPath());
end

function validate_saved_output(saved)
required = {'input','parameters','operatorMetadata','weaklyNonlinear'};
for index = 1:numel(required)
    if ~isfield(saved,required{index})
        error('vi_verify_g21_numerics:MissingField', ...
            'Saved output is missing %s.',required{index});
    end
end
if ~isfield(saved.weaklyNonlinear,'modes') || ...
        ~iscell(saved.weaklyNonlinear.modes) || ...
        ~isfield(saved.weaklyNonlinear,'g')
    error('vi_verify_g21_numerics:MissingCoefficients', ...
        'Saved output lacks the retained modes or cubic matrix.');
end
end
