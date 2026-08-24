function [selectedModes, selectedOperatingModes, report] = ...
        vi_select_slow_operating_modes(parameters, requestedModes, ...
        numberOfModes, temporalCutoff, rootOptions, userOptions)
%VI_SELECT_SLOW_OPERATING_MODES Pick retained branches with slow envelopes.
%
% The full two-mode WNL coefficient engine is expensive and is meaningful
% only for operating-point modes whose continuous growth/decay is slow
% compared with the fast forcing. This helper scans reduced cylinder
% branches, keeps the first requested branch by default, and replaces later
% fast branches by the slowest admissible alternatives. Already-computed
% retained modes are reused so their selected Floquet roots remain stable.

if nargin < 5 || isempty(rootOptions)
    rootOptions = struct();
end
if nargin < 6 || isempty(userOptions)
    userOptions = struct();
end
options = default_options(userOptions, requestedModes);
validateattributes(numberOfModes, {'numeric'}, ...
    {'scalar','integer','>=',1,'<=',numel(requestedModes)});
validateattributes(temporalCutoff, {'numeric'}, ...
    {'scalar','integer','nonnegative'});

scanRootOptions = rootOptions;
scanRootOptions.verbose = false;
scanRootOptions.suppressSolverOutput = true;

candidates = build_candidate_modes(requestedModes, options);
candidateReports = repmat(candidate_report_template(), ...
    numel(candidates), 1);
operatingModes = cell(numel(candidates), 1);

if options.verbose
    fprintf(['\nSlow operating-mode scan at a/g0 = %.12g ', ...
        '(limit |Re(lambda)|/omegaStar <= %.3g)\n'], ...
        parameters.aAnalysis, options.maximumSlowRateRatio);
end

for candidateIndex = 1:numel(candidates)
    mode = candidates(candidateIndex);
    candidateReports(candidateIndex).mode = mode;
    candidateReports(candidateIndex).label = mode.label;
    try
        precomputedMode = find_precomputed_operating_mode( ...
            options.precomputedOperatingModes, mode);
        usedPrecomputed = ~isempty(precomputedMode);
        if usedPrecomputed
            operatingMode = precomputedMode;
        else
            operatingMode = vi_operating_point_floquet_mode( ...
                parameters, mode, temporalCutoff, scanRootOptions, []);
        end
        if ~isfield(operatingMode, 'label') || ...
                isempty(operatingMode.label)
            operatingMode.label = mode.label;
        end
        ratio = abs(real(operatingMode.lambda)) / ...
            max(parameters.omegaStar, eps);
        operatingMode.slowRateRatio = ratio;
        operatingModes{candidateIndex} = operatingMode;
        candidateReports(candidateIndex).success = true;
        candidateReports(candidateIndex).precomputed = usedPrecomputed;
        candidateReports(candidateIndex).lambda = operatingMode.lambda;
        candidateReports(candidateIndex).s = operatingMode.s;
        candidateReports(candidateIndex).slowRateRatio = ratio;
        candidateReports(candidateIndex).residual = ...
            operating_mode_residual(operatingMode);
        candidateReports(candidateIndex).passesSlowGate = ...
            ratio <= options.maximumSlowRateRatio;
        if options.verbose
            sourceSuffix = '';
            if usedPrecomputed
                sourceSuffix = ' (precomputed)';
            end
            fprintf(['  %-24s lambda = %.6g%+.6gi, ratio = ', ...
                '%.4g%s\n'], ...
                operatingMode.label, real(operatingMode.lambda), ...
                imag(operatingMode.lambda), ratio, sourceSuffix);
        end
    catch scanError
        candidateReports(candidateIndex).errorIdentifier = ...
            scanError.identifier;
        candidateReports(candidateIndex).errorMessage = ...
            scanError.message;
        if options.verbose
            fprintf('  %-24s failed: %s\n', mode.label, ...
                scanError.message);
        end
    end
end

[selectedIndices, reason] = choose_candidates( ...
    candidates, candidateReports, requestedModes, numberOfModes, options);
selectedModes = requestedModes(1:numberOfModes);
selectedOperatingModes = cell(numberOfModes, 1);
success = numel(selectedIndices) == numberOfModes;
if success
    selectedModes = candidates(selectedIndices);
    for modeIndex = 1:numberOfModes
        selectedOperatingModes{modeIndex} = ...
            operatingModes{selectedIndices(modeIndex)};
    end
end

report = struct();
report.requestedModes = requestedModes(1:numberOfModes);
report.candidates = candidateReports;
report.success = success;
report.reason = reason;
report.selectedIndices = selectedIndices;
report.selectedModes = selectedModes;
report.maximumSlowRateRatio = options.maximumSlowRateRatio;
report.preserveFirstMode = options.preserveFirstMode;
report.changed = success && ...
    ~same_mode_sequence(requestedModes(1:numberOfModes), selectedModes);
report.autoSelectionAllowsModeRecovery = ...
    options.autoSelectionAllowsModeRecovery;

if options.verbose && success
    fprintf('  selected retained modes:\n');
    for modeIndex = 1:numberOfModes
        fprintf('    %d: %s\n', modeIndex, ...
            selectedOperatingModes{modeIndex}.label);
    end
elseif options.verbose
    fprintf('  no slow replacement selected: %s\n', reason);
end
end

function options = default_options(userOptions, requestedModes)
options = struct();
options.maximumSlowRateRatio = 0.10;
options.preserveFirstMode = true;
options.candidateAzimuthalOrders = [];
options.candidateRadialIndices = [];
options.candidateQuasifrequencies = [];
options.verbose = true;
options.autoSelectionAllowsModeRecovery = true;
options.precomputedOperatingModes = {};
names = fieldnames(userOptions);
for nameIndex = 1:numel(names)
    options.(names{nameIndex}) = userOptions.(names{nameIndex});
end
if isempty(options.candidateAzimuthalOrders)
    options.candidateAzimuthalOrders = unique([requestedModes.m]);
end
if isempty(options.candidateRadialIndices)
    options.candidateRadialIndices = unique([requestedModes.radialIndex]);
end
if isempty(options.candidateQuasifrequencies)
    options.candidateQuasifrequencies = unique([requestedModes.s]);
end
validateattributes(options.maximumSlowRateRatio, {'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(options.preserveFirstMode, {'logical','numeric'}, ...
    {'scalar'});
validateattributes(options.autoSelectionAllowsModeRecovery, ...
    {'logical','numeric'}, {'scalar'});
validateattributes(options.candidateAzimuthalOrders, {'numeric'}, ...
    {'vector','integer','nonnegative','finite','nonempty'});
validateattributes(options.candidateRadialIndices, {'numeric'}, ...
    {'vector','integer','positive','finite','nonempty'});
validateattributes(options.candidateQuasifrequencies, {'numeric'}, ...
    {'vector','real','finite','nonempty'});
if ~(isempty(options.precomputedOperatingModes) || ...
        iscell(options.precomputedOperatingModes) || ...
        isstruct(options.precomputedOperatingModes))
    error('vi_select_slow_operating_modes:BadPrecomputedOperatingModes', ...
        'precomputedOperatingModes must be empty, a cell array, or struct array.');
end
options.preserveFirstMode = logical(options.preserveFirstMode);
options.autoSelectionAllowsModeRecovery = ...
    logical(options.autoSelectionAllowsModeRecovery);
end

function candidates = build_candidate_modes(requestedModes, options)
candidates = requestedModes(:).';
for candidateIndex = 1:numel(candidates)
    if ~isfield(candidates(candidateIndex), 'label') || ...
            isempty(candidates(candidateIndex).label)
        candidates(candidateIndex).label = ...
            mode_label(candidates(candidateIndex));
    end
end
for m = options.candidateAzimuthalOrders(:).'
    for radialIndex = options.candidateRadialIndices(:).'
        for s = options.candidateQuasifrequencies(:).'
            mode = candidates(1);
            mode.m = m;
            mode.radialIndex = radialIndex;
            mode.s = s;
            mode.label = mode_label(mode);
            if isempty(find_equivalent_mode(candidates, mode))
                candidates(end+1) = mode; %#ok<AGROW>
            end
        end
    end
end
end

function label = mode_label(mode)
sWrapped = wnl_wrap_quasifrequency(mode.s);
if abs(sWrapped) < 1.0e-12
    suffix = 'harmonic';
elseif abs(sWrapped-0.5) < 1.0e-12
    suffix = 'subharmonic';
else
    suffix = sprintf('s%g', sWrapped);
    suffix = strrep(suffix, '.', 'p');
    suffix = strrep(suffix, '-', 'm');
end
label = sprintf('m%d_l%d_%s', mode.m, mode.radialIndex, suffix);
end

function report = candidate_report_template()
report = struct('mode',struct(),'label','','success',false, ...
    'precomputed',false,'lambda',NaN,'s',NaN,'slowRateRatio',Inf, ...
    'residual',Inf,'passesSlowGate',false,'errorIdentifier','', ...
    'errorMessage','');
end

function operatingMode = find_precomputed_operating_mode( ...
        precomputedModes, target)
operatingMode = [];
if isempty(precomputedModes)
    return;
end
if iscell(precomputedModes)
    for modeIndex = 1:numel(precomputedModes)
        candidate = precomputedModes{modeIndex};
        if operating_mode_matches(candidate, target)
            operatingMode = candidate;
            return;
        end
    end
elseif isstruct(precomputedModes)
    for modeIndex = 1:numel(precomputedModes)
        candidate = precomputedModes(modeIndex);
        if operating_mode_matches(candidate, target)
            operatingMode = candidate;
            return;
        end
    end
end
end

function tf = operating_mode_matches(operatingMode, target)
tf = false;
if ~isstruct(operatingMode) || isempty(fieldnames(operatingMode)) || ...
        ~isfield(operatingMode, 'lambda') || ...
        isempty(operatingMode.lambda)
    return;
end
candidateMode = struct();
if isfield(operatingMode, 'problem') && ...
        isstruct(operatingMode.problem)
    problem = operatingMode.problem;
    if isfield(problem, 'm') && isfield(problem, 'radialIndex')
        candidateMode.m = problem.m;
        candidateMode.radialIndex = problem.radialIndex;
        if isfield(operatingMode, 's') && ~isempty(operatingMode.s)
            candidateMode.s = operatingMode.s;
        elseif isfield(problem, 'targetS')
            candidateMode.s = problem.targetS;
        else
            candidateMode.s = target.s;
        end
    end
end
if isfield(candidateMode, 'm')
    tf = equivalent_mode(candidateMode, target);
elseif isfield(operatingMode, 'label') && ...
        isfield(target, 'label') && ...
        strcmp(operatingMode.label, target.label)
    tf = true;
end
end

function residual = operating_mode_residual(operatingMode)
residual = Inf;
if ~isstruct(operatingMode)
    return;
end
if isfield(operatingMode, 'diagnostics') && ...
        isstruct(operatingMode.diagnostics) && ...
        isfield(operatingMode.diagnostics, 'vectorResidual') && ...
        ~isempty(operatingMode.diagnostics.vectorResidual)
    residual = operatingMode.diagnostics.vectorResidual;
elseif isfield(operatingMode, 'rootSearch') && ...
        isstruct(operatingMode.rootSearch) && ...
        isfield(operatingMode.rootSearch, 'diagnostics') && ...
        isstruct(operatingMode.rootSearch.diagnostics) && ...
        isfield(operatingMode.rootSearch.diagnostics, ...
        'relativeSingularResidual') && ...
        ~isempty( ...
        operatingMode.rootSearch.diagnostics.relativeSingularResidual)
    residual = ...
        operatingMode.rootSearch.diagnostics.relativeSingularResidual;
end
if ~(isnumeric(residual) && isscalar(residual) && isfinite(residual))
    residual = Inf;
end
end

function [selectedIndices, reason] = choose_candidates( ...
        candidates, candidateReports, requestedModes, numberOfModes, options)
selectedIndices = zeros(0, 1);
reason = '';
if options.preserveFirstMode
    firstIndex = find_equivalent_mode(candidates, requestedModes(1));
    if isempty(firstIndex) || ~candidateReports(firstIndex).success
        reason = 'the first requested branch could not be evaluated';
        return;
    end
    if ~candidateReports(firstIndex).passesSlowGate
        reason = ['the first requested branch itself fails the slow-', ...
            'envelope gate'];
        return;
    end
    selectedIndices(end+1, 1) = firstIndex;
end

available = find([candidateReports.success] & ...
    [candidateReports.passesSlowGate]).';
for existingIndex = selectedIndices(:).'
    available = available(arrayfun(@(idx) ...
        ~equivalent_mode(candidates(idx), candidates(existingIndex)), ...
        available));
end
[~, order] = sort([candidateReports(available).slowRateRatio], ...
    'ascend');
available = available(order);
needed = numberOfModes-numel(selectedIndices);
if numel(available) < needed
    reason = sprintf(['only %d slow candidate(s) were found, but %d ', ...
        'additional retained mode(s) are needed'], ...
        numel(available), needed);
    return;
end
selectedIndices = [selectedIndices; available(1:needed)];
reason = 'selected slow operating-point modes';
end

function index = find_equivalent_mode(modes, target)
index = [];
for modeIndex = 1:numel(modes)
    if equivalent_mode(modes(modeIndex), target)
        index = modeIndex;
        return;
    end
end
end

function tf = equivalent_mode(a, b)
tf = a.m == b.m && a.radialIndex == b.radialIndex && ...
    abs(wnl_wrap_quasifrequency(a.s-b.s)) <= ...
    1.0e-12*max([1,abs(a.s),abs(b.s)]);
end

function tf = same_mode_sequence(a, b)
tf = numel(a) == numel(b);
if ~tf
    return;
end
for modeIndex = 1:numel(a)
    if ~equivalent_mode(a(modeIndex), b(modeIndex))
        tf = false;
        return;
    end
end
end
