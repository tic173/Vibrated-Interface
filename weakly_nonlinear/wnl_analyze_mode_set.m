function result = wnl_analyze_mode_set(model, specs, userOpts)
%WNL_ANALYZE_MODE_SET Self/cross coefficients for retained Floquet modes.
%
% Each spec may contain fields spec.direct and spec.left. For modes sharing
% a Floquet block, supplying these vectors is recommended so the intended
% null vector is selected.

if nargin < 3
    userOpts = struct();
end
opts = wnl_options(userOpts);
analysisWallClock = tic;
if ~iscell(specs)
    specs = num2cell(specs);
end
nModes = numel(specs);

if nModes > 1 && (~isempty(opts.direct) || ~isempty(opts.left))
    error('wnl_analyze_mode_set:GlobalModeVector', ...
        ['For multiple modes, put direct and left vectors in each spec ', ...
         'instead of opts.direct or opts.left.']);
end

modes = cell(nModes, 1);
bars = cell(nModes, 1);
referenceConsistency = cell(nModes,1);
modeRecoveryWallClock = tic;
for j = 1:nModes
    modes{j} = wnl_compute_mode(model, specs{j}, opts);
    bars{j} = wnl_conjugate_mode(model, modes{j}, opts);
    wnl_assert_mode_converged(modes{j}, opts);
    wnl_assert_mode_converged(bars{j}, opts);
    referenceConsistency{j} = ...
        wnl_assert_reference_eigenvalue_consistent(modes{j},opts);
    coefficientOpts = opts;
    coefficientOpts.modeResidualTolerance = ...
        opts.coefficientModeResidualTolerance;
    wnl_assert_mode_converged(modes{j}, coefficientOpts);
    % Conjugate direct modes enter every physical nonlinear field and use
    % the strict coefficient gate. Conjugate adjoints are only used to
    % detect/border a matching quadratic block and have already passed the
    % ordinary linear gate. A match is a quadratic-resonance diagnostic,
    % not part of the accepted nonresonant cubic reduction.
    coefficientBarOpts = coefficientOpts;
    coefficientBarOpts.checkAdjointModeResidual = false;
    wnl_assert_mode_converged(bars{j}, coefficientBarOpts);
end
modeRecoverySeconds = toc(modeRecoveryWallClock);

neutralModes = cell(2 * nModes, 1);
for j = 1:nModes
    neutralModes{2*j - 1} = modes{j};
    neutralModes{2*j} = bars{j};
end
neutralModes = wnl_unique_modes(neutralModes);

self = cell(nModes, 1);
cross = cell(nModes, nModes);
g = complex(nan(nModes, nModes));
valid = true(nModes, nModes);
coefficientWallClock = tic;
for a = 1:nModes
    self{a} = wnl_self_coefficient(model, modes{a}, bars{a}, ...
        neutralModes, opts);
    g(a, a) = self{a}.g;
    valid(a, a) = self{a}.validCubicScaling;
end

reusedForcedSolveCount = 0;
conjugateDifferenceReused = false;
twoModeSharingEnabled = nModes == 2 && ...
    opts.reuseTwoModeForcedFields;
if twoModeSharingEnabled
    % The two cross equations share four second-order fields:
    %   q_BbarB and q_AbarA are the already-computed self means;
    %   q_AB=q_BA is one common sum field;
    %   q_BbarA is the physical conjugate of q_AbarB when the difference
    %   block is nonresonant.  The conjugated candidate is checked with the
    %   full operator before reuse.  Thus the usual ten forced solves are
    %   reduced to six (or seven if the difference block is bordered).
    precomputed12 = struct('qBbarB',self{2}.qAbarA);
    cross{1,2} = wnl_cross_coefficient(model,modes{1},modes{2}, ...
        bars{2},neutralModes,opts,precomputed12);
    g(1,2) = cross{1,2}.g;
    valid(1,2) = cross{1,2}.validCubicScaling;

    precomputed21 = struct('qBbarB',self{1}.qAbarA);
    if ~isempty(cross{1,2}.qAB)
        precomputed21.qAB = cross{1,2}.qAB;
    end
    if ~isempty(cross{1,2}.qAbarB)
        [reverseDifference,conjugateDifferenceReused] = ...
            wnl_conjugate_forced_solution(model,cross{1,2}.qAbarB, ...
            [modes{2}.spec.label,'_AbarB'],opts);
        if conjugateDifferenceReused
            precomputed21.qAbarB = reverseDifference;
        end
    end
    cross{2,1} = wnl_cross_coefficient(model,modes{2},modes{1}, ...
        bars{1},neutralModes,opts,precomputed21);
    g(2,1) = cross{2,1}.g;
    valid(2,1) = cross{2,1}.validCubicScaling;
    reusedForcedSolveCount = count_reused(cross{1,2}) + ...
        count_reused(cross{2,1});
elseif nModes > 1
    % General fallback for callers retaining more than two modes.
    % The editable cylinder driver currently permits one or two modes.
    for a = 1:nModes
        for b = 1:nModes
            if a == b
                continue;
            end
            cross{a, b} = wnl_cross_coefficient(model, ...
                modes{a}, modes{b}, bars{b}, neutralModes, opts);
            g(a, b) = cross{a, b}.g;
            valid(a, b) = cross{a, b}.validCubicScaling;
        end
    end
end
baselineForcedSolveCount = 2*nModes+3*nModes*(nModes-1);
actualForcedSolveCount = 2*nModes;
for rowIndex = 1:nModes
    for columnIndex = 1:nModes
        if rowIndex ~= columnIndex && ~isempty(cross{rowIndex,columnIndex})
            actualForcedSolveCount = actualForcedSolveCount+ ...
                count_executed(cross{rowIndex,columnIndex});
        end
    end
end
skippedForcedSolveCount = max(baselineForcedSolveCount- ...
    actualForcedSolveCount-reusedForcedSolveCount,0);
if opts.verbose && twoModeSharingEnabled
    fprintf(['WNL two-mode forced fields: %d solve(s) executed, %d ', ...
        'shared field(s) reused, %d downstream solve(s) skipped by ', ...
        'fail-fast (baseline %d).\n'],actualForcedSolveCount, ...
        reusedForcedSolveCount,skippedForcedSolveCount, ...
        baselineForcedSolveCount);
end
coefficientSeconds = toc(coefficientWallClock);

mu = complex(nan(nModes, 1));
if isfield(model, 'detuning') && isa(model.detuning, 'function_handle') ...
        && ~isempty(opts.detuning)
    for j = 1:nModes
        direction = opts.detuning;
        if iscell(opts.detuning)
            direction = opts.detuning{j};
        end
        forcing = model.detuning(modes{j}.field, ...
            direction, modes{j}.spec);
        mu(j) = (modes{j}.left' * forcing(:)) / ...
            modes{j}.normalization;
    end
end

quadraticResonances = wnl_find_quadratic_resonances(model, ...
    modes, neutralModes, opts.resonanceTolerance);

result = struct();
result.modes = modes;
result.conjugateModes = bars;
result.neutralModes = neutralModes;
result.linearCoefficients = cellfun(@(mode) ...
    wnl_spec_lambda(mode.spec), modes);
result.referenceType = reference_type(result.linearCoefficients);
result.mu = mu;
result.g = g;
result.validCubicScaling = valid;
result.forcedSolvesValid = cellfun(@(entry) ...
    isempty(entry) || entry.forcedSolvesValid, cross);
result.forcedSolvesExploratoryUsable = cellfun(@(entry) ...
    isempty(entry) || entry.forcedSolvesExploratoryUsable, cross);
for j = 1:nModes
    result.forcedSolvesValid(j,j) = self{j}.forcedSolvesValid;
    result.forcedSolvesExploratoryUsable(j,j) = ...
        self{j}.forcedSolvesExploratoryUsable;
end
result.self = self;
result.cross = cross;
result.quadraticResonances = quadraticResonances;
result.referenceEigenvalueConsistency = referenceConsistency;
result.optimization = struct( ...
    'twoModeSharedForcedFields',twoModeSharingEnabled, ...
    'reusedForcedSolveCount',reusedForcedSolveCount, ...
    'conjugateDifferenceReused',conjugateDifferenceReused, ...
    'skippedForcedSolveCount',skippedForcedSolveCount, ...
    'baselineForcedSolveCount',baselineForcedSolveCount, ...
    'actualForcedSolveCount',actualForcedSolveCount);
result.timing = struct('modeRecoverySeconds',modeRecoverySeconds, ...
    'coefficientSeconds',coefficientSeconds, ...
    'totalSeconds',toc(analysisWallClock));
if opts.verbose
    fprintf(['WNL mode-set timing: recovery %.3f s, coefficients ', ...
        '%.3f s, total %.3f s\n'],modeRecoverySeconds, ...
        coefficientSeconds,result.timing.totalSeconds);
end
if strcmp(result.referenceType, 'operating-point Floquet reduction')
    result.amplitudeEquation = ...
        'dA_j/dt = lambda_j*A_j + sum_k g(j,k)*A_j*abs(A_k)^2';
else
    result.amplitudeEquation = ...
        'dA_j/dT = mu_j*A_j + sum_k g(j,k)*A_j*abs(A_k)^2';
end
end

function count = count_reused(crossResult)
count = 0;
if ~isstruct(crossResult) || ...
        ~isfield(crossResult,'reusedForcedFields')
    return;
end
values = struct2cell(crossResult.reusedForcedFields);
count = sum(cellfun(@(value) logical(value),values));
end

function count = count_executed(crossResult)
count = 0;
if ~isstruct(crossResult)
    return;
end
fieldNames = {'qBbarB','qAbarB','qAB'};
for fieldIndex = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIndex};
    if isfield(crossResult,fieldName) && ...
            ~isempty(crossResult.(fieldName)) && ...
            isfield(crossResult,'reusedForcedFields') && ...
            ~crossResult.reusedForcedFields.(fieldName)
        count = count+1;
    end
end
end

function value = reference_type(linearCoefficients)
if any(abs(linearCoefficients) > 1.0e-13)
    value = 'operating-point Floquet reduction';
else
    value = 'common-neutral-point reduction';
end
end
