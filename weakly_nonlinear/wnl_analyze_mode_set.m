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
modeSetValidity = wnl_validate_cubic_mode_set(specs);
coefficientComputed = requested_coefficients( ...
    opts.coefficientPairs,nModes);

if nModes > 1 && (~isempty(opts.direct) || ~isempty(opts.left))
    error('wnl_analyze_mode_set:GlobalModeVector', ...
        ['For multiple modes, put direct and left vectors in each spec ', ...
         'instead of opts.direct or opts.left.']);
end

modes = cell(nModes, 1);
bars = cell(nModes, 1);
coordinateTypes = cell(nModes,1);
realCoordinateProjection = cell(nModes,1);
referenceConsistency = cell(nModes,1);
modeRecoveryWallClock = tic;
for j = 1:nModes
    modes{j} = wnl_compute_mode(model, specs{j}, opts);
    bars{j} = wnl_conjugate_mode(model, modes{j}, opts);
    coordinateTypes{j} = wnl_mode_coordinate_type( ...
        modes{j},opts.realModeEigenvalueTolerance, ...
        opts.realModeMinimumConjugateOverlap);
    if strcmp(coordinateTypes{j},'real')
        [modes{j},realCoordinateProjection{j}] = ...
            wnl_enforce_real_mode(model,modes{j},bars{j},opts);
        bars{j} = wnl_conjugate_mode(model,modes{j},opts);
        bars{j}.coordinateType = 'real';
    else
        modes{j}.coordinateType = 'complex';
        bars{j}.coordinateType = 'complex';
        realCoordinateProjection{j} = struct('applied',false);
    end
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
valid = false(nModes, nModes);
selfNeeded = diag(coefficientComputed).';
for targetIndex = 1:nModes
    for sourceIndex = 1:nModes
        if coefficientComputed(targetIndex,sourceIndex) && ...
                targetIndex ~= sourceIndex && ...
                strcmp(coordinateTypes{targetIndex},'complex') && ...
                strcmp(coordinateTypes{sourceIndex},'complex') && ...
                wnl_phase_sensitive_coupling_allowed( ...
                modes{targetIndex},modes{sourceIndex})
            selfNeeded(sourceIndex) = true;
        end
    end
end
coefficientWallClock = tic;
for a = 1:nModes
    if ~selfNeeded(a)
        continue;
    end
    self{a} = wnl_self_coefficient(model, modes{a}, bars{a}, ...
        neutralModes, opts);
    if coefficientComputed(a,a)
        g(a, a) = self{a}.g;
        valid(a, a) = self{a}.validCubicScaling;
    end
end

reusedForcedSolveCount = 0;
conjugateDifferenceReused = false;
twoModeSharingEnabled = nModes == 2 && ...
    opts.reuseTwoModeForcedFields && all(coefficientComputed(:));
if twoModeSharingEnabled
    % The two cross equations share four second-order fields:
    %   q_BbarB and q_AbarA are the already-computed self means;
    %   q_AB=q_BA is one common sum field;
    %   q_BbarA is the physical conjugate of q_AbarB when the difference
    %   block is nonresonant.  The conjugated candidate is checked with the
    %   full operator before reuse.  Thus the usual ten forced solves are
    %   reduced to six (or seven if the difference block is bordered).
    precomputed12 = source_self_precomputed( ...
        self{2},coordinateTypes{2});
    cross{1,2} = wnl_cross_coefficient(model,modes{1},modes{2}, ...
        bars{2},neutralModes,opts,precomputed12);
    g(1,2) = cross{1,2}.g;
    valid(1,2) = cross{1,2}.validCubicScaling;

    precomputed21 = source_self_precomputed( ...
        self{1},coordinateTypes{1});
    if ~isempty(cross{1,2}.qAB)
        precomputed21.qAB = cross{1,2}.qAB;
    end
    reverseDifferenceSource = [];
    if strcmp(coordinateTypes{1},'complex') && ...
            strcmp(coordinateTypes{2},'complex') && ...
            isfield(cross{1,2},'qAbarB') && ...
            ~isempty(cross{1,2}.qAbarB)
        reverseDifferenceSource = cross{1,2}.qAbarB;
    elseif strcmp(coordinateTypes{1},'complex') && ...
            strcmp(coordinateTypes{2},'real') && ...
            isfield(cross{1,2},'qAB') && ...
            ~isempty(cross{1,2}.qAB)
        reverseDifferenceSource = cross{1,2}.qAB;
    end
    if ~isempty(reverseDifferenceSource)
        [reverseDifference,conjugateDifferenceReused] = ...
            wnl_conjugate_forced_solution(model,reverseDifferenceSource, ...
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
            if a == b || ~coefficientComputed(a,b)
                continue;
            end
            cross{a, b} = wnl_cross_coefficient(model, ...
                modes{a}, modes{b}, bars{b}, neutralModes, opts);
            g(a, b) = cross{a, b}.g;
            valid(a, b) = cross{a, b}.validCubicScaling;
        end
    end
end

phaseSensitive = cell(nModes,nModes);
phaseSensitiveG = complex(zeros(nModes,nModes));
for targetIndex = 1:nModes
    for sourceIndex = 1:nModes
        if ~coefficientComputed(targetIndex,sourceIndex) || ...
                targetIndex == sourceIndex || ...
                ~strcmp(coordinateTypes{targetIndex},'complex') || ...
                ~strcmp(coordinateTypes{sourceIndex},'complex') || ...
                ~wnl_phase_sensitive_coupling_allowed( ...
                modes{targetIndex},modes{sourceIndex})
            continue;
        end
        reverseCross = cross{sourceIndex,targetIndex};
        mixedField = [];
        if isstruct(reverseCross) && ...
                isfield(reverseCross,'qAbarB')
            mixedField = reverseCross.qAbarB;
        end
        phaseSensitive{targetIndex,sourceIndex} = ...
            wnl_phase_sensitive_cross_coefficient( ...
            model,modes{targetIndex},modes{sourceIndex}, ...
            bars{targetIndex},self{sourceIndex}.qAA, ...
            mixedField,opts);
        phaseSensitiveG(targetIndex,sourceIndex) = ...
            phaseSensitive{targetIndex,sourceIndex}.g;
        valid(targetIndex,sourceIndex) = ...
            valid(targetIndex,sourceIndex) && ...
            phaseSensitive{targetIndex,sourceIndex}.validCubicScaling;
    end
end
baselineForcedSolveCount = 0;
actualForcedSolveCount = 0;
for modeIndex = 1:nModes
    if selfNeeded(modeIndex)
        baselineForcedSolveCount = baselineForcedSolveCount+ ...
            1+strcmp(coordinateTypes{modeIndex},'complex');
    end
    actualForcedSolveCount = actualForcedSolveCount+ ...
        count_self_executed(self{modeIndex});
end
for rowIndex = 1:nModes
    for columnIndex = 1:nModes
        if rowIndex ~= columnIndex && ...
                coefficientComputed(rowIndex,columnIndex)
            baselineForcedSolveCount = baselineForcedSolveCount+ ...
                2+strcmp(coordinateTypes{columnIndex},'complex');
        end
    end
end
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
if opts.verbose && (twoModeSharingEnabled || ...
        ~isempty(opts.coefficientPairs))
    if isempty(opts.coefficientPairs)
        prefix = 'WNL two-mode forced fields';
    else
        prefix = 'WNL targeted coefficient forced fields';
    end
    fprintf(['%s: %d solve(s) executed, %d shared field(s) reused, ', ...
        '%d downstream solve(s) skipped by fail-fast (baseline %d).\n'], ...
        prefix,actualForcedSolveCount,reusedForcedSolveCount, ...
        skippedForcedSolveCount,baselineForcedSolveCount);
end
coefficientSeconds = toc(coefficientWallClock);

mu = complex(nan(nModes, 1));
muReality = cell(nModes,1);
if isfield(model, 'detuning') && isa(model.detuning, 'function_handle') ...
        && ~isempty(opts.detuning)
    for j = 1:nModes
        direction = opts.detuning;
        if iscell(opts.detuning)
            direction = opts.detuning{j};
        end
        forcing = model.detuning(modes{j}.field, ...
            direction, modes{j}.spec);
        rawMu = (modes{j}.left' * forcing(:)) / ...
            modes{j}.normalization;
        [mu(j),muReality{j}] = wnl_project_modal_coefficient( ...
            rawMu,modes{j},opts);
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
result.gPhysicalPeak = wnl_physical_cubic_coefficients(g,modes);
result.phaseSensitiveG = phaseSensitiveG;
result.phaseSensitiveGPhysicalPeak = ...
    wnl_physical_cubic_coefficients(phaseSensitiveG,modes);
result.phaseSensitive = phaseSensitive;
result.hasPhaseSensitiveCoupling = any(phaseSensitiveG(:) ~= 0);
result.modeSetValidity = modeSetValidity;
result.coordinateTypes = coordinateTypes;
result.amplitudeScalesToPeakZetaOverH = ...
    wnl_mode_amplitude_scale(modes);
result.realCoordinateProjection = realCoordinateProjection;
result.detuningCoefficientReality = muReality;
result.coefficientConventionVersion = ...
    'V2-real-self-conjugate-physical-peak';
result.coefficientComputed = coefficientComputed;
result.validCubicScaling = valid;
result.forcedSolvesValid = false(nModes,nModes);
result.forcedSolvesExploratoryUsable = false(nModes,nModes);
for targetIndex = 1:nModes
    for sourceIndex = 1:nModes
        if targetIndex ~= sourceIndex && ...
                coefficientComputed(targetIndex,sourceIndex) && ...
                ~isempty(cross{targetIndex,sourceIndex})
            result.forcedSolvesValid(targetIndex,sourceIndex) = ...
                cross{targetIndex,sourceIndex}.forcedSolvesValid;
            result.forcedSolvesExploratoryUsable( ...
                targetIndex,sourceIndex) = ...
                cross{targetIndex,sourceIndex}. ...
                forcedSolvesExploratoryUsable;
        end
    end
end
for j = 1:nModes
    if coefficientComputed(j,j)
        result.forcedSolvesValid(j,j) = self{j}.forcedSolvesValid;
        result.forcedSolvesExploratoryUsable(j,j) = ...
            self{j}.forcedSolvesExploratoryUsable;
    end
end
result.self = self;
result.cross = cross;
result.quadraticResonances = quadraticResonances;
result.referenceEigenvalueConsistency = referenceConsistency;
result.optimization = struct( ...
    'targetedCoefficientEvaluation',~isempty(opts.coefficientPairs), ...
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
        ['dA_j/dt = lambda_j*A_j + sum_k ', ...
         '[g(j,k)*A_j*abs(A_k)^2 + h(j,k)*conj(A_j)*A_k^2]'];
else
    result.amplitudeEquation = ...
        ['dA_j/dT = mu_j*A_j + sum_k ', ...
         '[g(j,k)*A_j*abs(A_k)^2 + h(j,k)*conj(A_j)*A_k^2]'];
end
end

function requested = requested_coefficients(pairs,nModes)
if isempty(pairs)
    requested = true(nModes,nModes);
    return;
end
validateattributes(pairs,{'numeric'}, ...
    {'2d','integer','positive','finite','ncols',2});
if any(pairs(:) > nModes)
    error('wnl_analyze_mode_set:CoefficientPairIndex', ...
        ['opts.coefficientPairs contains a mode index larger than the ', ...
         'number of retained modes (%d).'],nModes);
end
requested = false(nModes,nModes);
for pairIndex = 1:size(pairs,1)
    requested(pairs(pairIndex,1),pairs(pairIndex,2)) = true;
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
if ~isfield(crossResult,'reusedForcedFields')
    return;
end
fieldNames = fieldnames(crossResult.reusedForcedFields);
for fieldIndex = 1:numel(fieldNames)
    fieldName = fieldNames{fieldIndex};
    if isfield(crossResult,fieldName) && ...
            ~isempty(crossResult.(fieldName)) && ...
            ~crossResult.reusedForcedFields.(fieldName)
        count = count+1;
    end
end
end

function count = count_self_executed(selfResult)
count = 0;
for name = {'qAA','qAbarA'}
    if isstruct(selfResult) && isfield(selfResult,name{1}) && ...
            ~isempty(selfResult.(name{1}))
        count = count+1;
    end
end
end

function precomputed = source_self_precomputed(selfResult,coordinateType)
if strcmp(coordinateType,'real')
    precomputed = struct('qBB',selfResult.qAA);
else
    precomputed = struct('qBbarB',selfResult.qAbarA);
end
end

function value = reference_type(linearCoefficients)
if any(abs(linearCoefficients) > 1.0e-13)
    value = 'operating-point Floquet reduction';
else
    value = 'common-neutral-point reduction';
end
end
