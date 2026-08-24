function result = wnl_cross_coefficient(model, modeA, modeB, ...
    modeBBar, neutralModes, userOpts, precomputed)
%WNL_CROSS_COEFFICIENT Coefficient of A*|B|^2 in the A equation.
%
% g_AB = <psi_A,
%   2*C(phi_A,q_BbarB)
% + 2*C(phi_B,q_AbarB)
% + 2*C(phibar_B,q_AB)
% + 6*D(phi_A,phi_B,phibar_B)>.

if nargin < 5 || isempty(neutralModes)
    neutralModes = wnl_unique_modes({modeA, modeB, modeBBar});
end
if nargin < 6
    userOpts = struct();
end
if nargin < 7 || isempty(precomputed)
    precomputed = struct();
end
opts = wnl_options(userOpts);
reused = struct('qBbarB',false,'qAbarB',false,'qAB',false);

specBbarB = wnl_combine_spec(model, ...
    {modeB.spec, modeBBar.spec}, [1, 1], ...
    [modeB.spec.label, '_BbarB']);
if has_precomputed(precomputed,'qBbarB')
    qBbarB = validated_precomputed(precomputed.qBbarB, ...
        specBbarB,'qBbarB');
    reused.qBbarB = true;
else
    forcingBbarB = 2.0 * wnl_apply_quadratic(model, ...
        modeB.field, modeBBar.field, specBbarB);
    qBbarB = wnl_solve_forced(model, specBbarB, forcingBbarB, ...
        neutralModes, opts);
    reused.qBbarB = false;
end
if opts.forcedFailFast && forced_field_stops_branch(qBbarB,opts)
    result = stopped_result(qBbarB,[],[],reused, ...
        'mean field qBbarB');
    if opts.verbose
        fprintf('WNL cross %s <- |%s|^2 stopped: %s\n', ...
            modeA.spec.label,modeB.spec.label,result.message);
    end
    return;
end

specAbarB = wnl_combine_spec(model, ...
    {modeA.spec, modeBBar.spec}, [1, 1], ...
    [modeA.spec.label, '_AbarB']);
if has_precomputed(precomputed,'qAbarB')
    qAbarB = validated_precomputed(precomputed.qAbarB, ...
        specAbarB,'qAbarB');
    reused.qAbarB = true;
else
    forcingAbarB = 2.0 * wnl_apply_quadratic(model, ...
        modeA.field, modeBBar.field, specAbarB);
    qAbarB = wnl_solve_forced(model, specAbarB, forcingAbarB, ...
        neutralModes, opts);
    reused.qAbarB = false;
end
if opts.forcedFailFast && forced_field_stops_branch(qAbarB,opts)
    result = stopped_result(qBbarB,qAbarB,[],reused, ...
        'difference field qAbarB');
    if opts.verbose
        fprintf('WNL cross %s <- |%s|^2 stopped: %s\n', ...
            modeA.spec.label,modeB.spec.label,result.message);
    end
    return;
end

specAB = wnl_combine_spec(model, ...
    {modeA.spec, modeB.spec}, [1, 1], ...
    [modeA.spec.label, '_AB']);
if has_precomputed(precomputed,'qAB')
    qAB = validated_precomputed(precomputed.qAB,specAB,'qAB');
    reused.qAB = true;
else
    forcingAB = 2.0 * wnl_apply_quadratic(model, ...
        modeA.field, modeB.field, specAB);
    qAB = wnl_solve_forced(model, specAB, forcingAB, ...
        neutralModes, opts);
    reused.qAB = false;
end
if opts.verbose && any(structfun(@(value) value,reused))
    fprintf(['WNL cross %s <- |%s|^2 reused forced fields: ', ...
        'mean=%d, difference=%d, sum=%d\n'], ...
        modeA.spec.label,modeB.spec.label,reused.qBbarB, ...
        reused.qAbarB,reused.qAB);
end

quadraticResonance = qBbarB.quadraticResonance || ...
    qAbarB.quadraticResonance || qAB.quadraticResonance;
forcedSolvesValid = qBbarB.valid && qAbarB.valid && qAB.valid;
forcedSolvesExploratoryUsable = forced_field_usable(qBbarB) && ...
    forced_field_usable(qAbarB) && forced_field_usable(qAB);
if ~forcedSolvesValid && opts.stopOnUnconvergedForcedSolve
    result = base_result();
    result.validCubicScaling = false;
    result.forcedSolvesValid = false;
    result.forcedSolvesExploratoryUsable = false;
    result.qBbarB = qBbarB;
    result.qAbarB = qAbarB;
    result.qAB = qAB;
    result.reusedForcedFields = reused;
    result.message = ['A required O(A_j A_k) forced field failed the ', ...
        'forcing-relative residual gate. The cross-cubic coefficient ', ...
        'is withheld.'];
    return;
end
if quadraticResonance && opts.stopOnQuadraticResonance
    result = base_result();
    result.validCubicScaling = false;
    result.quadraticResonance = true;
    result.qBbarB = qBbarB;
    result.qAbarB = qAbarB;
    result.qAB = qAB;
    result.reusedForcedFields = reused;
    result.message = ['A combination field is quadratically resonant. ', ...
        'The nonresonant cross-Landau system is not the correct scaling.'];
    return;
end

termMean = 2.0 * wnl_apply_quadratic(model, ...
    modeA.field, qBbarB.field, modeA.spec);
termDifference = 2.0 * wnl_apply_quadratic(model, ...
    modeB.field, qAbarB.field, modeA.spec);
termSum = 2.0 * wnl_apply_quadratic(model, ...
    modeBBar.field, qAB.field, modeA.spec);
termDirect = 6.0 * wnl_apply_cubic(model, ...
    modeA.field, modeB.field, modeBBar.field, modeA.spec);
forcing = termMean + termDifference + termSum + termDirect;
g = (modeA.left' * forcing(:)) / modeA.normalization;

result = base_result();
result.validCubicScaling = ~quadraticResonance && forcedSolvesValid;
result.quadraticResonance = quadraticResonance;
result.forcedSolvesValid = forcedSolvesValid;
result.forcedSolvesExploratoryUsable = ...
    forcedSolvesExploratoryUsable;
result.g = g;
result.qBbarB = qBbarB;
result.qAbarB = qAbarB;
result.qAB = qAB;
result.reusedForcedFields = reused;
result.termMean = termMean;
result.termDifference = termDifference;
result.termSum = termSum;
result.termDirectCubic = termDirect;
result.cubicForcing = forcing;
result.message = '';
end

function result = base_result()
result = struct();
result.validCubicScaling = true;
result.g = NaN;
result.qBbarB = [];
result.qAbarB = [];
result.qAB = [];
result.termMean = [];
result.termDifference = [];
result.termSum = [];
result.termDirectCubic = [];
result.cubicForcing = [];
result.message = '';
result.quadraticResonance = false;
result.forcedSolvesValid = true;
result.forcedSolvesExploratoryUsable = true;
result.reusedForcedFields = struct('qBbarB',false, ...
    'qAbarB',false,'qAB',false);
end

function tf = forced_field_stops_branch(solution,opts)
tf = (~solution.valid && opts.stopOnUnconvergedForcedSolve) || ...
    (solution.quadraticResonance && opts.stopOnQuadraticResonance);
end

function result = stopped_result(qBbarB,qAbarB,qAB,reused,fieldName)
result = base_result();
result.validCubicScaling = false;
result.qBbarB = qBbarB;
result.qAbarB = qAbarB;
result.qAB = qAB;
result.reusedForcedFields = reused;
available = {qBbarB,qAbarB,qAB};
forcedValid = true;
quadraticResonance = false;
for fieldIndex = 1:numel(available)
    if isempty(available{fieldIndex})
        continue;
    end
    forcedValid = forcedValid && available{fieldIndex}.valid;
    quadraticResonance = quadraticResonance || ...
        available{fieldIndex}.quadraticResonance;
end
result.forcedSolvesValid = forcedValid;
result.forcedSolvesExploratoryUsable = all(cellfun( ...
    @forced_field_usable,available(~cellfun(@isempty,available)))) && ...
    all(~cellfun(@isempty,available));
result.quadraticResonance = quadraticResonance;
if quadraticResonance
    result.message = sprintf([ ...
        'Fail-fast after %s: the block is quadratically resonant. ', ...
        'Later combination fields were not solved because they cannot ', ...
        'restore the nonresonant cubic scaling.'],fieldName);
else
    result.message = sprintf([ ...
        'Fail-fast after %s: its forcing-relative residual failed the ', ...
        'required gate. Later combination fields were not solved ', ...
        'because they cannot restore this cross coefficient.'],fieldName);
end
end

function tf = has_precomputed(precomputed,name)
tf = isstruct(precomputed) && isfield(precomputed,name) && ...
    ~isempty(precomputed.(name));
end

function solution = validated_precomputed(solution,specExpected,name)
if ~isstruct(solution) || ~isfield(solution,'field') || ...
        ~isstruct(solution.field) || ~isfield(solution.field,'spec')
    error('wnl_cross_coefficient:BadPrecomputedField', ...
        'precomputed.%s must be a WNL forced-solution structure.',name);
end
if ~wnl_equivalent_spec(solution.field.spec,specExpected)
    error('wnl_cross_coefficient:PrecomputedBlockMismatch', ...
        'precomputed.%s occupies the wrong Floquet block.',name);
end
end

function tf = forced_field_usable(solution)
tf = isstruct(solution) && isfield(solution,'valid') && ...
    logical(solution.valid);
if isstruct(solution) && isfield(solution,'exploratoryUsable')
    tf = logical(solution.exploratoryUsable);
end
end
