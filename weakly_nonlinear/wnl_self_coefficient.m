function result = wnl_self_coefficient(model, mode, modeBar, neutralModes, userOpts)
%WNL_SELF_COEFFICIENT Compute mean, second-harmonic, and cubic feedback.
%
% The returned coefficient is
%   g = <psi, 2*C(phi,q_AbarA) + 2*C(phibar,q_AA)
%             + 3*D(phi,phi,phibar)>.

if nargin < 4 || isempty(neutralModes)
    neutralModes = wnl_unique_modes({mode, modeBar});
end
if nargin < 5
    userOpts = struct();
end
opts = wnl_options(userOpts);

specAA = wnl_combine_spec(model, {mode.spec, mode.spec}, ...
    [1, 1], [mode.spec.label, '_AA']);
forcingAA = wnl_apply_quadratic(model, ...
    mode.field, mode.field, specAA);
qAA = wnl_solve_forced(model, specAA, forcingAA, ...
    neutralModes, opts);

specAbarA = wnl_combine_spec(model, ...
    {mode.spec, modeBar.spec}, [1, 1], ...
    [mode.spec.label, '_AbarA']);
forcingAbarA = 2.0 * wnl_apply_quadratic(model, ...
    mode.field, modeBar.field, specAbarA);
qAbarA = wnl_solve_forced(model, specAbarA, forcingAbarA, ...
    neutralModes, opts);

quadraticResonance = qAA.quadraticResonance || ...
    qAbarA.quadraticResonance;
forcedSolvesValid = qAA.valid && qAbarA.valid;
forcedSolvesExploratoryUsable = forced_field_usable(qAA) && ...
    forced_field_usable(qAbarA);
if ~forcedSolvesValid && opts.stopOnUnconvergedForcedSolve
    result = base_result();
    result.validCubicScaling = false;
    result.forcedSolvesValid = false;
    result.forcedSolvesExploratoryUsable = false;
    result.qAA = qAA;
    result.qAbarA = qAbarA;
    result.message = ['A required O(A^2) forced field failed the ', ...
        'forcing-relative residual gate. The cubic coefficient is ', ...
        'withheld because projecting an inaccurate slaved field can ', ...
        'produce a spurious Landau coefficient.'];
    return;
end
if quadraticResonance && opts.stopOnQuadraticResonance
    result = base_result();
    result.validCubicScaling = false;
    result.quadraticResonance = true;
    result.g = NaN;
    result.qAA = qAA;
    result.qAbarA = qAbarA;
    result.message = ['A second-order forcing projects onto a retained ', ...
        'retained mode. Use the corresponding quadratic amplitude ', ...
        'scaling and derive the quadratic ', ...
        'amplitude system before cubic saturation.'];
    return;
end

termMean = 2.0 * wnl_apply_quadratic(model, ...
    mode.field, qAbarA.field, mode.spec);
termSecond = 2.0 * wnl_apply_quadratic(model, ...
    modeBar.field, qAA.field, mode.spec);
termDirect = 3.0 * wnl_apply_cubic(model, ...
    mode.field, mode.field, modeBar.field, mode.spec);
cubicForcing = termMean + termSecond + termDirect;
g = (mode.left' * cubicForcing(:)) / mode.normalization;

result = base_result();
result.validCubicScaling = ~quadraticResonance && forcedSolvesValid;
result.quadraticResonance = quadraticResonance;
result.forcedSolvesValid = forcedSolvesValid;
result.forcedSolvesExploratoryUsable = ...
    forcedSolvesExploratoryUsable;
result.g = g;
result.qAA = qAA;
result.qAbarA = qAbarA;
result.termMean = termMean;
result.termSecondHarmonic = termSecond;
result.termDirectCubic = termDirect;
result.cubicForcing = cubicForcing;
result.message = '';
end

function result = base_result()
result = struct();
result.validCubicScaling = true;
result.g = NaN;
result.qAA = [];
result.qAbarA = [];
result.termMean = [];
result.termSecondHarmonic = [];
result.termDirectCubic = [];
result.cubicForcing = [];
result.message = '';
result.quadraticResonance = false;
result.forcedSolvesValid = true;
result.forcedSolvesExploratoryUsable = true;
end

function tf = forced_field_usable(solution)
tf = isstruct(solution) && isfield(solution,'valid') && ...
    logical(solution.valid);
if isstruct(solution) && isfield(solution,'exploratoryUsable')
    tf = logical(solution.exploratoryUsable);
end
end
