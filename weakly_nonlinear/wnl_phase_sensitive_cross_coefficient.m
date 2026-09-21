function result = wnl_phase_sensitive_cross_coefficient( ...
        model,targetMode,sourceMode,targetModeBar, ...
        qSourceSource,qTargetBarSource,userOpts)
%WNL_PHASE_SENSITIVE_CROSS_COEFFICIENT Coefficient of bar(A)*B^2.
%
% For two distinct complex modes, the phase-sensitive monomial is allowed
% when 2*m_B-m_A=m_A and 2*s_B-s_A=s_A modulo one. Its coefficient is
%   h_AB=<psi_A,2*C(phi_B,q_barAB)+2*C(bar(phi_A),q_BB)
%                    +3*D(bar(phi_A),phi_B,phi_B)>.

if nargin < 7
    userOpts = struct();
end
opts = wnl_options(userOpts);
result = base_result();
result.allowedBySymmetry = wnl_phase_sensitive_coupling_allowed( ...
    targetMode,sourceMode);
if ~result.allowedBySymmetry
    result.g = 0;
    result.gUnprojected = 0;
    return;
end
if isempty(qSourceSource) || isempty(qTargetBarSource)
    result.validCubicScaling = false;
    result.forcedSolvesValid = false;
    result.forcedSolvesExploratoryUsable = false;
    result.message = ['A required second-order field for the ', ...
        'phase-sensitive cubic coefficient is unavailable.'];
    return;
end

result.forcedSolvesValid = ...
    qSourceSource.valid && qTargetBarSource.valid;
result.forcedSolvesExploratoryUsable = ...
    forced_field_usable(qSourceSource) && ...
    forced_field_usable(qTargetBarSource);
result.quadraticResonance = qSourceSource.quadraticResonance || ...
    qTargetBarSource.quadraticResonance;
if (~result.forcedSolvesValid && opts.stopOnUnconvergedForcedSolve) || ...
        (result.quadraticResonance && opts.stopOnQuadraticResonance)
    result.validCubicScaling = false;
    result.message = ['The phase-sensitive cubic coefficient was ', ...
        'withheld because a required second-order field is invalid or ', ...
        'quadratically resonant.'];
    return;
end

termMixed = 2.0*wnl_apply_quadratic(model,sourceMode.field, ...
    qTargetBarSource.field,targetMode.spec);
termSecond = 2.0*wnl_apply_quadratic(model,targetModeBar.field, ...
    qSourceSource.field,targetMode.spec);
termDirect = 3.0*wnl_apply_cubic(model,targetModeBar.field, ...
    sourceMode.field,sourceMode.field,targetMode.spec);
forcing = termMixed+termSecond+termDirect;
rawG = (targetMode.left'*forcing(:))/targetMode.normalization;
[g,reality] = wnl_project_modal_coefficient(rawG,targetMode,opts);

result.g = g;
result.gUnprojected = rawG;
result.coefficientReality = reality;
result.validCubicScaling = ~result.quadraticResonance && ...
    result.forcedSolvesValid && reality.accepted;
result.qSourceSource = qSourceSource;
result.qTargetBarSource = qTargetBarSource;
result.termMixed = termMixed;
result.termSecondHarmonic = termSecond;
result.termDirectCubic = termDirect;
result.cubicForcing = forcing;
end

function result = base_result()
result = struct('allowedBySymmetry',false,'validCubicScaling',true, ...
    'forcedSolvesValid',true,'forcedSolvesExploratoryUsable',true, ...
    'quadraticResonance',false,'g',0,'gUnprojected',0, ...
    'coefficientReality',struct(),'qSourceSource',[], ...
    'qTargetBarSource',[],'termMixed',[], ...
    'termSecondHarmonic',[],'termDirectCubic',[], ...
    'cubicForcing',[],'message','');
end

function tf = forced_field_usable(solution)
tf = isstruct(solution) && isfield(solution,'valid') && ...
    logical(solution.valid);
if isstruct(solution) && isfield(solution,'exploratoryUsable')
    tf = logical(solution.exploratoryUsable);
end
end
