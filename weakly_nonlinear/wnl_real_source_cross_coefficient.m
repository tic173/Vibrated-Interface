function result = wnl_real_source_cross_coefficient(model,modeA,modeB, ...
        neutralModes,userOpts,precomputed)
%WNL_REAL_SOURCE_CROSS_COEFFICIENT Coefficient of A*B^2 for real B.
%
% The target A may be real or complex, while B is a self-conjugate real
% coordinate. The required fields and projection are
%   A(2*lambda_B) q_BB = C(phi_B,phi_B),
%   A(lambda_A+lambda_B) q_AB = 2*C(phi_A,phi_B),
%   g_AB = <psi_A, 2*C(phi_A,q_BB) + 2*C(phi_B,q_AB)
%                    + 3*D(phi_A,phi_B,phi_B)>.

if nargin < 4 || isempty(neutralModes)
    neutralModes = wnl_unique_modes({modeA,modeB});
end
if nargin < 5
    userOpts = struct();
end
if nargin < 6 || isempty(precomputed)
    precomputed = struct();
end
opts = wnl_options(userOpts);
reused = struct('qBB',false,'qAB',false);
targetIsReal = strcmp(wnl_mode_coordinate_type( ...
    modeA,opts.realModeEigenvalueTolerance),'real');

specBB = wnl_combine_spec(model,{modeB.spec,modeB.spec},[1,1], ...
    [modeB.spec.label,'_BB']);
if has_precomputed(precomputed,'qBB')
    qBB = validated_precomputed(precomputed.qBB,specBB,'qBB');
    reused.qBB = true;
else
    forcingBB = wnl_apply_quadratic( ...
        model,modeB.field,modeB.field,specBB);
    qBB = wnl_solve_forced(model,specBB,forcingBB,neutralModes,opts);
end
if opts.forcedFailFast && forced_field_stops_branch(qBB,opts)
    result = stopped_result(qBB,[],reused,'real source field qBB');
    return;
end

specAB = wnl_combine_spec(model,{modeA.spec,modeB.spec},[1,1], ...
    [modeA.spec.label,'_AB']);
if has_precomputed(precomputed,'qAB')
    qAB = validated_precomputed(precomputed.qAB,specAB,'qAB');
    reused.qAB = true;
else
    forcingAB = 2.0*wnl_apply_quadratic( ...
        model,modeA.field,modeB.field,specAB);
    qAB = wnl_solve_forced(model,specAB,forcingAB,neutralModes,opts);
end

quadraticResonance = qBB.quadraticResonance || qAB.quadraticResonance;
forcedSolvesValid = qBB.valid && qAB.valid;
forcedSolvesExploratoryUsable = ...
    forced_field_usable(qBB) && forced_field_usable(qAB);
if ~forcedSolvesValid && opts.stopOnUnconvergedForcedSolve
    result = stopped_result(qBB,qAB,reused,'real mixed field qAB');
    return;
end
if quadraticResonance && opts.stopOnQuadraticResonance
    result = stopped_result(qBB,qAB,reused,'quadratically resonant field');
    result.quadraticResonance = true;
    result.message = ['A real-source combination field is quadratically ', ...
        'resonant. The cubic-only amplitude system is invalid.'];
    return;
end

[qBBField,qBBReality] = ...
    wnl_project_self_conjugate_forced_solution(model,qBB,opts);
[qABField,qABReality] = project_forced_if_real( ...
    model,qAB,targetIsReal,opts);
termMean = 2.0*wnl_apply_quadratic( ...
    model,modeA.field,qBBField,modeA.spec);
termMixed = 2.0*wnl_apply_quadratic( ...
    model,modeB.field,qABField,modeA.spec);
termDirect = 3.0*wnl_apply_cubic( ...
    model,modeA.field,modeB.field,modeB.field,modeA.spec);
[termMeanField,termMeanReality] = project_if_real(model, ...
    wnl_make_field(modeA.spec,termMean),targetIsReal,opts);
[termMixedField,termMixedReality] = project_if_real(model, ...
    wnl_make_field(modeA.spec,termMixed),targetIsReal,opts);
[termDirectField,termDirectReality] = project_if_real(model, ...
    wnl_make_field(modeA.spec,termDirect),targetIsReal,opts);
termMean = termMeanField.coeff;
termMixed = termMixedField.coeff;
termDirect = termDirectField.coeff;
forcing = termMean+termMixed+termDirect;
[forcingField,forcingReality] = project_if_real(model, ...
    wnl_make_field(modeA.spec,forcing),targetIsReal,opts);
forcing = forcingField.coeff;
rawG = (modeA.left'*forcing(:))/modeA.normalization;
[g,reality] = wnl_project_modal_coefficient(rawG,modeA,opts);

result = base_result();
result.validCubicScaling = ~quadraticResonance && ...
    forcedSolvesValid && qBBReality.accepted && qABReality.accepted && ...
    termMeanReality.accepted && termMixedReality.accepted && ...
    termDirectReality.accepted && ...
    forcingReality.accepted && reality.accepted;
result.quadraticResonance = quadraticResonance;
result.forcedSolvesValid = forcedSolvesValid;
result.forcedSolvesExploratoryUsable = forcedSolvesExploratoryUsable;
result.g = g;
result.gUnprojected = rawG;
result.coefficientReality = reality;
result.qBB = qBB;
result.qBbarB = qBB;
result.qAB = qAB;
result.reusedForcedFields = reused;
result.termMean = termMean;
result.termMixed = termMixed;
result.termDirectCubic = termDirect;
result.cubicForcing = forcing;
result.fieldReality = struct('qBB',qBBReality,'qAB',qABReality, ...
    'termMean',termMeanReality,'termMixed',termMixedReality, ...
    'termDirectCubic',termDirectReality, ...
    'cubicForcing',forcingReality);
if ~reality.accepted
    result.message = sprintf([ ...
        'A modal coefficient had imaginary part %.3e, above the ', ...
        'allowed %.3e. The coefficient was withheld.'], ...
        abs(imag(rawG)),reality.allowedImaginaryPart);
elseif ~qBBReality.accepted || ~qABReality.accepted || ...
        ~termMeanReality.accepted || ~termMixedReality.accepted || ...
        ~termDirectReality.accepted || ...
        ~forcingReality.accepted
    result.message = sprintf([ ...
        'A real forced field had conjugacy defect %.3e, above the ', ...
        'allowed %.3e. The coefficient was withheld.'], ...
        max([qBBReality.relativeDefect,qABReality.relativeDefect, ...
        termMeanReality.relativeDefect,termMixedReality.relativeDefect, ...
        termDirectReality.relativeDefect,forcingReality.relativeDefect]), ...
        opts.realFieldConjugacyTolerance);
end
end

function [field,diagnostics] = project_if_real( ...
        model,field,targetIsReal,opts)
if targetIsReal
    [field,diagnostics] = ...
        wnl_project_self_conjugate_field(model,field,opts);
else
    diagnostics = struct('applied',false,'accepted',true, ...
        'relativeDefect',0);
end
end

function [field,diagnostics] = project_forced_if_real( ...
        model,solution,targetIsReal,opts)
if targetIsReal
    [field,diagnostics] = ...
        wnl_project_self_conjugate_forced_solution( ...
        model,solution,opts);
else
    field = solution.field;
    diagnostics = struct('applied',false,'accepted',true, ...
        'relativeDefect',0);
end
end

function result = base_result()
result = struct('validCubicScaling',true,'g',NaN, ...
    'gUnprojected',NaN,'qBB',[],'qBbarB',[], ...
    'qAbarB',[],'qAB',[],'termMean',[],'termDifference',[], ...
    'termSum',[],'termMixed',[],'termDirectCubic',[], ...
    'cubicForcing',[],'message','','quadraticResonance',false, ...
    'forcedSolvesValid',true,'forcedSolvesExploratoryUsable',true, ...
    'reusedForcedFields',struct('qBB',false,'qAB',false), ...
    'coordinateConvention','real-source', ...
    'coefficientReality',struct(),'fieldReality',struct());
end

function tf = forced_field_stops_branch(solution,opts)
tf = (~solution.valid && opts.stopOnUnconvergedForcedSolve) || ...
    (solution.quadraticResonance && opts.stopOnQuadraticResonance);
end

function result = stopped_result(qBB,qAB,reused,fieldName)
result = base_result();
result.validCubicScaling = false;
result.qBB = qBB;
result.qBbarB = qBB;
result.qAB = qAB;
result.reusedForcedFields = reused;
available = {qBB,qAB};
available = available(~cellfun(@isempty,available));
result.forcedSolvesValid = all(cellfun(@(q) q.valid,available));
result.forcedSolvesExploratoryUsable = ...
    all(cellfun(@forced_field_usable,available));
result.quadraticResonance = ...
    any(cellfun(@(q) q.quadraticResonance,available));
if result.quadraticResonance
    reason = 'the block is quadratically resonant';
else
    reason = 'its forcing-relative residual failed the required gate';
end
result.message = sprintf('Fail-fast after %s: %s.',fieldName,reason);
end

function tf = has_precomputed(precomputed,name)
tf = isstruct(precomputed) && isfield(precomputed,name) && ...
    ~isempty(precomputed.(name));
end

function solution = validated_precomputed(solution,specExpected,name)
if ~isstruct(solution) || ~isfield(solution,'field') || ...
        ~isstruct(solution.field) || ~isfield(solution.field,'spec')
    error('wnl_real_source_cross_coefficient:BadPrecomputedField', ...
        'precomputed.%s must be a WNL forced-solution structure.',name);
end
if ~wnl_equivalent_spec(solution.field.spec,specExpected)
    error('wnl_real_source_cross_coefficient:PrecomputedBlockMismatch', ...
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
