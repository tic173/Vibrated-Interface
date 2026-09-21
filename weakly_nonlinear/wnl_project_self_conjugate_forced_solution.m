function [projected,diagnostics] = ...
        wnl_project_self_conjugate_forced_solution( ...
        model,solution,userOpts)
%WNL_PROJECT_SELF_CONJUGATE_FORCED_SOLUTION Project and recheck a real field.
%
% Primitive-variable DAE solves can carry an anti-real pressure-nullspace
% component even when the forcing and physical solution are real. A large
% raw conjugacy defect is therefore not, by itself, grounds for rejecting a
% coefficient. When needed, recheck the projected field in the original
% unequilibrated Floquet equations and accept only if that residual passes.

if nargin < 3
    userOpts = struct();
end
opts = wnl_options(userOpts);
required = {'field','forcing','valid','quadraticResonance'};
for fieldIndex = 1:numel(required)
    if ~isfield(solution,required{fieldIndex})
        error('wnl_project_self_conjugate_forced_solution:BadSolution', ...
            'solution.%s is required.',required{fieldIndex});
    end
end

[projected,diagnostics] = wnl_project_self_conjugate_field( ...
    model,solution.field,opts);
diagnostics.rawConjugacyAccepted = diagnostics.accepted;
diagnostics.projectedEquationResidual = NaN;
diagnostics.projectedRelativeEquationResidual = NaN;
diagnostics.projectedEquationTolerance = ...
    opts.forcedSolveResidualTolerance;
diagnostics.acceptedByProjectedEquation = false;
diagnostics.acceptanceReason = 'raw conjugacy defect passed';
if diagnostics.accepted
    return;
end

diagnostics.acceptanceReason = 'raw conjugacy defect failed';
if ~logical(solution.valid) || logical(solution.quadraticResonance)
    return;
end
if isfield(solution,'lambda') && ~isempty(solution.lambda)
    diagnostics.acceptanceReason = [ ...
        'projected residual recheck is unavailable for a bordered field'];
    return;
end
forcing = solution.forcing(:);
if any(~isfinite(real(forcing))) || any(~isfinite(imag(forcing)))
    diagnostics.acceptanceReason = 'the stored forcing is nonfinite';
    return;
end

block = model.block(solution.field.spec);
projectedVector = wnl_field_vector(projected);
equationResidual = norm(block.A*projectedVector-forcing);
relativeResidual = equationResidual/max(norm(forcing),eps);
diagnostics.projectedEquationResidual = equationResidual;
diagnostics.projectedRelativeEquationResidual = relativeResidual;
diagnostics.acceptedByProjectedEquation = ...
    isfinite(relativeResidual) && ...
    relativeResidual <= opts.forcedSolveResidualTolerance && ...
    diagnostics.projectedRealityResidual <= ...
    opts.realFieldConjugacyTolerance;
diagnostics.accepted = diagnostics.acceptedByProjectedEquation;
if diagnostics.accepted
    diagnostics.acceptanceReason = [ ...
        'projected field passed the original Floquet equation'];
else
    diagnostics.acceptanceReason = [ ...
        'projected field missed the original Floquet equation gate'];
end
end
