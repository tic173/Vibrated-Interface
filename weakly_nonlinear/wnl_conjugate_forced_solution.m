function [solutionBar, reusable] = wnl_conjugate_forced_solution( ...
        model, solution, label, userOpts)
%WNL_CONJUGATE_FORCED_SOLUTION Reuse a nonresonant forced solve by symmetry.
%
% If A*q=f has no neutral-mode border, complex conjugation together with
% the Floquet index map gives the forced field in the conjugate block.  The
% conjugated candidate is checked against the actual conjugate operator and
% forcing before it is marked valid.  Bordered/resonant solutions are not
% reused because their projection coordinates require an explicit mapping.

if nargin < 3 || isempty(label)
    label = [solution.spec.label, '_bar'];
end
if nargin < 4
    userOpts = struct();
end
opts = wnl_options(userOpts);
solutionBar = [];
reusable = false;
required = {'spec','field','forcing','lambda','projections','valid'};
for fieldIndex = 1:numel(required)
    if ~isfield(solution,required{fieldIndex})
        error('wnl_conjugate_forced_solution:MissingField', ...
            'The source forced solution is missing field %s.', ...
            required{fieldIndex});
    end
end
if ~isempty(solution.lambda) || ~isempty(solution.projections)
    return;
end

fieldBar = wnl_conjugate_field(model,solution.field,label);
forcingField = wnl_make_field(solution.spec,solution.forcing);
forcingBarField = wnl_conjugate_field(model,forcingField, ...
    [label, '_forcing']);
forcingBarField.spec.label = fieldBar.spec.label;
forcingBar = wnl_field_vector(forcingBarField);
vectorBar = wnl_field_vector(fieldBar);
blockBar = model.block(fieldBar.spec);
equationResidual = norm(blockBar.A*vectorBar-forcingBar);
forcingNorm = norm(forcingBar);
relativeEquationResidual = equationResidual/max(forcingNorm,eps);

solutionBar = solution;
solutionBar.spec = fieldBar.spec;
solutionBar.field = fieldBar;
solutionBar.vector = vectorBar;
solutionBar.forcing = forcingBar;
solutionBar.projections = complex(zeros(0,1));
solutionBar.lambda = complex(zeros(0,1));
solutionBar.quadraticResonance = false;
solutionBar.unborderedResidual = equationResidual;
solutionBar.equationResidual = equationResidual;
solutionBar.constraintResidual = 0;
solutionBar.forcingNorm = forcingNorm;
solutionBar.relativeEquationResidual = relativeEquationResidual;
solutionBar.relativeConstraintResidual = 0;
solutionBar.valid = logical(solution.valid) && ...
    isfinite(relativeEquationResidual) && ...
    relativeEquationResidual <= opts.forcedSolveResidualTolerance;
sourceExploratoryUsable = logical(solution.valid);
if isfield(solution,'exploratoryUsable')
    sourceExploratoryUsable = logical(solution.exploratoryUsable);
end
exploratoryTolerance = opts.forcedSolveResidualTolerance;
if ~isempty(opts.forcedExploratoryResidualTolerance)
    exploratoryTolerance = opts.forcedExploratoryResidualTolerance;
end
solutionBar.exploratoryResidualTolerance = exploratoryTolerance;
solutionBar.exploratoryUsable = sourceExploratoryUsable && ...
    isfinite(relativeEquationResidual) && ...
    relativeEquationResidual <= exploratoryTolerance;
solutionBar.solveDiagnostics = struct('flag',0,'totalIterations',0, ...
    'restarts',0,'relativeResidual',relativeEquationResidual, ...
    'reusedConjugate',true);
solutionBar.gmresDiagnostics = struct('attempted',false, ...
    'accepted',false,'flag',NaN,'iterations',0, ...
    'relativeResidual',NaN);
solutionBar.residualReport = [];
solutionBar.solveSeconds = 0;
solutionBar.reusedConjugate = true;
solutionBar.reuseSourceLabel = solution.spec.label;
reusable = solutionBar.exploratoryUsable;
if opts.verbose
    fprintf(['WNL conjugate forced-field reuse %s <- %s: ', ...
        'forcing-relative residual %.3e, reusable/exact = %d/%d\n'], ...
        solutionBar.spec.label,solution.spec.label, ...
        relativeEquationResidual,reusable,solutionBar.valid);
end
end
