function discretization = vi_cylinder_wnl_discretization(parameters)
%VI_CYLINDER_WNL_DISCRETIZATION Fixed-domain cylinder grids and state map.
%
% Each azimuthal/temporal coefficient uses the primitive state
%   [u_r^d,u_theta^d,w^d,p^d,u_r^l,u_theta^l,w^l,p^l,zeta].
% Fluid arrays are Nr-by-Nz and zeta is Nr-by-1. The same layout is used
% for every azimuthal block, which is required by the convolution engine.
% Nz is either one global Chebyshev grid or the concatenation of several
% Chebyshev elements. Multi-domain grids retain both copies of every
% artificial interface; the linear operator supplies their matching rows.
% Radial values remain collocated, but Dr/Drr may be generated from a
% Bessel-enriched space containing the retained cylindrical mode branches.

required = {'R0', 'C', 'Bd', 'At', 'eta', 'omegaStar', ...
    'g_sgn', 'aCritical', 'phase', 'numerics'};
for fieldIndex = 1:numel(required)
    assert(isfield(parameters, required{fieldIndex}), ...
        'parameters.%s is required.', required{fieldIndex});
end

numerics = parameters.numerics;
requiredNumerics = {'Nr'};
for fieldIndex = 1:numel(requiredNumerics)
    assert(isfield(numerics, requiredNumerics{fieldIndex}), ...
        'parameters.numerics.%s is required.', ...
        requiredNumerics{fieldIndex});
end

nr = numerics.Nr;
validateattributes(nr, {'numeric'}, ...
    {'scalar', 'integer', '>=', 6});

if isfield(numerics, 'Ntheta')
    ntheta = numerics.Ntheta;
else
    ntheta = 32;
end
validateattributes(ntheta, {'numeric'}, ...
    {'scalar', 'integer', '>=', 12});

if isfield(numerics, 'quadraticStep')
    quadraticStep = numerics.quadraticStep;
else
    quadraticStep = 2.0e-4;
end
if isfield(numerics, 'cubicStep')
    cubicStep = numerics.cubicStep;
else
    cubicStep = 2.0e-3;
end
validateattributes(quadraticStep, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'});
validateattributes(cubicStep, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'});

radial = vi_cylinder_radial_grid(parameters);
r = radial.r;
verticalD = vi_cylinder_vertical_grid(parameters, 'lower');
verticalL = vi_cylinder_vertical_grid(parameters, 'upper');
nzD = verticalD.numberOfPoints;
nzL = verticalL.numberOfPoints;

nD = nr*nzD;
nL = nr*nzL;
next = 0;
layout.d.ur = next + (1:nD); next = next+nD;
layout.d.ut = next + (1:nD); next = next+nD;
layout.d.w = next + (1:nD); next = next+nD;
layout.d.p = next + (1:nD); next = next+nD;
layout.l.ur = next + (1:nL); next = next+nL;
layout.l.ut = next + (1:nL); next = next+nL;
layout.l.w = next + (1:nL); next = next+nL;
layout.l.p = next + (1:nL); next = next+nL;
layout.zeta = next + (1:nr); next = next+nr;
layout.ndof = next;
layout.nr = nr;
layout.nzD = nzD;
layout.nzL = nzL;
layout.nD = nD;
layout.nL = nL;

contactLine = 'free';
if isfield(parameters, 'boundary') && ...
        isfield(parameters.boundary, 'contactLine')
    contactLine = lower(char(parameters.boundary.contactLine));
end
if ~ismember(contactLine, {'free', 'pinned'})
    error('vi_cylinder_wnl_discretization:ContactLine', ...
        'boundary.contactLine must be ''free'' or ''pinned''.');
end

rhoD = 1;
rhoL = (1-parameters.At)/(1+parameters.At);
muD = 1;
muL = parameters.eta;
hAtwood = 2*parameters.At/(1+parameters.At);

discretization = struct();
discretization.parameters = parameters;
discretization.layout = layout;
discretization.r = r;
discretization.zD = verticalD.z;
discretization.zL = verticalL.z;
discretization.Dr = sparse(radial.D);
discretization.Drr = sparse(radial.D2);
discretization.radial = radial;
discretization.DzD = sparse(verticalD.D);
discretization.DzzD = sparse(verticalD.D2);
discretization.DzL = sparse(verticalL.D);
discretization.DzzL = sparse(verticalL.D2);
discretization.verticalD = verticalD;
discretization.verticalL = verticalL;
discretization.denseGaugeVerticalIndex = ...
    gauge_vertical_index(verticalD, -0.5);
discretization.ntheta = ntheta;
discretization.quadraticStep = quadraticStep;
discretization.cubicStep = cubicStep;
discretization.contactLine = contactLine;
discretization.rhoD = rhoD;
discretization.rhoL = rhoL;
discretization.muD = muD;
discretization.muL = muL;
discretization.hAtwood = hAtwood;
discretization.ndof = layout.ndof;
discretization.sidewallConditions = ...
    {'u_r=0', 'd_r(r*u_theta)=0', 'd_r(w)=0'};
discretization.useParallelNonlinearActions = false;
discretization.reportTiming = false;
discretization.nonlinearTemporalOversampling = 2.0;
% The implementation removes the exact flat-interface linear remainder
% before differencing. The remaining rational ALE and surface terms can
% still suffer small-step roundoff or large-step truncation. The adaptive
% stencil evaluates a short step-doubling sequence and uses Richardson
% agreement to select a cancellation-safe result. These execution fields do
% not alter the linear operator or the recovery-cache signature.
discretization.adaptiveDirectionalSteps = true;
discretization.quadraticStepMultipliers = [1,2,4,8];
discretization.cubicStepMultipliers = [1,2,4];
if isfield(parameters,'execution')
    if isfield(parameters.execution,'useParallelNonlinearActions')
        discretization.useParallelNonlinearActions = logical( ...
            parameters.execution.useParallelNonlinearActions);
    end
    if isfield(parameters.execution,'reportTiming')
        discretization.reportTiming = logical( ...
            parameters.execution.reportTiming);
    end
    if isfield(parameters.execution,'nonlinearTemporalOversampling')
        discretization.nonlinearTemporalOversampling = ...
            parameters.execution.nonlinearTemporalOversampling;
    end
    if isfield(parameters.execution,'adaptiveDirectionalSteps')
        discretization.adaptiveDirectionalSteps = logical( ...
            parameters.execution.adaptiveDirectionalSteps);
    end
    if isfield(parameters.execution,'quadraticStepMultipliers')
        discretization.quadraticStepMultipliers = ...
            parameters.execution.quadraticStepMultipliers;
    end
    if isfield(parameters.execution,'cubicStepMultipliers')
        discretization.cubicStepMultipliers = ...
            parameters.execution.cubicStepMultipliers;
    end
end
validateattributes(discretization.nonlinearTemporalOversampling, ...
    {'numeric'},{'scalar','real','>=',1,'finite'});
validateattributes(discretization.adaptiveDirectionalSteps, ...
    {'logical'},{'scalar'});
discretization.quadraticStepMultipliers = validate_step_multipliers( ...
    discretization.quadraticStepMultipliers,'quadraticStepMultipliers');
discretization.cubicStepMultipliers = validate_step_multipliers( ...
    discretization.cubicStepMultipliers,'cubicStepMultipliers');
if ~discretization.adaptiveDirectionalSteps
    discretization.quadraticStepMultipliers = 1;
    discretization.cubicStepMultipliers = 1;
end
end

function values = validate_step_multipliers(values,name)
validateattributes(values,{'numeric'}, ...
    {'vector','real','positive','finite','nonempty'});
values = unique(values(:).','stable');
if any(diff(values) <= 0)
    error('vi_cylinder_wnl_discretization:StepMultipliers', ...
        '%s must be strictly increasing.',name);
end
if numel(values) > 1 && numel(values) < 3
    error('vi_cylinder_wnl_discretization:StepMultiplierCount', ...
        ['%s needs at least three entries for a Richardson agreement ', ...
         'test, or exactly one entry to disable step selection.'],name);
end
if numel(values) > 1 && any(abs( ...
        values(2:end)./values(1:end-1)-2) > 100*eps)
    error('vi_cylinder_wnl_discretization:StepMultiplierRatio', ...
        '%s must form a factor-of-two sequence.',name);
end
end

function index = gauge_vertical_index(grid, targetCoordinate)
% Keep the pressure gauge away from physical and element boundaries. This
% prevents it from replacing one of the multi-domain matching conditions.
excluded = unique([1; grid.numberOfPoints; grid.interfacePairs(:)]);
candidates = setdiff((1:grid.numberOfPoints).', excluded);
if isempty(candidates)
    error('vi_cylinder_wnl_discretization:NoGaugePoint', ...
        'The lower vertical grid has no interior non-interface gauge point.');
end
[~, nearest] = min(abs(grid.z(candidates)-targetCoordinate));
index = candidates(nearest);
end
