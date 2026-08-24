function action = vi_cylinder_wnl_nonlinear_action(order, ...
    discretization, vectors, specs, specOut, fourierShift, frequencies)
%VI_CYLINDER_WNL_NONLINEAR_ACTION Directional C or D cylinder action.
%
% The nonlinear remainder of the exact fixed-domain ALE residual is
% evaluated for complex Fourier directions and differentiated by symmetric
% polarization. Terms whose second and third Frechet derivatives are
% analytically zero are removed before differencing; this avoids a
% roundoff/step residual floor without changing C or D. The retained
% remainder contains bulk convection and metric terms, nonlinear
% incompressibility and kinematics, nonlinear surface traction geometry,
% capillarity, and vibration acting on the tilted interface normal.

validateattributes(order, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2, '<=', 3});
if ~iscell(vectors)
    vectors = num2cell(vectors, 1);
end
if ~iscell(specs)
    specs = num2cell(specs);
end
assert(numel(vectors) == order && numel(specs) == order, ...
    'There must be one vector and one specification per direction.');
if all(cellfun(@(value) isstruct(value) && ...
        isfield(value, 'coeff') && isfield(value, 'spec'), vectors))
    action = field_level_action(order, discretization, vectors, specOut);
    return;
end
frequencies = frequencies(:).';
assert(numel(frequencies) == order, ...
    'There must be one temporal frequency per direction.');

if order == 2
    step = discretization.quadraticStep;
    rpp = residual_for_weights(discretization, vectors, specs, ...
        step*[1, 1], frequencies, specOut, fourierShift);
    rpm = residual_for_weights(discretization, vectors, specs, ...
        step*[1, -1], frequencies, specOut, fourierShift);
    rmp = residual_for_weights(discretization, vectors, specs, ...
        step*[-1, 1], frequencies, specOut, fourierShift);
    rmm = residual_for_weights(discretization, vectors, specs, ...
        step*[-1, -1], frequencies, specOut, fourierShift);
    secondDerivative = (rpp-rpm-rmp+rmm)/(4*step^2);
    action = -0.5*secondDerivative;
else
    step = discretization.cubicStep;
    thirdDerivative = complex(zeros(discretization.ndof, 1));
    signValues = [-1, 1];
    for sx = signValues
        for sy = signValues
            for sz = signValues
                signs = [sx, sy, sz];
                residual = residual_for_weights(discretization, ...
                    vectors, specs, step*signs, frequencies, ...
                    specOut, fourierShift);
                thirdDerivative = thirdDerivative + ...
                    prod(signs)*residual;
            end
        end
    end
    thirdDerivative = thirdDerivative/(8*step^3);
    action = -thirdDerivative/6;
end
end


function action = field_level_action(order, d, fields, specOut)
actionTimer = tic;
maximumFrequency = 0;
for fieldIndex = 1:numel(fields)
    frequency = fields{fieldIndex}.spec.n + fields{fieldIndex}.spec.s;
    maximumFrequency = maximumFrequency + max(abs(frequency));
end
maximumFrequency = maximumFrequency+1; % vibration shift
ntau = max(32,2*ceil( ...
    d.nonlinearTemporalOversampling*maximumFrequency)+4);
if mod(ntau, 2) ~= 0
    ntau = ntau+1;
end
tauValues = (0:ntau-1)*(2*pi/ntau);

if order == 2
    baseStep = d.quadraticStep;
    stepMultipliers = d.quadraticStepMultipliers;
else
    baseStep = d.cubicStep;
    stepMultipliers = d.cubicStepMultipliers;
end

useParallel = d.useParallelNonlinearActions && ...
    license('test','Distrib_Computing_Toolbox') && ...
    ~isempty(gcp('nocreate'));
rawActions = cell(numel(stepMultipliers),1);
for stepIndex = 1:numel(stepMultipliers)
    step = baseStep*stepMultipliers(stepIndex);
    [rawActions{stepIndex},parallelStepCompleted] = ...
        field_level_action_at_step( ...
        order,d,fields,specOut,step,tauValues,useParallel);
    if useParallel && ~parallelStepCompleted
        % A process pool can expire during a long preceding recovery or a
        % worker can abort under memory pressure. The failed step was
        % recomputed serially; keep later steps serial as well.
        useParallel = false;
    end
end
[action,stepDiagnostics] = select_directional_action( ...
    rawActions,stepMultipliers);
if d.reportTiming
    fprintf(['WNL nonlinear action order %d, m=%+d: %d temporal ', ...
        'samples x %d step(s), parallel=%d, %.3f s\n'], ...
        order,specOut.m,ntau,numel(stepMultipliers),useParallel, ...
        toc(actionTimer));
    if stepDiagnostics.adaptive
        fprintf(['  directional-step Richardson selection: base %.3e, ', ...
            'multipliers %s, selected range [%g,%g], estimated ', ...
            'relative disagreement %.3e\n'],baseStep, ...
            mat2str(stepMultipliers), ...
            stepDiagnostics.selectedMultiplierRange(1), ...
            stepDiagnostics.selectedMultiplierRange(2), ...
            stepDiagnostics.relativeDisagreement);
    end
end
end

function [action,parallelCompleted] = field_level_action_at_step( ...
        order,d,fields,specOut,step,tauValues,useParallel)
ntau = numel(tauValues);
action = complex(zeros(d.ndof,numel(specOut.n)));
parallelCompleted = false;
if useParallel
    % Every temporal collocation sample is independent. Matrix-valued
    % reduction avoids storing an ndof-by-nt-by-ntau temporary array.
    try
        parfor tauIndex = 1:ntau
            contribution = field_level_tau_contribution(order,d,fields, ...
                specOut,step,tauValues(tauIndex),ntau);
            action = action+contribution;
        end
        parallelCompleted = true;
        return;
    catch parallelError
        warning('vi_cylinder_wnl_nonlinear_action:ParallelFallback', ...
            ['Parallel nonlinear sampling failed (%s). The current ', ...
             'directional step will be recomputed serially and the ', ...
             'remaining steps will stay serial.'],parallelError.message);
        action = complex(zeros(d.ndof,numel(specOut.n)));
    end
end
for tauIndex = 1:ntau
    action = action+field_level_tau_contribution(order,d,fields, ...
        specOut,step,tauValues(tauIndex),ntau);
end
end

function [action,information] = select_directional_action( ...
        rawActions,stepMultipliers)
numberOfSteps = numel(rawActions);
information = struct('adaptive',numberOfSteps > 1, ...
    'selectedRichardsonPair',0,'selectedMultiplierRange', ...
    [stepMultipliers(1),stepMultipliers(1)], ...
    'relativeDisagreement',NaN,'richardsonDifferences',zeros(0,1));
if numberOfSteps == 1
    action = rawActions{1};
    return;
end

% Every symmetric mixed second- or third-derivative stencil has an O(h^2)
% leading truncation error. Adjacent steps differ by a factor of two, so
% (4*D(h)-D(2h))/3 removes that term. Agreement between neighboring
% Richardson estimates identifies the range in which neither small-step
% cancellation nor large-step truncation dominates.
richardson = cell(numberOfSteps-1,1);
for stepIndex = 1:numberOfSteps-1
    richardson{stepIndex} = ...
        (4*rawActions{stepIndex}-rawActions{stepIndex+1})/3;
end
differences = inf(numberOfSteps-2,1);
for pairIndex = 1:numberOfSteps-2
    first = richardson{pairIndex};
    second = richardson{pairIndex+1};
    if directional_action_is_finite(first) && ...
            directional_action_is_finite(second)
        differences(pairIndex) = norm(first-second,'fro') / ...
            max([norm(first,'fro'),norm(second,'fro'),eps]);
    end
end
[bestDifference,bestPair] = min(differences);
if isempty(bestDifference) || ~isfinite(bestDifference)
    error('vi_cylinder_wnl_nonlinear_action:DirectionalStepFailure', ...
        ['No finite neighboring Richardson estimates were obtained. ', ...
         'Inspect the ALE step multipliers and mode normalization.']);
end
action = 0.5*(richardson{bestPair}+richardson{bestPair+1});
information.selectedRichardsonPair = bestPair;
information.selectedMultiplierRange = ...
    [stepMultipliers(bestPair),stepMultipliers(bestPair+2)];
information.relativeDisagreement = bestDifference;
information.richardsonDifferences = differences;
end

function tf = directional_action_is_finite(value)
tf = all(isfinite(real(value(:)))) && ...
    all(isfinite(imag(value(:))));
end

function contribution = field_level_tau_contribution( ...
        order,d,fields,specOut,step,tau,ntau)
spatialVectors = cell(order,1);
timeVectors = cell(order,1);
for fieldIndex = 1:order
    frequency = fields{fieldIndex}.spec.n(:)+ ...
        fields{fieldIndex}.spec.s;
    phase = exp(1i*frequency*tau);
    spatialVectors{fieldIndex} = fields{fieldIndex}.coeff*phase;
    timeVectors{fieldIndex} = fields{fieldIndex}.coeff * ...
        ((wnl_spec_lambda(fields{fieldIndex}.spec)+ ...
        1i*d.parameters.omegaStar*frequency).*phase);
end
gravity = d.parameters.g_sgn+d.parameters.aCritical * ...
    cos(tau+d.parameters.phase);
if order == 2
    rpp = snapshot_residual(d,spatialVectors,timeVectors,fields, ...
        step*[1,1],specOut.m,gravity);
    rpm = snapshot_residual(d,spatialVectors,timeVectors,fields, ...
        step*[1,-1],specOut.m,gravity);
    rmp = snapshot_residual(d,spatialVectors,timeVectors,fields, ...
        step*[-1,1],specOut.m,gravity);
    rmm = snapshot_residual(d,spatialVectors,timeVectors,fields, ...
        step*[-1,-1],specOut.m,gravity);
    nonlinearSnapshot = -0.5*(rpp-rpm-rmp+rmm)/(4*step^2);
else
    thirdDerivative = complex(zeros(d.ndof,1));
    signValues = [-1,1];
    for sx = signValues
        for sy = signValues
            for sz = signValues
                signs = [sx,sy,sz];
                value = snapshot_residual(d,spatialVectors,timeVectors, ...
                    fields,step*signs,specOut.m,gravity);
                thirdDerivative = thirdDerivative+prod(signs)*value;
            end
        end
    end
    nonlinearSnapshot = -thirdDerivative/(48*step^3);
end
outputPhase = exp(-1i*(specOut.n+specOut.s)*tau)/ntau;
contribution = nonlinearSnapshot*outputPhase;
end

function residual = snapshot_residual(d, spatialVectors, timeVectors, ...
    fields, weights, mOut, gravity)
specs = cellfun(@(field) field.spec, fields, 'UniformOutput', false);
magnitudes = cellfun(@(s) abs(s.m), specs);
maximumMode = max([magnitudes(:); abs(mOut); 1]);
ntheta = max(d.ntheta, 4*maximumMode+8);
if mod(ntheta, 2) ~= 0
    ntheta = ntheta+1;
end
theta = reshape((0:ntheta-1)*(2*pi/ntheta), 1, 1, []);
[state, stateTime] = zero_physical_state(d, ntheta);
for directionIndex = 1:numel(spatialVectors)
    coefficient = unpack_coefficient(d, spatialVectors{directionIndex});
    timeCoefficient = unpack_coefficient(d, timeVectors{directionIndex});
    phase = exp(1i*specs{directionIndex}.m*theta);
    state = add_modal_state(state, coefficient, ...
        weights(directionIndex), phase);
    stateTime = add_modal_state(stateTime, timeCoefficient, ...
        weights(directionIndex), phase);
end
residual = complete_residual_coefficient(d, state, stateTime, ...
    theta, mOut, gravity, false);
end

function residualCoefficient = residual_for_weights(d, vectors, specs, ...
    weights, frequencies, specOut, fourierShift)
magnitudes = cellfun(@(s) abs(s.m), specs);
maximumMode = max([magnitudes(:); abs(specOut.m); 1]);
ntheta = max(d.ntheta, 4*maximumMode+8);
if mod(ntheta, 2) ~= 0
    ntheta = ntheta+1;
end
theta = reshape((0:ntheta-1)*(2*pi/ntheta), 1, 1, []);

[state, stateTime] = zero_physical_state(d, ntheta);
for directionIndex = 1:numel(vectors)
    coefficient = unpack_coefficient(d, vectors{directionIndex});
    phase = exp(1i*specs{directionIndex}.m*theta);
    timeFactor = wnl_spec_lambda(specs{directionIndex}) + ...
        1i*d.parameters.omegaStar*frequencies(directionIndex);
    state = add_modal_state(state, coefficient, ...
        weights(directionIndex), phase);
    stateTime = add_modal_state(stateTime, coefficient, ...
        weights(directionIndex)*timeFactor, phase);
end

if fourierShift == 0
    gravityCoefficient = d.parameters.g_sgn;
    residualCoefficient = complete_residual_coefficient( ...
        d, state, stateTime, theta, specOut.m, ...
        gravityCoefficient, false);
elseif abs(fourierShift) == 1
    gravityCoefficient = 0.5*d.parameters.aCritical * ...
        exp(1i*fourierShift*d.parameters.phase);
    residualCoefficient = complete_residual_coefficient( ...
        d, state, stateTime, theta, specOut.m, ...
        gravityCoefficient, true);
else
    residualCoefficient = complex(zeros(d.ndof, 1));
end
end

function residual = complete_residual_coefficient(d, state, stateTime, ...
    theta, mOut, gravityCoefficient, forcingOnly)
layout = d.layout;
ntheta = numel(theta);
if forcingOnly
    residual = complex(zeros(d.ndof, 1));
else
    [residualD, derivativesD] = layer_residual(d, state.d, ...
        stateTime.d, state.zeta, stateTime.zeta, 'd');
    [residualL, derivativesL] = layer_residual(d, state.l, ...
        stateTime.l, state.zeta, stateTime.zeta, 'l');
    coefficientD = project_layer(residualD, theta, mOut);
    coefficientL = project_layer(residualL, theta, mOut);
    residual = pack_layer_residuals(d, coefficientD, coefficientL);
end

% Exact interface geometry and traction are needed for both the autonomous
% residual and the +/-1 vibration Fourier coefficients.
zetaR = apply_r(state.zeta, d.Dr);
zetaTheta = theta_derivative(state.zeta);
zetaThetaOverR = divide_by_r(zetaTheta, d.r);
slopeMagnitudeSquared = zetaR.^2 + zetaThetaOverR.^2;
normalDenominator = sqrt(1+slopeMagnitudeSquared);
nR = -zetaR./normalDenominator;
nTheta = -zetaThetaOverR./normalDenominator;
nZMinusOne = -slopeMagnitudeSquared ./ ...
    (normalDenominator.*(1+normalDenominator));

if forcingOnly
    scalarNormalForce = d.hAtwood*gravityCoefficient*state.zeta;
    tractionR = scalarNormalForce.*nR;
    tractionTheta = scalarNormalForce.*nTheta;
    % Remove the flat-interface linear vibration term exactly. Only the
    % tilted-normal correction belongs to C/D.
    tractionZ = scalarNormalForce.*nZMinusOne;
    kinematic = complex(zeros(size(state.zeta)));
else
    interfaceD = interface_values(d, state.d, derivativesD, 'd');
    interfaceL = interface_values(d, state.l, derivativesL, 'l');
    tractionD = stress_times_normal_remainder(d,interfaceD,state.d, ...
        d.DzD,derivativesD.geometry,d.layout.nzD, ...
        nR,nTheta,nZMinusOne,d.muD);
    tractionL = stress_times_normal_remainder(d,interfaceL,state.l, ...
        d.DzL,derivativesL.geometry,1,nR,nTheta,nZMinusOne,d.muL);

    surfaceFactor = sqrt(1+slopeMagnitudeSquared);
    linearCurvature = apply_r(zetaR,d.Dr) + ...
        divide_by_r(zetaR,d.r) + ...
        divide_by_r(theta_derivative(zetaThetaOverR),d.r);
    inverseSurfaceRemainder = -slopeMagnitudeSquared ./ ...
        (surfaceFactor.*(1+surfaceFactor));
    curvatureRadialRemainder = zetaR.*inverseSurfaceRemainder;
    curvatureThetaRemainder = ...
        zetaThetaOverR.*inverseSurfaceRemainder;
    curvatureRemainder = ...
        apply_r(curvatureRadialRemainder,d.Dr) + ...
        divide_by_r(curvatureRadialRemainder,d.r) + ...
        divide_by_r(theta_derivative(curvatureThetaRemainder),d.r);
    linearScalarNormalForce = ...
        d.hAtwood*gravityCoefficient*state.zeta - ...
        linearCurvature/d.parameters.Bd;
    nonlinearScalarNormalForce = ...
        -curvatureRemainder/d.parameters.Bd;
    scalarNormalForce = ...
        linearScalarNormalForce+nonlinearScalarNormalForce;
    tractionR = tractionD.r-tractionL.r + ...
        scalarNormalForce.*nR;
    tractionTheta = tractionD.theta-tractionL.theta + ...
        scalarNormalForce.*nTheta;
    tractionZ = tractionD.z-tractionL.z + ...
        scalarNormalForce.*nZMinusOne + ...
        nonlinearScalarNormalForce;
    % zeta_t-w is the exactly linear part of the kinematic condition and
    % therefore has zero C/D action. Evaluate only the nonlinear geometric
    % transport here; retaining the linear part in a small-step stencil
    % contaminates the mean/combination forcing through cancellation.
    kinematic = interfaceD.ur.*zetaR + ...
        interfaceD.ut.*zetaThetaOverR;
end

tractionRCoefficient = fourier_project(tractionR, theta, mOut);
tractionThetaCoefficient = fourier_project( ...
    tractionTheta, theta, mOut);
tractionZCoefficient = fourier_project(tractionZ, theta, mOut);
kinematicCoefficient = fourier_project(kinematic, theta, mOut);

% Linear boundary rows have exactly zero quadratic and cubic derivatives.
% Set those rows to zero before replacing the interior interface rows.
residual = zero_linear_boundary_rows(residual, d, mOut);
for ir = 2:layout.nr-1
    nodeD = grid_index(ir, layout.nzD, layout.nr);
    nodeL = grid_index(ir, 1, layout.nr);
    % The ALE map places both fluids on the same fixed interface, so all
    % three velocity-continuity equations are exactly linear. Their
    % quadratic and cubic actions are identically zero—not small finite-
    % difference numbers.
    residual(layout.d.ur(nodeD)) = 0;
    residual(layout.d.ut(nodeD)) = 0;
    residual(layout.d.w(nodeD)) = 0;
    residual(layout.l.ur(nodeL)) = tractionRCoefficient(ir);
    residual(layout.l.ut(nodeL)) = tractionThetaCoefficient(ir);
    residual(layout.l.w(nodeL)) = tractionZCoefficient(ir);
    if ~(mOut == 0 && ir == volume_constraint_index(layout.nr))
        residual(layout.zeta(ir)) = kinematicCoefficient(ir);
    end
end

% The arrays above retain a singleton vertical dimension. Make the output
% shape explicit and guard against accidental theta-grid leakage.
residual = reshape(residual, d.ndof, 1);
assert(numel(residual) == d.ndof && ntheta >= 1);
end

function [residual, derivatives] = layer_residual(d, state, stateTime, ...
    zeta, zetaTime, layerName)
if strcmp(layerName, 'd')
    Dz = d.DzD;
    z = d.zD;
    chi = 1+z;
    chiPrime = 1;
    rho = d.rhoD;
    viscosity = d.muD*d.parameters.C;
else
    Dz = d.DzL;
    z = d.zL;
    chi = 1-z;
    chiPrime = -1;
    rho = d.rhoL;
    viscosity = d.muL*d.parameters.C;
end

geometry = make_geometry(d, Dz, chi, chiPrime, zeta, zetaTime);
derivatives = velocity_derivatives(d, Dz, geometry, state);
urTimeRemainder = nonlinear_time_derivative( ...
    Dz,geometry,state.ur);
utTimeRemainder = nonlinear_time_derivative( ...
    Dz,geometry,state.ut);
wTimeRemainder = nonlinear_time_derivative( ...
    Dz,geometry,state.w);

urConvective = state.ur.*derivatives.urR + ...
    state.ut.*derivatives.urTheta + ...
    state.w.*derivatives.urZ - divide_by_r(state.ut.^2, d.r);
utConvective = state.ur.*derivatives.utR + ...
    state.ut.*derivatives.utTheta + ...
    state.w.*derivatives.utZ + ...
    divide_by_r(state.ur.*state.ut, d.r);
wConvective = state.ur.*derivatives.wR + ...
    state.ut.*derivatives.wTheta + ...
    state.w.*derivatives.wZ;

lapUrRemainder = nonlinear_scalar_laplacian( ...
    d,Dz,geometry,state.ur) - 2*divide_by_r( ...
    nonlinear_theta_derivative(d,Dz,geometry,state.ut),d.r);
lapUtRemainder = nonlinear_scalar_laplacian( ...
    d,Dz,geometry,state.ut) + 2*divide_by_r( ...
    nonlinear_theta_derivative(d,Dz,geometry,state.ur),d.r);
lapWRemainder = nonlinear_scalar_laplacian( ...
    d,Dz,geometry,state.w);

pressureRemainderR = nonlinear_r_derivative( ...
    d,Dz,geometry,state.p);
pressureRemainderTheta = nonlinear_theta_derivative( ...
    d,Dz,geometry,state.p);
pressureRemainderZ = nonlinear_z_derivative(Dz,geometry,state.p);

% Return the nonlinear ALE remainder, not the complete residual.  The
% flat-interface time, pressure, viscous, and incompressibility terms are
% exactly linear and have zero C/D derivative. The helpers below construct
% operator differences analytically rather than subtracting two large
% derivative arrays, preventing a second cancellation inside this remainder.
residual.ur = rho*(urTimeRemainder+urConvective) + ...
    pressureRemainderR-viscosity*lapUrRemainder;
residual.ut = rho*(utTimeRemainder+utConvective) + ...
    pressureRemainderTheta-viscosity*lapUtRemainder;
residual.w = rho*(wTimeRemainder+wConvective) + ...
    pressureRemainderZ-viscosity*lapWRemainder;
residual.p = nonlinear_r_derivative(d,Dz,geometry,state.ur) + ...
    nonlinear_theta_derivative(d,Dz,geometry,state.ut) + ...
    nonlinear_z_derivative(Dz,geometry,state.w);
derivatives.geometry = geometry;
end

function remainder = nonlinear_scalar_laplacian(d,Dz,geometry,field)
flatR = apply_r(field,d.Dr);
flatTheta = divide_by_r(theta_derivative(field),d.r);
flatZ = apply_z(field,Dz);
deltaR = nonlinear_r_derivative(d,Dz,geometry,field);
deltaTheta = nonlinear_theta_derivative(d,Dz,geometry,field);
deltaZ = nonlinear_z_derivative(Dz,geometry,field);
remainder = apply_r(deltaR,d.Dr) + ...
    nonlinear_r_derivative(d,Dz,geometry,flatR) + ...
    nonlinear_r_derivative(d,Dz,geometry,deltaR) + ...
    divide_by_r(deltaR,d.r) + ...
    divide_by_r(theta_derivative(deltaTheta),d.r) + ...
    nonlinear_theta_derivative(d,Dz,geometry,flatTheta) + ...
    nonlinear_theta_derivative(d,Dz,geometry,deltaTheta) + ...
    apply_z(deltaZ,Dz) + ...
    nonlinear_z_derivative(Dz,geometry,flatZ) + ...
    nonlinear_z_derivative(Dz,geometry,deltaZ);
end

function geometry = make_geometry(d, Dz, chi, chiPrime, ...
    zeta, zetaTime)
nr = d.layout.nr;
nz = numel(chi);
ntheta = size(zeta, 3);
zetaR = apply_r(zeta, d.Dr);
zetaTheta = theta_derivative(zeta);
zeta3 = reshape(zeta, nr, 1, ntheta);
zetaTime3 = reshape(zetaTime, nr, 1, ntheta);
zetaR3 = reshape(zetaR, nr, 1, ntheta);
zetaTheta3 = reshape(zetaTheta, nr, 1, ntheta);
chi3 = reshape(chi, 1, nz, 1);
J = 1+chiPrime*zeta3;
if any(abs(J(:)) < 0.25)
    error('vi_cylinder_wnl_nonlinear_action:InvalidALEMap', ...
        'The directional perturbation made the ALE Jacobian too small.');
end
geometry.Dz = Dz;
geometry.chi3 = chi3;
geometry.J = J;
geometry.jacobianPerturbation = chiPrime*zeta3;
geometry.zetaTime3 = zetaTime3;
geometry.zetaR3 = zetaR3;
geometry.zetaTheta3 = zetaTheta3;
end

function derivatives = velocity_derivatives(d, Dz, geometry, state)
derivatives.urR = physical_r_derivative(d, Dz, geometry, state.ur);
derivatives.urTheta = physical_theta_derivative( ...
    d, Dz, geometry, state.ur);
derivatives.urZ = physical_z_derivative(Dz, geometry, state.ur);
derivatives.utR = physical_r_derivative(d, Dz, geometry, state.ut);
derivatives.utTheta = physical_theta_derivative( ...
    d, Dz, geometry, state.ut);
derivatives.utZ = physical_z_derivative(Dz, geometry, state.ut);
derivatives.wR = physical_r_derivative(d, Dz, geometry, state.w);
derivatives.wTheta = physical_theta_derivative( ...
    d, Dz, geometry, state.w);
derivatives.wZ = physical_z_derivative(Dz, geometry, state.w);
end

function value = physical_r_derivative(d, Dz, geometry, field)
fieldZ = apply_z(field, Dz);
value = apply_r(field, d.Dr) - ...
    geometry.chi3.*geometry.zetaR3./geometry.J.*fieldZ;
end

function value = physical_theta_derivative(d, Dz, geometry, field)
fieldZ = apply_z(field, Dz);
value = divide_by_r(theta_derivative(field) - ...
    geometry.chi3.*geometry.zetaTheta3./geometry.J.*fieldZ, d.r);
end

function value = physical_z_derivative(Dz, geometry, field)
value = apply_z(field, Dz)./geometry.J;
end

function value = nonlinear_r_derivative(~,Dz,geometry,field)
fieldZ = apply_z(field, Dz);
value = -geometry.chi3.*geometry.zetaR3./geometry.J.*fieldZ;
end

function value = nonlinear_theta_derivative(d,Dz,geometry,field)
fieldZ = apply_z(field,Dz);
value = divide_by_r(-geometry.chi3.*geometry.zetaTheta3 ./ ...
    geometry.J.*fieldZ,d.r);
end

function value = nonlinear_z_derivative(Dz,geometry,field)
fieldZ = apply_z(field,Dz);
value = -geometry.jacobianPerturbation./geometry.J.*fieldZ;
end

function value = nonlinear_time_derivative(Dz,geometry,field)
fieldZ = apply_z(field,Dz);
value = -geometry.chi3.*geometry.zetaTime3./geometry.J.*fieldZ;
end

function interface = interface_values(d, state, derivatives, layerName)
if strcmp(layerName, 'd')
    iz = d.layout.nzD;
else
    iz = 1;
end
interface.ur = state.ur(:, iz, :);
interface.ut = state.ut(:, iz, :);
interface.w = state.w(:, iz, :);
interface.p = state.p(:, iz, :);
interface.urR = derivatives.urR(:, iz, :);
interface.urTheta = derivatives.urTheta(:, iz, :);
interface.urZ = derivatives.urZ(:, iz, :);
interface.utR = derivatives.utR(:, iz, :);
interface.utTheta = derivatives.utTheta(:, iz, :);
interface.utZ = derivatives.utZ(:, iz, :);
interface.wR = derivatives.wR(:, iz, :);
interface.wTheta = derivatives.wTheta(:, iz, :);
interface.wZ = derivatives.wZ(:, iz, :);
end

function traction = stress_times_normal(d, interface, ...
    nR, nTheta, nZ, viscosityRatio)
viscosity = viscosityRatio*d.parameters.C;
Trr = -interface.p + 2*viscosity*interface.urR;
TthetaTheta = -interface.p + 2*viscosity*( ...
    interface.utTheta+divide_by_r(interface.ur, d.r));
Tzz = -interface.p + 2*viscosity*interface.wZ;
TrTheta = viscosity*(interface.urTheta + ...
    interface.utR-divide_by_r(interface.ut, d.r));
TrZ = viscosity*(interface.urZ+interface.wR);
TthetaZ = viscosity*(interface.utZ+interface.wTheta);
traction.r = Trr.*nR + TrTheta.*nTheta + TrZ.*nZ;
traction.theta = TrTheta.*nR + TthetaTheta.*nTheta + ...
    TthetaZ.*nZ;
traction.z = TrZ.*nR + TthetaZ.*nTheta + Tzz.*nZ;
end

function traction = stress_times_normal_remainder( ...
        d,interface,state,Dz,geometry,iz,nR,nTheta,nZMinusOne, ...
        viscosityRatio)
% Evaluate T*n-(T_flat*e_z) without subtracting two complete tractions.
% T*(n-e_z) is already nonlinear. The remaining vertical-column correction
% contains only metric derivative remainders, which are formed explicitly.
traction = stress_times_normal( ...
    d,interface,nR,nTheta,nZMinusOne,viscosityRatio);
viscosity = viscosityRatio*d.parameters.C;
deltaUrZ = nonlinear_z_derivative(Dz,geometry,state.ur);
deltaUtZ = nonlinear_z_derivative(Dz,geometry,state.ut);
deltaWZ = nonlinear_z_derivative(Dz,geometry,state.w);
deltaWR = nonlinear_r_derivative(d,Dz,geometry,state.w);
deltaWTheta = nonlinear_theta_derivative(d,Dz,geometry,state.w);
traction.r = traction.r + viscosity*( ...
    deltaUrZ(:,iz,:)+deltaWR(:,iz,:));
traction.theta = traction.theta + viscosity*( ...
    deltaUtZ(:,iz,:)+deltaWTheta(:,iz,:));
traction.z = traction.z + 2*viscosity*deltaWZ(:,iz,:);
end

function residual = pack_layer_residuals(d, dense, light)
layout = d.layout;
residual = complex(zeros(d.ndof, 1));
residual(layout.d.ur) = dense.ur(:);
residual(layout.d.ut) = dense.ut(:);
residual(layout.d.w) = dense.w(:);
residual(layout.d.p) = dense.p(:);
residual(layout.l.ur) = light.ur(:);
residual(layout.l.ut) = light.ut(:);
residual(layout.l.w) = light.w(:);
residual(layout.l.p) = light.p(:);
end

function coefficient = project_layer(residual, theta, mOut)
coefficient.ur = fourier_project(residual.ur, theta, mOut);
coefficient.ut = fourier_project(residual.ut, theta, mOut);
coefficient.w = fourier_project(residual.w, theta, mOut);
coefficient.p = fourier_project(residual.p, theta, mOut);
end

function coefficient = fourier_project(field, theta, mOut)
phase = exp(-1i*mOut*theta);
coefficient = sum(field.*phase, 3)/numel(theta);
end

function residual = zero_linear_boundary_rows(residual, d, mOut)
layout = d.layout;
nr = layout.nr;
% Horizontal rigid walls, excluding radial corners.
for ir = 2:nr-1
    nodeD = grid_index(ir, 1, nr);
    nodeL = grid_index(ir, layout.nzL, nr);
    residual([layout.d.ur(nodeD), layout.d.ut(nodeD), ...
        layout.d.w(nodeD)]) = 0;
    residual([layout.l.ur(nodeL), layout.l.ut(nodeL), ...
        layout.l.w(nodeL)]) = 0;
end
% Axis has four linear regularity rows; the sidewall has three linear
% velocity rows while its divergence row remains active.
for iz = 1:layout.nzD
    axis = grid_index(1, iz, nr);
    side = grid_index(nr, iz, nr);
    residual([layout.d.ur(axis), layout.d.ut(axis), ...
        layout.d.w(axis), layout.d.p(axis)]) = 0;
    residual([layout.d.ur(side), layout.d.ut(side), ...
        layout.d.w(side)]) = 0;
end
for iz = 1:layout.nzL
    axis = grid_index(1, iz, nr);
    side = grid_index(nr, iz, nr);
    residual([layout.l.ur(axis), layout.l.ut(axis), ...
        layout.l.w(axis), layout.l.p(axis)]) = 0;
    residual([layout.l.ur(side), layout.l.ut(side), ...
        layout.l.w(side)]) = 0;
end

% Artificial vertical interfaces carry only linear same-fluid matching
% conditions. Zero the C0 velocity/pressure rows on the left copy and the
% C1 velocity rows on the right copy. The right-copy pressure row remains
% the nonlinear incompressibility equation, exactly as in the linear
% operator.
residual = zero_vertical_matching_rows( ...
    residual, layout.d, d.verticalD.interfacePairs, nr);
residual = zero_vertical_matching_rows( ...
    residual, layout.l, d.verticalL.interfacePairs, nr);
residual(layout.zeta([1, nr])) = 0;
if mOut == 0
    residual(layout.zeta(volume_constraint_index(nr))) = 0;
    gaugeIr = max(2, min(nr-1, ceil(nr/2)));
    gaugeIz = d.denseGaugeVerticalIndex;
    gauge = grid_index(gaugeIr, gaugeIz, nr);
    residual(layout.d.p(gauge)) = 0;
end
end

function residual = zero_vertical_matching_rows( ...
    residual, layer, interfacePairs, nr)
for interfaceIndex = 1:size(interfacePairs, 1)
    leftIz = interfacePairs(interfaceIndex, 1);
    rightIz = interfacePairs(interfaceIndex, 2);
    for ir = 2:nr-1
        leftNode = grid_index(ir, leftIz, nr);
        rightNode = grid_index(ir, rightIz, nr);
        residual([layer.ur(leftNode), layer.ut(leftNode), ...
            layer.w(leftNode), layer.p(leftNode), ...
            layer.ur(rightNode), layer.ut(rightNode), ...
            layer.w(rightNode)]) = 0;
    end
end
end

function [state, stateTime] = zero_physical_state(d, ntheta)
nr = d.layout.nr;
state.d = zero_layer(nr, d.layout.nzD, ntheta);
state.l = zero_layer(nr, d.layout.nzL, ntheta);
state.zeta = complex(zeros(nr, 1, ntheta));
stateTime = state;
end

function layer = zero_layer(nr, nz, ntheta)
zero = complex(zeros(nr, nz, ntheta));
layer.ur = zero;
layer.ut = zero;
layer.w = zero;
layer.p = zero;
end

function coefficient = unpack_coefficient(d, vector)
layout = d.layout;
vector = vector(:);
if numel(vector) ~= d.ndof
    error('vi_cylinder_wnl_nonlinear_action:StateSize', ...
        'Expected %d spatial unknowns, received %d.', ...
        d.ndof, numel(vector));
end
coefficient.d.ur = reshape(vector(layout.d.ur), ...
    layout.nr, layout.nzD);
coefficient.d.ut = reshape(vector(layout.d.ut), ...
    layout.nr, layout.nzD);
coefficient.d.w = reshape(vector(layout.d.w), ...
    layout.nr, layout.nzD);
coefficient.d.p = reshape(vector(layout.d.p), ...
    layout.nr, layout.nzD);
coefficient.l.ur = reshape(vector(layout.l.ur), ...
    layout.nr, layout.nzL);
coefficient.l.ut = reshape(vector(layout.l.ut), ...
    layout.nr, layout.nzL);
coefficient.l.w = reshape(vector(layout.l.w), ...
    layout.nr, layout.nzL);
coefficient.l.p = reshape(vector(layout.l.p), ...
    layout.nr, layout.nzL);
coefficient.zeta = reshape(vector(layout.zeta), layout.nr, 1);
end

function state = add_modal_state(state, coefficient, weight, phase)
fields = {'ur', 'ut', 'w', 'p'};
for fieldIndex = 1:numel(fields)
    name = fields{fieldIndex};
    state.d.(name) = state.d.(name) + ...
        weight*coefficient.d.(name).*phase;
    state.l.(name) = state.l.(name) + ...
        weight*coefficient.l.(name).*phase;
end
state.zeta = state.zeta + weight*coefficient.zeta.*phase;
end

function derivative = apply_r(field, Dr)
fieldSize = size(field);
derivative = reshape(Dr*reshape(field, fieldSize(1), []), fieldSize);
end

function derivative = apply_z(field, Dz)
fieldSize = size(field);
permuted = permute(field, [2, 1, 3]);
permuted = reshape(permuted, fieldSize(2), []);
permuted = Dz*permuted;
derivative = ipermute(reshape(permuted, ...
    fieldSize(2), fieldSize(1), fieldSize(3)), [2, 1, 3]);
end

function derivative = theta_derivative(field)
ntheta = size(field, 3);
if mod(ntheta, 2) == 0
    modes = [0:(ntheta/2-1), -ntheta/2:-1];
else
    half = (ntheta-1)/2;
    modes = [0:half, -half:-1];
end
modes = reshape(modes, 1, 1, []);
derivative = ifft(1i*modes.*fft(field, [], 3), [], 3);
end

function quotient = divide_by_r(field, r)
% Regularized cylindrical division. Values at r=0 are obtained by
% polynomial extrapolation of the already-divided interior values. This is
% the appropriate l'Hopital limit for regular Fourier fields and prevents
% an arbitrary axis value from contaminating Chebyshev differentiation.
quotient = complex(zeros(size(field), 'like', field));
shape = size(field);
interior = reshape(field(2:end, :, :), numel(r)-1, []);
quotientInterior = interior./r(2:end);
quotient(2:end, :, :) = reshape(quotientInterior, ...
    [numel(r)-1, shape(2:end)]);
sampleIndices = 2:min(4, numel(r));
sampleR = r(sampleIndices);
weights = ones(numel(sampleIndices), 1);
for j = 1:numel(sampleIndices)
    for k = 1:numel(sampleIndices)
        if k ~= j
            weights(j) = weights(j) * ...
                (-sampleR(k))/(sampleR(j)-sampleR(k));
        end
    end
end
sampleValues = reshape(quotient(sampleIndices, :, :), ...
    numel(sampleIndices), []);
axisValue = weights.'*sampleValues;
quotient(1, :, :) = reshape(axisValue, [1, shape(2:end)]);
end

function node = grid_index(ir, iz, nr)
node = ir+(iz-1)*nr;
end

function index = volume_constraint_index(nr)
index = max(2, min(nr-1, ceil(nr/2)));
end
