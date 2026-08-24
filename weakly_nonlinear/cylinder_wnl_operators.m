function [operators, metadata] = cylinder_wnl_operators(parameters)
%CYLINDER_WNL_OPERATORS Physical WNL operators for the two-fluid cylinder.
%
% This factory discretizes the same model used by the reduced linear code:
% two viscous layers of nondimensional depth one, no slip at z=+/-1,
% legacy free-slide sidewall velocity conditions, full velocity and stress
% continuity at z=zeta, and either a free (d_r zeta=0) or pinned (zeta=0)
% contact line. One or more matched Chebyshev elements may be used in each
% vertical layer. The primitive variables and incompressibility constraint
% are retained, so bulk and boundary nonlinear forcing enter the adjoint
% solvability projection.

discretization = vi_cylinder_wnl_discretization(parameters);
blockCache = containers.Map('KeyType', 'char', 'ValueType', 'any');

operators = struct();
operators.B0 = @mass_callback;
operators.Lhat = @linear_fourier_callback;
operators.C = @quadratic_callback;
operators.D = @cubic_callback;
operators.CField = @quadratic_field_callback;
operators.DField = @cubic_field_callback;
operators.P = @detuning_callback;
operators.solveForced = @forced_solve_callback;
operators.conjugateAdjointRows = @conjugate_adjoint_rows_callback;
operators.nonlinearFourierShifts = [-1, 0, 1];

metadata = struct();
metadata.ndof = discretization.ndof;
metadata.discretization = discretization;
metadata.layout = discretization.layout;
metadata.contactLine = discretization.contactLine;
metadata.sidewallConditions = discretization.sidewallConditions;
metadata.radialGrid = discretization.radial;
metadata.verticalGrid.lower = discretization.verticalD;
metadata.verticalGrid.upper = discretization.verticalL;
metadata.normalizeDirect = @normalize_direct_mode;
metadata.zetaOverHPerUnitAmplitude = 1.0;
metadata.description = [ ...
    'Primitive-variable optionally Bessel-enriched radial / multi-domain ', ...
    'Chebyshev vertical / Fourier ALE discretization with the linear ', ...
    'solver''s wall, interface, and contact-line conditions. Artificial ', ...
    'vertical element boundaries enforce C0 primitive-variable and C1 ', ...
    'velocity matching.'];

    function B0 = mass_callback(spec)
        blocks = get_blocks(spec.m);
        B0 = blocks.B0;
    end

    function L = linear_fourier_callback(spec, k)
        blocks = get_blocks(spec.m);
        if k == 0
            L = blocks.L0;
        elseif k == 1
            L = blocks.Lplus;
        elseif k == -1
            L = blocks.Lminus;
        else
            L = sparse(discretization.ndof, discretization.ndof);
        end
    end

    function value = quadratic_callback(a, b, specA, specB, ...
            specOut, k, nuA, nuB)
        assert(specOut.m == specA.m+specB.m, ...
            'Quadratic azimuthal selection rule was violated.');
        value = vi_cylinder_wnl_nonlinear_action(2, ...
            discretization, {a, b}, {specA, specB}, specOut, ...
            k, [nuA, nuB]);
    end

    function value = cubic_callback(a, b, c, specA, specB, specC, ...
            specOut, k, nuA, nuB, nuC)
        assert(specOut.m == specA.m+specB.m+specC.m, ...
            'Cubic azimuthal selection rule was violated.');
        value = vi_cylinder_wnl_nonlinear_action(3, ...
            discretization, {a, b, c}, {specA, specB, specC}, ...
            specOut, k, [nuA, nuB, nuC]);
    end

    function value = quadratic_field_callback(fieldA, fieldB, specOut)
        value = vi_cylinder_wnl_nonlinear_action(2, ...
            discretization, {fieldA, fieldB}, ...
            {fieldA.spec, fieldB.spec}, specOut, [], []);
    end

    function value = cubic_field_callback(fieldA, fieldB, fieldC, specOut)
        value = vi_cylinder_wnl_nonlinear_action(3, ...
            discretization, {fieldA, fieldB, fieldC}, ...
            {fieldA.spec, fieldB.spec, fieldC.spec}, specOut, [], []);
    end

    function P = detuning_callback(spec, nOut, nIn, direction)
        blocks = get_blocks(spec.m);
        scale = detuning_scale(direction);
        k = nOut-nIn;
        if k == 1
            P = scale*blocks.dLdaPlus;
        elseif k == -1
            P = scale*blocks.dLdaMinus;
        else
            P = sparse(discretization.ndof, discretization.ndof);
        end
    end

    function result = forced_solve_callback(spec,forcing,opts)
        blocks = get_blocks(spec.m);
        result = vi_cylinder_wnl_forced_schur_solve( ...
            blocks,discretization.parameters.omegaStar, ...
            spec,forcing,opts);
    end

    function blocks = get_blocks(m)
        key = sprintf('%+d', m);
        if isKey(blockCache, key)
            blocks = blockCache(key);
            return;
        end
        [blocks.B0, blocks.L0, blocks.Lplus, blocks.Lminus, ...
            blocks.dLdaPlus, blocks.dLdaMinus, blocks.diagnostics] = ...
            vi_cylinder_wnl_linear_operators(discretization, m);
        blockCache(key) = blocks;
    end

    function normalized = normalize_direct_mode(field)
        zeta = field.coeff(discretization.layout.zeta, :);
        numberOfTimes = max(64, 8*numel(field.spec.n));
        tau = (0:numberOfTimes-1)*(2*pi/numberOfTimes);
        frequency = field.spec.n(:)+field.spec.s;
        temporal = exp(1i*frequency*tau);
        signal = zeta*temporal;
        [peak, peakIndex] = max(abs(signal(:)));
        if peak <= eps*max(norm(zeta), 1)
            error('cylinder_wnl_operators:ZeroInterfaceMode', ...
                ['The selected full-state mode has zero interface ', ...
                 'displacement. Supply a directSeed for the requested ', ...
                 '(m,radialIndex,s) branch and check the constrained ', ...
                 'full-residual reconstruction diagnostics.']);
        end
        phase = exp(-1i*angle(signal(peakIndex)));
        normalized = field;
        normalized.coeff = phase*field.coeff/peak;
    end

    function transformed = conjugate_adjoint_rows_callback( ...
            conjugatedField, sourceSpec)
        % Direct states use one primitive-variable ordering for m and -m.
        % Equation rows differ at r=0: the regularity combinations
        % u_r+i*u_theta and u_r-i*u_theta exchange their stored row
        % locations under complex conjugation. An algebraic adjoint lives
        % in equation-row space and therefore needs this extra permutation.
        if ~isstruct(conjugatedField) || ...
                ~isfield(conjugatedField,'coeff') || ...
                ~isfield(conjugatedField,'spec')
            error('cylinder_wnl_operators:ConjugateAdjointField', ...
                'The conjugated adjoint must be a WNL field structure.');
        end
        if conjugatedField.spec.m ~= -sourceSpec.m
            error('cylinder_wnl_operators:ConjugateAdjointMode', ...
                ['The conjugated adjoint has m=%d but the source mode ', ...
                 'requires m=%d.'],conjugatedField.spec.m,-sourceSpec.m);
        end
        transformed = conjugatedField;
        coefficients = transformed.coeff;
        coefficients = swap_axis_regularity_rows( ...
            coefficients,discretization.layout.d, ...
            discretization.layout.nzD,discretization.layout.nr);
        coefficients = swap_axis_regularity_rows( ...
            coefficients,discretization.layout.l, ...
            discretization.layout.nzL,discretization.layout.nr);
        transformed.coeff = coefficients;
    end
end

function coefficients = swap_axis_regularity_rows( ...
        coefficients,layer,nz,nr)
for verticalIndex = 1:nz
    axisNode = 1+(verticalIndex-1)*nr;
    firstRow = layer.ur(axisNode);
    secondRow = layer.ut(axisNode);
    savedFirst = coefficients(firstRow,:);
    coefficients(firstRow,:) = coefficients(secondRow,:);
    coefficients(secondRow,:) = savedFirst;
end
end

function scale = detuning_scale(direction)
if isnumeric(direction)
    validateattributes(direction, {'numeric'}, ...
        {'scalar', 'finite'});
    scale = direction;
elseif isstruct(direction) && isfield(direction, 'acceleration')
    scale = direction.acceleration;
    validateattributes(scale, {'numeric'}, {'scalar', 'finite'});
else
    error('cylinder_wnl_operators:DetuningDirection', ...
        ['The detuning direction must be numeric or a struct with an ', ...
         'acceleration field.']);
end
end
