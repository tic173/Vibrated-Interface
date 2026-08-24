function grid = vi_cylinder_radial_grid(parameters)
%VI_CYLINDER_RADIAL_GRID Chebyshev or Bessel-enriched radial operators.
%
% The Bessel-enriched option keeps Chebyshev--Lobatto evaluation nodes but
% replaces polynomial differentiation by differentiation in a problem-
% adapted radial space. The mandatory space contains, for every retained
% mode, J_m(beta*r) and the J_|m+/-1|(beta*r) functions carried by the
% horizontal velocity. Quadratic products and regular polynomials complete
% the square interpolation basis. This differentiates the selected linear
% Bessel branches analytically while preserving the existing pointwise ALE
% nonlinear residual.

numerics = parameters.numerics;
validateattributes(numerics.Nr, {'numeric'}, ...
    {'scalar', 'integer', '>=', 6});
nr = numerics.Nr;
[r, chebyshevD] = vi_chebyshev_lobatto(nr, [0, parameters.R0]);

options = default_options();
if isfield(numerics, 'radialGrid') && ~isempty(numerics.radialGrid)
    supplied = numerics.radialGrid;
    names = fieldnames(supplied);
    for nameIndex = 1:numel(names)
        options.(names{nameIndex}) = supplied.(names{nameIndex});
    end
end
requestedType = validatestring(options.type, ...
    {'chebyshev', 'besselEnriched'});

if strcmp(requestedType, 'chebyshev')
    grid = chebyshev_grid(r, chebyshevD, requestedType);
    return;
end

validateattributes(options.maximumProductOrder, {'numeric'}, ...
    {'scalar', 'integer', '>=', 1, '<=', 3});
validateattributes(options.maximumConditionNumber, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'});
validateattributes(options.independenceTolerance, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'});
validateattributes(options.fallbackToChebyshev, ...
    {'logical', 'numeric'}, {'scalar'});
if isnumeric(options.fallbackToChebyshev) && ...
        ~ismember(options.fallbackToChebyshev, [0, 1])
    error('vi_cylinder_radial_grid:BadFallbackFlag', ...
        'fallbackToChebyshev must be true/false or 1/0.');
end

try
    [values, first, second, labels, numberOfMandatory] = ...
        candidate_library(parameters, r, options.maximumProductOrder);
    if numberOfMandatory > nr
        error('vi_cylinder_radial_grid:TooFewRadialPoints', ...
            ['Nr=%d is smaller than the %d mandatory constant/Bessel ', ...
             'functions. Increase Nr or retain fewer target modes.'], ...
            nr, numberOfMandatory);
    end
    selected = select_independent_columns(values, numberOfMandatory, ...
        nr, options.independenceTolerance);
    selectedValues = values(:, selected);
    selectedFirst = first(:, selected);
    selectedSecond = second(:, selected);
    columnNorms = sqrt(sum(abs(selectedValues).^2, 1));
    selectedValues = selectedValues./columnNorms;
    selectedFirst = selectedFirst./columnNorms;
    selectedSecond = selectedSecond./columnNorms;
    conditionNumber = cond(full(selectedValues));
    if ~isfinite(conditionNumber) || ...
            conditionNumber > options.maximumConditionNumber
        error('vi_cylinder_radial_grid:IllConditionedBasis', ...
            ['The selected Bessel basis has condition number %.3e, ', ...
             'above the requested limit %.3e. Change Nr or the retained ', ...
             'radial branches.'], conditionNumber, ...
            options.maximumConditionNumber);
    end

    Dr = selectedFirst/selectedValues;
    Drr = selectedSecond/selectedValues;
    firstResidual = norm(Dr*selectedValues-selectedFirst, 'fro') / ...
        max(norm(selectedFirst, 'fro'), eps);
    secondResidual = norm(Drr*selectedValues-selectedSecond, 'fro') / ...
        max(norm(selectedSecond, 'fro'), eps);

    grid = struct();
    grid.r = r;
    grid.D = sparse(Dr);
    grid.D2 = sparse(Drr);
    grid.typeRequested = requestedType;
    grid.typeUsed = 'besselEnriched';
    grid.numberOfPoints = nr;
    grid.selectedBasisLabels = labels(selected);
    grid.numberOfMandatoryFunctions = numberOfMandatory;
    grid.conditionNumber = conditionNumber;
    grid.firstDerivativeBasisResidual = firstResidual;
    grid.secondDerivativeBasisResidual = secondResidual;
    grid.maximumProductOrder = options.maximumProductOrder;
    grid.fallbackReason = '';
catch basisError
    isExpectedBasisFailure = strncmp(basisError.identifier, ...
        'vi_cylinder_radial_grid:', ...
        numel('vi_cylinder_radial_grid:'));
    if ~logical(options.fallbackToChebyshev) || ...
            ~isExpectedBasisFailure
        rethrow(basisError);
    end
    warning('vi_cylinder_radial_grid:BesselFallback', ...
        ['Bessel-enriched differentiation was rejected: %s ', ...
         'Using the Chebyshev radial operator.'], basisError.message);
    grid = chebyshev_grid(r, chebyshevD, requestedType);
    grid.fallbackReason = basisError.message;
end
end

function options = default_options()
options.type = 'chebyshev';
options.maximumProductOrder = 2;
options.maximumConditionNumber = 1.0e10;
options.independenceTolerance = 1.0e-10;
options.fallbackToChebyshev = true;
end

function grid = chebyshev_grid(r, D, requestedType)
grid = struct();
grid.r = r;
grid.D = sparse(D);
grid.D2 = sparse(D*D);
grid.typeRequested = requestedType;
grid.typeUsed = 'chebyshev';
grid.numberOfPoints = numel(r);
grid.selectedBasisLabels = {};
grid.numberOfMandatoryFunctions = 0;
grid.conditionNumber = NaN;
grid.firstDerivativeBasisResidual = NaN;
grid.secondDerivativeBasisResidual = NaN;
grid.maximumProductOrder = 0;
grid.fallbackReason = '';
end

function [values, first, second, labels, numberOfMandatory] = ...
    candidate_library(parameters, r, maximumProductOrder)
if ~isfield(parameters, 'modes') || isempty(parameters.modes)
    error('vi_cylinder_radial_grid:MissingModes', ...
        ['parameters.modes is required for Bessel-enriched radial ', ...
         'differentiation.']);
end

values = ones(numel(r), 1);
first = zeros(numel(r), 1);
second = zeros(numel(r), 1);
labels = {'constant'};
atomValues = zeros(numel(r), 0);
atomFirst = zeros(numel(r), 0);
atomSecond = zeros(numel(r), 0);
atomLabels = cell(1, 0);
atomKeys = cell(1, 0);

for modeIndex = 1:numel(parameters.modes)
    mode = parameters.modes(modeIndex);
    if isfield(mode, 'betaStar') && ~isempty(mode.betaStar)
        beta = mode.betaStar;
    elseif isfield(mode, 'radialIndex')
        roots = bessel_derivative_root(mode.m, mode.radialIndex);
        beta = roots(mode.radialIndex)/parameters.R0;
    else
        error('vi_cylinder_radial_grid:MissingBeta', ...
            'Every retained mode needs betaStar or radialIndex.');
    end
    validateattributes(beta, {'numeric'}, ...
        {'scalar', 'real', 'positive', 'finite'});
    orders = unique(abs([mode.m-1, mode.m, mode.m+1]), 'stable');
    for orderIndex = 1:numel(orders)
        order = orders(orderIndex);
        key = sprintf('%d|%.14e', order, beta);
        if any(strcmp(atomKeys, key))
            continue;
        end
        [value, derivative, secondDerivative] = ...
            bessel_atom(order, beta, r);
        atomValues(:, end+1) = value; %#ok<AGROW>
        atomFirst(:, end+1) = derivative; %#ok<AGROW>
        atomSecond(:, end+1) = secondDerivative; %#ok<AGROW>
        atomLabels{end+1} = sprintf('J_%d(%.6g r)', order, beta); %#ok<AGROW>
        atomKeys{end+1} = key; %#ok<AGROW>
    end
end

values = [values, atomValues];
first = [first, atomFirst];
second = [second, atomSecond];
labels = [labels, atomLabels];
numberOfMandatory = size(values, 2);

if maximumProductOrder >= 2
    numberOfAtoms = size(atomValues, 2);
    for firstIndex = 1:numberOfAtoms
        for secondIndex = firstIndex:numberOfAtoms
            [product, productFirst, productSecond] = product_two( ...
                atomValues(:, firstIndex), atomFirst(:, firstIndex), ...
                atomSecond(:, firstIndex), ...
                atomValues(:, secondIndex), atomFirst(:, secondIndex), ...
                atomSecond(:, secondIndex));
            values(:, end+1) = product; %#ok<AGROW>
            first(:, end+1) = productFirst; %#ok<AGROW>
            second(:, end+1) = productSecond; %#ok<AGROW>
            labels{end+1} = sprintf('(%s)(%s)', ...
                atomLabels{firstIndex}, atomLabels{secondIndex}); %#ok<AGROW>
        end
    end
end

if maximumProductOrder >= 3
    numberOfAtoms = size(atomValues, 2);
    for firstIndex = 1:numberOfAtoms
        for secondIndex = firstIndex:numberOfAtoms
            [pair, pairFirst, pairSecond] = product_two( ...
                atomValues(:, firstIndex), atomFirst(:, firstIndex), ...
                atomSecond(:, firstIndex), ...
                atomValues(:, secondIndex), atomFirst(:, secondIndex), ...
                atomSecond(:, secondIndex));
            for thirdIndex = secondIndex:numberOfAtoms
                [product, productFirst, productSecond] = product_two( ...
                    pair, pairFirst, pairSecond, ...
                    atomValues(:, thirdIndex), atomFirst(:, thirdIndex), ...
                    atomSecond(:, thirdIndex));
                values(:, end+1) = product; %#ok<AGROW>
                first(:, end+1) = productFirst; %#ok<AGROW>
                second(:, end+1) = productSecond; %#ok<AGROW>
                labels{end+1} = sprintf('(%s)(%s)(%s)', ...
                    atomLabels{firstIndex}, atomLabels{secondIndex}, ...
                    atomLabels{thirdIndex}); %#ok<AGROW>
            end
        end
    end
end

% Regular powers guarantee a complete fallback space on the evaluation
% nodes. Greedy selection uses them only when they add more independent
% radial information than the available Bessel products.
x = r/parameters.R0;
for degree = 1:max(numel(r)-1, 1)
    value = x.^degree;
    derivative = (degree/parameters.R0)*x.^(degree-1);
    if degree == 1
        secondDerivative = zeros(size(x));
    else
        secondDerivative = degree*(degree-1)/parameters.R0^2 * ...
            x.^(degree-2);
    end
    values(:, end+1) = value; %#ok<AGROW>
    first(:, end+1) = derivative; %#ok<AGROW>
    second(:, end+1) = secondDerivative; %#ok<AGROW>
    labels{end+1} = sprintf('(r/R0)^%d', degree); %#ok<AGROW>
end
end

function selected = select_independent_columns( ...
    values, numberOfMandatory, numberToSelect, tolerance)
selected = 1:numberOfMandatory;
available = (numberOfMandatory+1):size(values, 2);
if numerical_rank(values(:, selected), tolerance) < numberOfMandatory
    error('vi_cylinder_radial_grid:DependentMandatoryBasis', ...
        ['The mandatory constant/Bessel functions are linearly dependent ', ...
         'on this radial grid. Change Nr.']);
end
while numel(selected) < numberToSelect
    normalizedSelected = normalize_columns(values(:, selected));
    [Q, ~] = qr(normalizedSelected, 0);
    bestScore = -Inf;
    bestAvailableIndex = NaN;
    for candidateIndex = 1:numel(available)
        candidate = values(:, available(candidateIndex));
        candidateNorm = norm(candidate);
        if candidateNorm <= eps
            continue;
        end
        candidate = candidate/candidateNorm;
        residual = candidate-Q*(Q'*candidate);
        score = norm(residual);
        if score > bestScore
            bestScore = score;
            bestAvailableIndex = candidateIndex;
        end
    end
    if ~isfinite(bestScore) || bestScore <= tolerance
        error('vi_cylinder_radial_grid:IncompleteBasis', ...
            ['Only %d independent radial functions could be selected ', ...
             'for Nr=%d.'], numel(selected), numberToSelect);
    end
    selected(end+1) = available(bestAvailableIndex); %#ok<AGROW>
    available(bestAvailableIndex) = [];
end
end

function rankValue = numerical_rank(matrix, tolerance)
singularValues = svd(normalize_columns(matrix), 'econ');
rankValue = nnz(singularValues > tolerance*max(singularValues));
end

function normalized = normalize_columns(matrix)
norms = sqrt(sum(abs(matrix).^2, 1));
norms(norms <= eps) = 1;
normalized = matrix./norms;
end

function [value, first, second] = bessel_atom(order, beta, r)
argument = beta*r;
value = besselj(order, argument);
first = 0.5*beta*(besselj(order-1, argument) - ...
    besselj(order+1, argument));
second = 0.25*beta^2*(besselj(order-2, argument) - ...
    2*besselj(order, argument) + besselj(order+2, argument));
end

function [value, first, second] = product_two( ...
    a, aFirst, aSecond, b, bFirst, bSecond)
value = a.*b;
first = aFirst.*b+a.*bFirst;
second = aSecond.*b+2*aFirst.*bFirst+a.*bSecond;
end
