function grid = vi_multidomain_chebyshev_grid(breaks, pointsPerElement)
%VI_MULTIDOMAIN_CHEBYSHEV_GRID Piecewise Chebyshev--Lobatto grid.
%
% Adjacent elements retain separate copies of their common endpoint.  The
% caller must impose matching conditions between the returned interface
% index pairs.  Differentiation matrices are block diagonal and therefore
% differentiate independently inside every spectral element.

breaks = breaks(:).';
pointsPerElement = pointsPerElement(:).';
if numel(breaks) < 2 || any(~isfinite(breaks)) || ...
        any(diff(breaks) <= 0)
    error('vi_multidomain_chebyshev_grid:BadBreaks', ...
        'breaks must be a finite strictly increasing vector.');
end
numberOfElements = numel(breaks)-1;
if isscalar(pointsPerElement)
    pointsPerElement = repmat(pointsPerElement, 1, numberOfElements);
end
if numel(pointsPerElement) ~= numberOfElements || ...
        any(pointsPerElement < 3) || ...
        any(pointsPerElement ~= round(pointsPerElement))
    error('vi_multidomain_chebyshev_grid:BadPointCounts', ...
        ['pointsPerElement must contain one integer >=3 for every ', ...
         'spectral element.']);
end

z = zeros(sum(pointsPerElement), 1);
elementIndices = cell(numberOfElements, 1);
derivativeBlocks = cell(numberOfElements, 1);
secondDerivativeBlocks = cell(numberOfElements, 1);
next = 0;
for elementIndex = 1:numberOfElements
    numberOfPoints = pointsPerElement(elementIndex);
    [elementZ, elementD] = vi_chebyshev_lobatto( ...
        numberOfPoints, breaks(elementIndex:elementIndex+1));
    indices = next+(1:numberOfPoints);
    z(indices) = elementZ;
    elementIndices{elementIndex} = indices;
    derivativeBlocks{elementIndex} = sparse(elementD);
    secondDerivativeBlocks{elementIndex} = sparse(elementD*elementD);
    next = next+numberOfPoints;
end

if numberOfElements > 1
    leftIndices = cumsum(pointsPerElement(1:end-1));
    rightIndices = leftIndices+1;
    interfacePairs = [leftIndices(:), rightIndices(:)];
    interfaceCoordinates = breaks(2:end-1).';
else
    interfacePairs = zeros(0, 2);
    interfaceCoordinates = zeros(0, 1);
end

grid = struct();
grid.z = z;
grid.D = blkdiag(derivativeBlocks{:});
grid.D2 = blkdiag(secondDerivativeBlocks{:});
grid.breaks = breaks;
grid.pointsPerElement = pointsPerElement;
grid.numberOfElements = numberOfElements;
grid.elementIndices = elementIndices;
grid.interfacePairs = interfacePairs;
grid.interfaceCoordinates = interfaceCoordinates;
grid.numberOfPoints = numel(z);
end
