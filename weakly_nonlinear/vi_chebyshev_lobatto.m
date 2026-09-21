function [x, D, quadratureWeights] = ...
        vi_chebyshev_lobatto(numberOfPoints, interval)
%VI_CHEBYSHEV_LOBATTO Grid, derivative, and Clenshaw--Curtis weights.
%
% [x,D,w] = vi_chebyshev_lobatto(N,[a,b]) returns N ascending points on
% [a,b]. D*f approximates df/dx and w'*f approximates its integral.

if nargin < 2 || isempty(interval)
    interval = [-1, 1];
end
validateattributes(numberOfPoints, {'numeric'}, ...
    {'scalar', 'integer', '>=', 2});
validateattributes(interval, {'numeric'}, ...
    {'vector', 'numel', 2, 'real', 'finite'});
assert(interval(2) > interval(1), ...
    'The right interval endpoint must exceed the left endpoint.');

n = numberOfPoints - 1;
j = (0:n).';
xReference = cos(pi*j/n);
c = [2; ones(n-1, 1); 2].*(-1).^j;
dX = xReference - xReference.';
DReference = (c*(1./c).')./(dX + eye(numberOfPoints));
DReference = DReference - diag(sum(DReference, 2));

quadratureReference = zeros(numberOfPoints,1);
if n == 1
    quadratureReference(:) = 1;
else
    interior = 2:n;
    theta = pi*j/n;
    values = ones(n-1,1);
    if mod(n,2) == 0
        quadratureReference([1,end]) = 1/(n^2-1);
        for k = 1:n/2-1
            values = values-2*cos(2*k*theta(interior))/(4*k^2-1);
        end
        values = values-cos(n*theta(interior))/(n^2-1);
    else
        quadratureReference([1,end]) = 1/n^2;
        for k = 1:(n-1)/2
            values = values-2*cos(2*k*theta(interior))/(4*k^2-1);
        end
    end
    quadratureReference(interior) = 2*values/n;
end

% Put the grid in ascending order and map it to the requested interval.
permutation = numberOfPoints:-1:1;
xReference = xReference(permutation);
DReference = DReference(permutation, permutation);
quadratureReference = quadratureReference(permutation);
scale = 2/(interval(2)-interval(1));
x = interval(1) + (interval(2)-interval(1)) * ...
    (xReference+1)/2;
D = scale*DReference;
quadratureWeights = ...
    ((interval(2)-interval(1))/2)*quadratureReference;
end
