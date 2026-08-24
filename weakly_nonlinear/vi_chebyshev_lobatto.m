function [x, D] = vi_chebyshev_lobatto(numberOfPoints, interval)
%VI_CHEBYSHEV_LOBATTO Chebyshev--Lobatto grid and first derivative.
%
% [x,D] = vi_chebyshev_lobatto(N,[a,b]) returns N ascending points on
% [a,b]. For values f sampled at x, D*f approximates df/dx.

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

% Put the grid in ascending order and map it to the requested interval.
permutation = numberOfPoints:-1:1;
xReference = xReference(permutation);
DReference = DReference(permutation, permutation);
scale = 2/(interval(2)-interval(1));
x = interval(1) + (interval(2)-interval(1)) * ...
    (xReference+1)/2;
D = scale*DReference;
end

