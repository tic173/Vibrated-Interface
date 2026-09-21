function dA = wnl_rhs_landau( ...
        ~, A, linearCoefficient, g, phaseSensitiveG)
%WNL_RHS_LANDAU Symmetry-complete coupled cubic amplitude equations.
%
% dA_j/dt = lambda_j*A_j + sum_k [
%   g(j,k)*A_j*|A_k|^2 + h(j,k)*conj(A_j)*A_k^2].
% The third argument is the actual linear coefficient.  At an operating
% point it is the exact Floquet exponent.  A neutral-point caller may still
% pass a detuning coefficient such as (a-a_c)*mu.

A = A(:);
linearCoefficient = linearCoefficient(:);
if nargin < 5 || isempty(phaseSensitiveG)
    phaseSensitiveG = complex(zeros(size(g)));
end
if ~isequal(size(g), [numel(A), numel(A)])
    error('wnl_rhs_landau:BadG', ...
        'g must be a square matrix with one row per amplitude.');
end
if ~isequal(size(phaseSensitiveG),size(g))
    error('wnl_rhs_landau:BadPhaseSensitiveG', ...
        'h must have the same size as g.');
end
dA = linearCoefficient.*A + A.*(g*abs(A).^2) + ...
    conj(A).*(phaseSensitiveG*(A.^2));
end
