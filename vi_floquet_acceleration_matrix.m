function B = vi_floquet_acceleration_matrix(numberOfHarmonics, phase)
%VI_FLOQUET_ACCELERATION_MATRIX Reduced cosine-forcing coupling matrix.
%
% The primitive-cylinder residual contains
%
%   + a*cos(tau+phase)*zeta.
%
% After the normal-stress equation is divided by half the density jump,
% its row-n contribution is
%
%   + a*exp(+i*phase)*zeta_(n-1)
%   + a*exp(-i*phase)*zeta_(n+1).
%
% The legacy reduced solvers write their matrix as diag(An)-a*B.  B is
% therefore the negative of the two off-diagonal coefficients above.  This
% convention makes the reduced harmonic vector use exactly the same time
% origin and forcing phase as vi_cylinder_wnl_linear_operators.

validateattributes(numberOfHarmonics, {'numeric'}, ...
    {'scalar', 'integer', 'positive', 'finite'});
validateattributes(phase, {'numeric'}, ...
    {'scalar', 'real', 'finite'});

lower = exp(1i*phase)*ones(numberOfHarmonics, 1);
upper = exp(-1i*phase)*ones(numberOfHarmonics, 1);
cosineCoupling = spdiags([lower, zeros(numberOfHarmonics, 1), upper], ...
    -1:1, numberOfHarmonics, numberOfHarmonics);
B = -full(cosineCoupling);
end
