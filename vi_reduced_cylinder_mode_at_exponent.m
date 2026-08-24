function [zeta, diagnostics] = vi_reduced_cylinder_mode_at_exponent( ...
    acceleration, omegaStar, R0, m, radialIndex, C, Bd, At, eta, N, ...
    modeType, gSign, exponent, phase)
%VI_REDUCED_CYLINDER_MODE_AT_EXPONENT Harmonic vector at a fixed exponent.
%
% Unlike the growth-rate root search, this routine does not move the
% exponent. It seeds the full WNL operator at either a neutral exponent or
% an exact off-neutral Floquet exponent, at the acceleration used by the
% full operator and with the same forcing phase.

if nargin < 14 || isempty(phase)
    phase = 0;
end

roots = bessel_derivative_root(m, radialIndex);
betaStar = roots(radialIndex)/R0;
if strcmpi(modeType, 'SH')
    numberOfHarmonics = 2*(N+1);
elseif strcmpi(modeType, 'H')
    numberOfHarmonics = 2*N+1;
else
    error('vi_reduced_cylinder_mode_at_exponent:ModeType', ...
        'modeType must be ''H'' or ''SH''.');
end

diagonal = vi_reduced_cylinder_coefficients(N, exponent, betaStar, ...
    omegaStar, At, eta, C, Bd, modeType, gSign);
forcing = vi_floquet_acceleration_matrix(numberOfHarmonics, phase);
modeMatrix = diag(diagonal)-acceleration*forcing;
[~, singularMatrix, rightVectors] = svd(modeMatrix, 'econ');
singularValues = diag(singularMatrix);
zeta = rightVectors(:, end);
zeta = zeta/max(norm(zeta), eps);

diagnostics = struct();
diagnostics.acceleration = acceleration;
diagnostics.exponent = exponent;
diagnostics.forcingPhase = phase;
diagnostics.betaStar = betaStar;
diagnostics.singularValues = singularValues;
diagnostics.relativeSingularResidual = min(singularValues)/ ...
    max(max(singularValues), eps);
diagnostics.vectorResidual = norm(modeMatrix*zeta)/ ...
    (max(norm(modeMatrix, 2), eps)*norm(zeta));
end
