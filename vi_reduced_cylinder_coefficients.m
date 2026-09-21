function An = vi_reduced_cylinder_coefficients(n, sStar, betaStar, ...
    omegaStar, At, eta, C, Bd, modeType, gSign, layerDepths)
%VI_REDUCED_CYLINDER_COEFFICIENTS Finite-depth separable Schur diagonal.
% Optional layerDepths = [lowerDepth upperDepth]/referenceDepth, default [1 1].
% betaStar is k*h for Cartesian Fourier or cylindrical Bessel modes.
%
% This is the single implementation used by both the neutral-acceleration
% and Floquet-growth solvers.  For every temporal harmonic it solves the
% vertical Stokes problem with no slip at both layer walls, velocity and
% tangential-stress continuity at z=0, and unit interface displacement.
%
% The original code switched to a semi-infinite lower layer when a large
% Stokes wavenumber made exp(q*beta) overflow.  That switch removed the
% bottom-wall conditions and made its mode incompatible with the full WNL
% operator.  The basis used here contains only decaying exponentials such
% as exp(-q*beta); it enforces both finite walls without overflow.

if nargin < 11 || isempty(layerDepths), layerDepths = [1 1]; end
validateattributes(layerDepths, {'numeric'}, ...
    {'real','finite','positive','vector','numel',2});

if strcmpi(modeType, 'SH')
    harmonicIndices = -n-1:n;
elseif strcmpi(modeType, 'H')
    harmonicIndices = -n:n;
else
    error('vi_reduced_cylinder_coefficients:ModeType', ...
        'modeType must be ''H'' or ''SH''.');
end

densityRatio = (1-At)/(1+At);
An = complex(zeros(numel(harmonicIndices), 1));

for harmonicPosition = 1:numel(harmonicIndices)
    temporalIndex = harmonicIndices(harmonicPosition);
    lambda = sStar+1i*temporalIndex*omegaStar;

    % The exactly static harmonic has no velocity field.  Treat it
    % analytically because qDense=qLight=1 makes any four-exponential
    % representation linearly dependent.
    if lambda == 0
        An(harmonicPosition) = 2*(gSign + ...
            betaStar^2*(1+At)/(2*Bd*At));
        continue;
    end

    qDense = sqrt(1+lambda/(C*betaStar^2));
    qLight = sqrt(1+lambda*densityRatio/(eta*C*betaStar^2));
    [denseFirstDerivative, denseThirdOverBetaSquared, ...
        lightThirdOverBetaSquared] = interface_derivatives( ...
        lambda, betaStar, qDense, qLight, eta, layerDepths);

    % From vertical and horizontal momentum,
    %   p = mu*w'''/beta^2-(rho*lambda/beta^2+mu)*w'.
    % Substitution into the normal-traction jump, divided by half the
    % density jump, gives the diagonal below.  Both fluids' w''' terms are
    % retained separately.
    An(harmonicPosition) = ...
        2*denseFirstDerivative*(lambda/betaStar^2 + ...
        3*C*(1-eta)*(1+At)/(2*At)) ...
        - C*(1+At)/At*(denseThirdOverBetaSquared - ...
        eta*lightThirdOverBetaSquared) ...
        + 2*(gSign + betaStar^2*(1+At)/(2*Bd*At));
end
end

function [denseFirst, denseThirdOverBetaSquared, ...
    lightThirdOverBetaSquared] = interface_derivatives( ...
    lambda, beta, qDense, qLight, eta, layerDepths)
% Shift every growing exponential to its adjacent wall, so basis values
% stay bounded at BOTH ends even for large beta or unequal depths.
% "dense" and "light" below are legacy names for lower and upper fluid.
hLower = layerDepths(1); hUpper = layerDepths(2);
denseRates = beta*[1, -1, qDense, -qDense];
lightRates = beta*[1, -1, -qLight, qLight];
denseAtInterface = [1, exp(-beta*hLower), 1, exp(-qDense*beta*hLower)];
denseAtBottom = [exp(-beta*hLower), 1, exp(-qDense*beta*hLower), 1];
lightAtInterface = [exp(-beta*hUpper), 1, 1, exp(-qLight*beta*hUpper)];
lightAtTop = [1, exp(-beta*hUpper), exp(-qLight*beta*hUpper), 1];

denseDerivativeAtInterface = denseRates.*denseAtInterface;
lightDerivativeAtInterface = lightRates.*lightAtInterface;
denseSecondAtInterface = denseRates.^2.*denseAtInterface;
lightSecondAtInterface = lightRates.^2.*lightAtInterface;

leftHandSide = complex(zeros(8, 8));
% Kinematic condition and velocity continuity.
leftHandSide(1, 1:4) = denseAtInterface;
leftHandSide(2, 1:4) = denseAtInterface;
leftHandSide(2, 5:8) = -lightAtInterface;
leftHandSide(3, 1:4) = denseDerivativeAtInterface;
leftHandSide(3, 5:8) = -lightDerivativeAtInterface;
% Tangential traction: (wD''+beta^2*wD)
%                    -eta*(wL''+beta^2*wL)=0.
leftHandSide(4, 1:4) = denseSecondAtInterface + ...
    beta^2*denseAtInterface;
leftHandSide(4, 5:8) = -eta*(lightSecondAtInterface + ...
    beta^2*lightAtInterface);
% No penetration and no slip at the two horizontal walls.
leftHandSide(5, 1:4) = denseAtBottom;
leftHandSide(6, 1:4) = denseRates.*denseAtBottom;
leftHandSide(7, 5:8) = lightAtTop;
leftHandSide(8, 5:8) = lightRates.*lightAtTop;
rightHandSide = [lambda;0;0;0;0;0;0;0];

% Equilibrate rows and columns before the small dense solve.  This is
% important when |q*beta| is hundreds or thousands, even though the basis
% values themselves no longer overflow.
rowNorm = sqrt(sum(abs(leftHandSide).^2, 2));
rowScale = 1./max(rowNorm, eps);
scaledMatrix = rowScale.*leftHandSide;
scaledRightHandSide = rowScale.*rightHandSide;
columnNorm = sqrt(sum(abs(scaledMatrix).^2, 1)).';
columnScale = 1./max(columnNorm, eps);
scaledMatrix = scaledMatrix.*columnScale.';
scaledCoefficients = scaledMatrix\scaledRightHandSide;
coefficients = columnScale.*scaledCoefficients;

denseCoefficients = coefficients(1:4).';
lightCoefficients = coefficients(5:8).';
denseFirst = sum(denseCoefficients.*denseDerivativeAtInterface);
denseThirdOverBetaSquared = sum(denseCoefficients.* ...
    (denseRates.^3.*denseAtInterface))/beta^2;
lightThirdOverBetaSquared = sum(lightCoefficients.* ...
    (lightRates.^3.*lightAtInterface))/beta^2;
end
