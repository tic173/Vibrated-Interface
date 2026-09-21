function [vector,diagnostics] = vi_axisymmetric_separable_direct( ...
        spec,zetaTemporalCoefficients,metadata,parameters)
%VI_AXISYMMETRIC_SEPARABLE_DIRECT Lift an exact m=0 Bessel Floquet mode.
%
% For a free contact line, the true and legacy stress-free sidewall models
% are identical on an axisymmetric mode: u_theta=0, while zeta, w, and p are
% proportional to J_0(beta*r), and u_r is proportional to J_1(beta*r).
% This routine evaluates the same stable finite-depth vertical Stokes basis
% used by vi_reduced_cylinder_coefficients and returns the corresponding full
% primitive-variable Floquet vector.

assert(isstruct(spec) && spec.m == 0, ...
    'vi_axisymmetric_separable_direct requires spec.m=0.');
assert(isfield(spec,'betaStar') && ~isempty(spec.betaStar), ...
    'spec.betaStar is required.');
assert(isfield(metadata,'contactLine') && ...
    strcmpi(metadata.contactLine,'free'), ...
    'The analytic separable lift requires a free contact line.');
requiredParameters = {'omegaStar','C','At','eta'};
for fieldIndex = 1:numel(requiredParameters)
    assert(isfield(parameters,requiredParameters{fieldIndex}), ...
        'parameters.%s is required.',requiredParameters{fieldIndex});
end

zetaTemporalCoefficients = zetaTemporalCoefficients(:).';
if numel(zetaTemporalCoefficients) ~= numel(spec.n)
    error('vi_axisymmetric_separable_direct:TemporalSize', ...
        'One interface coefficient is required for every Floquet harmonic.');
end
layout = metadata.layout;
r = metadata.discretization.r(:);
zDense = metadata.discretization.zD(:).';
zLight = metadata.discretization.zL(:).';
beta = spec.betaStar;
bessel0 = besselj(0,beta*r);
bessel1 = besselj(1,beta*r);

coefficients = complex(zeros(metadata.ndof,numel(spec.n)));
verticalConditionNumbers = zeros(numel(spec.n),1);
for harmonicIndex = 1:numel(spec.n)
    temporalRate = wnl_spec_lambda(spec) + 1i*parameters.omegaStar* ...
        (spec.n(harmonicIndex)+spec.s);
    [dense,light,verticalConditionNumbers(harmonicIndex)] = ...
        vertical_profiles(temporalRate,beta,parameters,zDense,zLight);
    amplitude = zetaTemporalCoefficients(harmonicIndex);
    coefficients(layout.d.ur,harmonicIndex) = ...
        reshape(amplitude*bessel1*dense.ur,[],1);
    coefficients(layout.d.w,harmonicIndex) = ...
        reshape(amplitude*bessel0*dense.w,[],1);
    coefficients(layout.d.p,harmonicIndex) = ...
        reshape(amplitude*bessel0*dense.p,[],1);
    coefficients(layout.l.ur,harmonicIndex) = ...
        reshape(amplitude*bessel1*light.ur,[],1);
    coefficients(layout.l.w,harmonicIndex) = ...
        reshape(amplitude*bessel0*light.w,[],1);
    coefficients(layout.l.p,harmonicIndex) = ...
        reshape(amplitude*bessel0*light.p,[],1);
    coefficients(layout.zeta,harmonicIndex) = amplitude*bessel0;

    % The primitive operator removes one common pressure gauge in every m=0
    % temporal block. Subtracting the same constant from both fluids leaves
    % momentum and the interfacial pressure jump unchanged.
    gaugeIr = metadata.discretization.denseGaugeRadialIndex;
    gaugeIz = metadata.discretization.denseGaugeVerticalIndex;
    densePressure = reshape(coefficients( ...
        layout.d.p,harmonicIndex),layout.nr,layout.nzD);
    lightPressure = reshape(coefficients( ...
        layout.l.p,harmonicIndex),layout.nr,layout.nzL);
    gaugeValue = densePressure(gaugeIr,gaugeIz);
    densePressure = densePressure-gaugeValue;
    lightPressure = lightPressure-gaugeValue;
    coefficients(layout.d.p,harmonicIndex) = densePressure(:);
    coefficients(layout.l.p,harmonicIndex) = lightPressure(:);
end

vector = coefficients(:);
diagnostics = struct();
diagnostics.method = ...
    'exact axisymmetric Bessel / finite-depth vertical Stokes lift';
diagnostics.verticalConditionNumbers = verticalConditionNumbers;
diagnostics.maximumVerticalConditionNumber = ...
    max(verticalConditionNumbers);
diagnostics.interfaceCorrectionRelativeL2 = 0;
diagnostics.exactlySeparable = true;
end

function [dense,light,conditionNumber] = vertical_profiles( ...
        lambda,beta,parameters,zDense,zLight)
densityRatio = (1-parameters.At)/(1+parameters.At);
qDense = sqrt(1+lambda/(parameters.C*beta^2));
qLight = sqrt(1+lambda*densityRatio / ...
    (parameters.eta*parameters.C*beta^2));
denseRates = beta*[1,-1,qDense,-qDense];
lightRates = beta*[1,-1,-qLight,qLight];
denseAtInterface = [1,1,1,exp(-qDense*beta)];
denseAtBottom = [exp(-beta),exp(beta),exp(-qDense*beta),1];
lightAtInterface = [1,1,1,exp(-qLight*beta)];
lightAtTop = [exp(beta),exp(-beta),exp(-qLight*beta),1];
denseFirstAtInterface = denseRates.*denseAtInterface;
lightFirstAtInterface = lightRates.*lightAtInterface;
denseSecondAtInterface = denseRates.^2.*denseAtInterface;
lightSecondAtInterface = lightRates.^2.*lightAtInterface;

matrix = complex(zeros(8,8));
matrix(1,1:4) = denseAtInterface;
matrix(2,1:4) = denseAtInterface;
matrix(2,5:8) = -lightAtInterface;
matrix(3,1:4) = denseFirstAtInterface;
matrix(3,5:8) = -lightFirstAtInterface;
matrix(4,1:4) = denseSecondAtInterface + beta^2*denseAtInterface;
matrix(4,5:8) = -parameters.eta*( ...
    lightSecondAtInterface+beta^2*lightAtInterface);
matrix(5,1:4) = denseAtBottom;
matrix(6,1:4) = denseRates.*denseAtBottom;
matrix(7,5:8) = lightAtTop;
matrix(8,5:8) = lightRates.*lightAtTop;
rhs = [lambda;0;0;0;0;0;0;0];
[scaledMatrix,scaledRhs,rowScale,columnScale] = ...
    equilibrate_small_system(matrix,rhs);
scaledCoefficients = scaledMatrix\scaledRhs;
verticalCoefficients = columnScale.*scaledCoefficients;
conditionNumber = cond(scaledMatrix);

denseBasis = [exp(beta*zDense);exp(-beta*zDense); ...
    exp(qDense*beta*zDense);exp(-qDense*beta*(zDense+1))];
lightBasis = [exp(beta*zLight);exp(-beta*zLight); ...
    exp(-qLight*beta*zLight);exp(qLight*beta*(zLight-1))];
denseCoefficients = verticalCoefficients(1:4).';
lightCoefficients = verticalCoefficients(5:8).';
denseW = denseCoefficients*denseBasis;
lightW = lightCoefficients*lightBasis;
denseFirst = (denseCoefficients.*denseRates)*denseBasis;
lightFirst = (lightCoefficients.*lightRates)*lightBasis;
denseThird = (denseCoefficients.*(denseRates.^3))*denseBasis;
lightThird = (lightCoefficients.*(lightRates.^3))*lightBasis;

dense.w = denseW;
dense.ur = -denseFirst/beta;
dense.p = parameters.C*denseThird/beta^2 - ...
    (lambda/beta^2+parameters.C)*denseFirst;
light.w = lightW;
light.ur = -lightFirst/beta;
light.p = parameters.eta*parameters.C*lightThird/beta^2 - ...
    (densityRatio*lambda/beta^2+ ...
    parameters.eta*parameters.C)*lightFirst;

% Keep the scales in diagnostics-friendly scope and silence an otherwise
% unused-output warning in older MATLAB releases.
assert(all(isfinite(rowScale)));
end

function [scaledMatrix,scaledRhs,rowScale,columnScale] = ...
        equilibrate_small_system(matrix,rhs)
rowNorm = sqrt(sum(abs(matrix).^2,2));
rowScale = 1./max(rowNorm,eps);
scaledMatrix = rowScale.*matrix;
scaledRhs = rowScale.*rhs;
columnNorm = sqrt(sum(abs(scaledMatrix).^2,1)).';
columnScale = 1./max(columnNorm,eps);
scaledMatrix = scaledMatrix.*columnScale.';
end
