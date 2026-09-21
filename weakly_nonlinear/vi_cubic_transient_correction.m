function result = vi_cubic_transient_correction(timeRequested, A0, lambda, ...
        g, limits, phaseSensitiveG)
%VI_CUBIC_TRANSIENT_CORRECTION First nonlinear correction to linear growth.
%
% A_L,j(t) = A_j(0)*exp(lambda_j*t)
% delta A_j^(3)(t) = A_L,j(t)*sum_k g(j,k)*|A_k(0)|^2*I_k(t),
% I_k(t) = integral_0^t exp(2*Re(lambda_k)*s) ds.
%
% The cubic forcing is evaluated on A_L, so this is an O(A^3) finite-time
% perturbation rather than a nonlinear saturation model. The returned record
% is truncated before amplitude, relative-correction, or rate-ratio limits
% are exceeded.

if nargin < 5 || isempty(limits)
    limits = struct();
end
if nargin < 6 || isempty(phaseSensitiveG)
    phaseSensitiveG = complex(zeros(size(g)));
end
timeRequested = timeRequested(:);
if numel(timeRequested) < 2 || any(~isfinite(timeRequested)) || ...
        any(diff(timeRequested) <= 0)
    error('vi_cubic_transient_correction:TimeVector', ...
        'timeRequested must contain at least two increasing finite values.');
end
A0 = A0(:);
lambda = lambda(:);
numberOfModes = numel(A0);
if numel(lambda) ~= numberOfModes || ...
        ~isequal(size(g), [numberOfModes, numberOfModes]) || ...
        ~isequal(size(phaseSensitiveG),size(g))
    error('vi_cubic_transient_correction:DimensionMismatch', ...
        ['A0, lambda, g, and phaseSensitiveG must describe the same ', ...
         'number of modes.']);
end
if any(~isfinite(real(A0))) || any(~isfinite(imag(A0))) || ...
        any(~isfinite(real(lambda))) || any(~isfinite(imag(lambda))) || ...
        any(~isfinite(real(g(:)))) || any(~isfinite(imag(g(:)))) || ...
        any(~isfinite(real(phaseSensitiveG(:)))) || ...
        any(~isfinite(imag(phaseSensitiveG(:))))
    error('vi_cubic_transient_correction:NonfiniteInput', ...
        'A0, lambda, g, and phaseSensitiveG must be finite.');
end

maximumAmplitude = get_limit(limits, 'maximumAmplitude', Inf, true);
maximumRelativeCorrection = get_limit( ...
    limits, 'maximumRelativeCorrection', Inf, true);
maximumRateRatio = get_limit( ...
    limits, 'maximumNonlinearToLinearRateRatio', Inf, true);
linearRateFloor = get_limit(limits, 'linearRateFloor', 1.0e-12, false);

elapsedTime = timeRequested-timeRequested(1);
linearAmplitudes = exp(elapsedTime*lambda.') .* A0.';
intensityIntegrals = zeros(numel(elapsedTime), numberOfModes);
for modeIndex = 1:numberOfModes
    twiceGrowthRate = 2*real(lambda(modeIndex));
    if abs(twiceGrowthRate)*max(elapsedTime) < 1.0e-6
        intensityIntegrals(:,modeIndex) = elapsedTime .* (1 + ...
            0.5*twiceGrowthRate*elapsedTime + ...
            (twiceGrowthRate*elapsedTime).^2/6);
    else
        intensityIntegrals(:,modeIndex) = ...
            expm1(twiceGrowthRate*elapsedTime)/twiceGrowthRate;
    end
end
weightedIntegrals = intensityIntegrals .* abs(A0.').^2;
cubicMultipliers = weightedIntegrals*g.';
crossSaturationCorrections = linearAmplitudes.*cubicMultipliers;
phaseSensitiveCorrections = complex(zeros(size(linearAmplitudes)));
for targetMode = 1:numberOfModes
    for sourceMode = 1:numberOfModes
        coefficient = phaseSensitiveG(targetMode,sourceMode)* ...
            conj(A0(targetMode))*A0(sourceMode)^2;
        if coefficient == 0
            continue;
        end
        exponent = conj(lambda(targetMode))+2*lambda(sourceMode)- ...
            lambda(targetMode);
        integral = exponential_integral(exponent,elapsedTime);
        phaseSensitiveCorrections(:,targetMode) = ...
            phaseSensitiveCorrections(:,targetMode) + ...
            exp(lambda(targetMode)*elapsedTime).*coefficient.*integral;
    end
end
cubicCorrections = ...
    crossSaturationCorrections+phaseSensitiveCorrections;
correctedAmplitudes = linearAmplitudes+cubicCorrections;

linearScale = max(abs(linearAmplitudes), 1.0e-14);
relativeCorrections = abs(cubicCorrections)./linearScale;
crossSaturationRhs = linearAmplitudes.* ...
    (abs(linearAmplitudes).^2*g.');
phaseSensitiveRhs = conj(linearAmplitudes).* ...
    ((linearAmplitudes.^2)*phaseSensitiveG.');
nonlinearRhs = crossSaturationRhs+phaseSensitiveRhs;
nonlinearRateCorrections = nonlinearRhs./ ...
    max(abs(linearAmplitudes),1.0e-14);
rateDenominator = max(abs(linearAmplitudes.*lambda.'), ...
    linearRateFloor*max(abs(linearAmplitudes),1.0e-14));
nonlinearToLinearRateRatios = abs(nonlinearRhs)./rateDenominator;

amplitudeMetric = max(abs(correctedAmplitudes), [], 2);
relativeCorrectionMetric = max(relativeCorrections, [], 2);
rateRatioMetric = max(nonlinearToLinearRateRatios, [], 2);
pointValid = amplitudeMetric <= maximumAmplitude & ...
    relativeCorrectionMetric <= maximumRelativeCorrection & ...
    rateRatioMetric <= maximumRateRatio;
firstRejectedIndex = find(~pointValid, 1, 'first');
if isempty(firstRejectedIndex)
    keepCount = numel(timeRequested);
    firstRejectedTime = NaN;
    stopReasons = {};
else
    keepCount = max(1, firstRejectedIndex-1);
    firstRejectedTime = timeRequested(firstRejectedIndex);
    stopReasons = {};
    if amplitudeMetric(firstRejectedIndex) > maximumAmplitude
        stopReasons{end+1} = 'amplitude limit';
    end
    if relativeCorrectionMetric(firstRejectedIndex) > ...
            maximumRelativeCorrection
        stopReasons{end+1} = 'relative cubic-correction limit';
    end
    if rateRatioMetric(firstRejectedIndex) > maximumRateRatio
        stopReasons{end+1} = ...
            'nonlinear-to-linear rate-ratio limit';
    end
end
keep = 1:keepCount;

result = struct();
result.model = 'first operating-point cubic transient correction';
result.time = timeRequested(keep);
result.linearAmplitudes = linearAmplitudes(keep,:);
result.cubicCorrections = cubicCorrections(keep,:);
result.crossSaturationCorrections = ...
    crossSaturationCorrections(keep,:);
result.phaseSensitiveCorrections = ...
    phaseSensitiveCorrections(keep,:);
result.correctedAmplitudes = correctedAmplitudes(keep,:);
result.relativeCorrections = relativeCorrections(keep,:);
result.nonlinearRateCorrections = nonlinearRateCorrections(keep,:);
result.nonlinearRhsTerms = nonlinearRhs(keep,:);
result.phaseSensitiveG = phaseSensitiveG;
result.nonlinearToLinearRateRatios = ...
    nonlinearToLinearRateRatios(keep,:);
result.initialConditionValid = pointValid(1);
result.returnedWindowValid = all(pointValid(keep));
result.completedRequestedWindow = isempty(firstRejectedIndex);
result.stopTime = timeRequested(keepCount);
result.firstRejectedTime = firstRejectedTime;
result.stopReasons = stopReasons;
result.maximumAmplitude = maximumAmplitude;
result.maximumRelativeCorrection = maximumRelativeCorrection;
result.maximumNonlinearToLinearRateRatio = maximumRateRatio;
result.linearRateFloor = linearRateFloor;
result.maximumReachedAmplitude = max(amplitudeMetric(keep));
result.maximumReachedRelativeCorrection = ...
    max(relativeCorrectionMetric(keep));
result.maximumReachedNonlinearToLinearRateRatio = ...
    max(rateRatioMetric(keep));
end

function integral = exponential_integral(exponent,time)
if abs(exponent)*max(time) < 1.0e-6
    integral = time.*(1+0.5*exponent*time+ ...
        (exponent*time).^2/6);
else
    integral = expm1(exponent*time)/exponent;
end
end

function value = get_limit(limits, name, defaultValue, requirePositive)
if isfield(limits, name) && ~isempty(limits.(name))
    value = limits.(name);
else
    value = defaultValue;
end
attributes = {'scalar','real','nonnegative'};
if requirePositive
    attributes = {'scalar','real','positive'};
end
validateattributes(value, {'numeric'}, attributes);
end
