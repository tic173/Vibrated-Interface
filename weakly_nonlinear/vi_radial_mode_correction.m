function analysis = vi_radial_mode_correction( ...
        mode,metadata,parameters,userSettings)
%VI_RADIAL_MODE_CORRECTION Separate a Floquet mode from its Bessel seed.
%
% For every temporal harmonic,
%
%   zeta_n(r) = b_n*J_m(beta*r) + c_n(r),
%
% where c_n is orthogonal to J_m(beta*r) in the cylindrical surface inner
% product integral_0^R conj(f)*g*r*dr. The returned selection flag controls
% surface reconstruction only. The full recovered state remains the state
% used by the linear operator and nonlinear solvability calculations.

if nargin < 4
    userSettings = struct();
end
assert(isstruct(mode) && isfield(mode,'field') && ...
    isstruct(mode.field) && isfield(mode.field,'coeff') && ...
    isfield(mode.field,'spec'), ...
    'mode.field with coeff and spec is required.');
assert(isfield(metadata,'layout') && isfield(metadata.layout,'zeta'), ...
    'metadata.layout.zeta is required.');
assert(isfield(metadata,'discretization') && ...
    isfield(metadata.discretization,'r'), ...
    'metadata.discretization.r is required.');

field = mode.field;
spec = field.spec;
if ~isfield(spec,'betaStar') || isempty(spec.betaStar)
    if isfield(mode,'spec') && isfield(mode.spec,'betaStar') && ...
            ~isempty(mode.spec.betaStar)
        betaStar = mode.spec.betaStar;
    else
        error('vi_radial_mode_correction:MissingBeta', ...
            'The recovered mode specification must contain betaStar.');
    end
else
    betaStar = spec.betaStar;
end
validateattributes(betaStar,{'numeric'}, ...
    {'scalar','real','positive','finite'});

r = metadata.discretization.r(:);
zetaRows = metadata.layout.zeta(:);
if size(field.coeff,1) < max(zetaRows) || ...
        size(field.coeff,2) ~= numel(spec.n)
    error('vi_radial_mode_correction:CoefficientSize', ...
        'The mode field is incompatible with the interface layout/spec.');
end
zetaCoefficients = field.coeff(zetaRows,:);
weights = cylindrical_weights(metadata,r);

besselShape = besselj(abs(spec.m),betaStar*r);
besselPeak = max(abs(besselShape));
if besselPeak <= eps*max(norm(besselShape),1)
    error('vi_radial_mode_correction:ZeroBesselShape', ...
        'The nominal Bessel shape is zero on the radial grid.');
end
besselShape = besselShape/besselPeak;
denominator = real(besselShape'*(weights.*besselShape));
if denominator <= eps*max(norm(weights),1)
    error('vi_radial_mode_correction:DegenerateInnerProduct', ...
        'The cylindrical Bessel inner product is numerically zero.');
end
besselTemporalCoefficients = ...
    (besselShape'*(weights.*zetaCoefficients))/denominator;
separableCoefficients = ...
    besselShape*besselTemporalCoefficients;
correctionCoefficients = zetaCoefficients-separableCoefficients;

fullNorm = weighted_frobenius_norm(zetaCoefficients,weights);
besselNorm = weighted_frobenius_norm(separableCoefficients,weights);
correctionNorm = weighted_frobenius_norm( ...
    correctionCoefficients,weights);
relativeCorrection = correctionNorm/max(fullNorm,eps);

harmonicFullNorms = weighted_column_norms(zetaCoefficients,weights);
harmonicCorrectionNorms = weighted_column_norms( ...
    correctionCoefficients,weights);
perHarmonicRelativeCorrection = zeros(size(harmonicFullNorms));
activeHarmonics = harmonicFullNorms > ...
    64*eps*max(max(harmonicFullNorms),1);
perHarmonicRelativeCorrection(activeHarmonics) = ...
    harmonicCorrectionNorms(activeHarmonics)./ ...
    harmonicFullNorms(activeHarmonics);

policy = validatestring(setting(userSettings, ...
    'radialCorrectionPolicy','auto'),{'auto','always','never'});
threshold = setting(userSettings, ...
    'radialCorrectionRelativeL2Threshold',1.0e-2);
validateattributes(threshold,{'numeric'}, ...
    {'scalar','real','nonnegative','finite'});
numberOfTemporalSamples = setting(userSettings, ...
    'radialCorrectionTemporalSamples',128);
validateattributes(numberOfTemporalSamples,{'numeric'}, ...
    {'scalar','integer','>=',16,'finite'});
exactlySeparableBoundaryModel = ...
    isfield(metadata,'contactLine') && ...
    strcmpi(metadata.contactLine,'free') && ...
    isfield(metadata,'sidewallTangentialCondition') && ...
    (spec.m == 0 || strcmp(metadata.sidewallTangentialCondition, ...
    'legacyFreeSlide'));

switch policy
    case 'always'
        includeCorrection = true;
    case 'never'
        includeCorrection = false;
    otherwise
        includeCorrection = ~exactlySeparableBoundaryModel && ...
            relativeCorrection >= threshold;
end
if includeCorrection
    selectedCoefficients = zetaCoefficients;
else
    selectedCoefficients = separableCoefficients;
end

tau = (0:numberOfTemporalSamples-1)*(2*pi/numberOfTemporalSamples);
frequency = spec.n(:)+spec.s;
phase = exp(1i*frequency*tau);
fullCarrier = zetaCoefficients*phase;
besselCarrier = separableCoefficients*phase;
correctionCarrier = correctionCoefficients*phase;
fullTimeNorm = weighted_column_norms(fullCarrier,weights);
correctionTimeNorm = weighted_column_norms(correctionCarrier,weights);
timeResolvedRelativeCorrection = zeros(size(fullTimeNorm));
activeTimes = fullTimeNorm > 64*eps*max(max(fullTimeNorm),1);
timeResolvedRelativeCorrection(activeTimes) = ...
    correctionTimeNorm(activeTimes)./fullTimeNorm(activeTimes);

[~,dominantHarmonicPosition] = max(harmonicFullNorms);
orthogonalityNumerator = abs( ...
    besselShape'*(weights.*correctionCoefficients));
orthogonalityScale = sqrt(denominator)*max( ...
    harmonicCorrectionNorms,eps);

analysis = struct();
analysis.definition = ...
    'zeta_n(r)=b_n*J_m(beta*r)+c_n(r), <J_m,c_n>_{r dr}=0';
analysis.label = spec.label;
analysis.m = spec.m;
analysis.betaStar = betaStar;
analysis.r = r;
if isfield(parameters,'R0') && ~isempty(parameters.R0)
    analysis.rOverR = r/parameters.R0;
else
    analysis.rOverR = r/max(r);
end
analysis.quadratureWeightsRdr = weights;
analysis.harmonicIndices = spec.n(:).';
analysis.frequencies = frequency(:).';
analysis.besselShape = besselShape;
analysis.besselTemporalCoefficients = besselTemporalCoefficients;
analysis.fullCoefficients = zetaCoefficients;
analysis.separableCoefficients = separableCoefficients;
analysis.correctionCoefficients = correctionCoefficients;
analysis.selectedCoefficients = selectedCoefficients;
analysis.fullWeightedL2 = fullNorm;
analysis.besselWeightedL2 = besselNorm;
analysis.correctionWeightedL2 = correctionNorm;
analysis.relativeCorrectionL2 = relativeCorrection;
analysis.correctionEnergyFraction = relativeCorrection^2;
analysis.harmonicFullWeightedL2 = harmonicFullNorms;
analysis.harmonicCorrectionWeightedL2 = harmonicCorrectionNorms;
analysis.perHarmonicRelativeCorrectionL2 = ...
    perHarmonicRelativeCorrection;
analysis.dominantHarmonicPosition = dominantHarmonicPosition;
analysis.dominantHarmonicIndex = spec.n(dominantHarmonicPosition);
analysis.maximumTimeResolvedRelativeCorrectionL2 = ...
    max(timeResolvedRelativeCorrection);
analysis.peakCorrectionFraction = max(abs(correctionCarrier(:))) / ...
    max(max(abs(fullCarrier(:))),eps);
analysis.orthogonalityResidual = max( ...
    orthogonalityNumerator./orthogonalityScale);
analysis.policy = policy;
analysis.relativeL2Threshold = threshold;
analysis.includeInSurfacePattern = includeCorrection;
analysis.exactlySeparableBoundaryModel = ...
    exactlySeparableBoundaryModel;
analysis.nonseparableCorrectionInterpretation = ...
    correction_interpretation(exactlySeparableBoundaryModel);
analysis.temporalPhaseSamples = tau;
analysis.fullCarrierSamples = fullCarrier;
analysis.besselCarrierSamples = besselCarrier;
analysis.correctionCarrierSamples = correctionCarrier;
analysis.timeResolvedRelativeCorrectionL2 = ...
    timeResolvedRelativeCorrection;
end

function value = correction_interpretation(exactlySeparableBoundaryModel)
if exactlySeparableBoundaryModel
    value = [ ...
        'numerical recovery/discretization departure from an exactly ', ...
        'separable Bessel interface; suppressed by the auto policy'];
else
    value = 'physical nonseparable radial correction in the selected model';
end
end

function weights = cylindrical_weights(metadata,r)
weights = [];
if isfield(metadata.discretization,'radial') && ...
        isfield(metadata.discretization.radial,'quadratureWeights')
    weights = metadata.discretization.radial.quadratureWeights(:);
end
if isempty(weights)
    weights = trapezoidal_weights(r);
end
if numel(weights) ~= numel(r)
    error('vi_radial_mode_correction:QuadratureSize', ...
        'The radial quadrature weights have the wrong length.');
end
weights = real(weights).*r;
negativeTolerance = 1.0e-12*max(max(abs(weights)),1);
if any(weights < -negativeTolerance)
    error('vi_radial_mode_correction:NegativeQuadrature', ...
        'The cylindrical radial quadrature contains negative weights.');
end
weights = max(weights,0);
end

function values = weighted_column_norms(matrix,weights)
values = sqrt(max(real(sum(conj(matrix).*(weights.*matrix),1)),0));
end

function value = weighted_frobenius_norm(matrix,weights)
value = norm(weighted_column_norms(matrix,weights));
end

function weights = trapezoidal_weights(r)
[sortedR,order] = sort(r);
sortedWeights = zeros(size(sortedR));
sortedWeights(1) = (sortedR(2)-sortedR(1))/2;
sortedWeights(end) = (sortedR(end)-sortedR(end-1))/2;
sortedWeights(2:end-1) = ...
    (sortedR(3:end)-sortedR(1:end-2))/2;
weights = zeros(size(r));
weights(order) = sortedWeights;
end

function value = setting(settings,name,defaultValue)
if isfield(settings,name) && ~isempty(settings.(name))
    value = settings.(name);
else
    value = defaultValue;
end
end
