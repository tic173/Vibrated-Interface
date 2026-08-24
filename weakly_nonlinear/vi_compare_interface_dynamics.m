function result = vi_compare_interface_dynamics(linearResult, ...
        wnlAmplitudeOverH, wnlResult, metadata, parameters, ...
        settings, amplitudeScale)
%VI_COMPARE_INTERFACE_DYNAMICS Reconstruct the same interface observable.
%
% The linear prediction is
%
%   zeta_L/h = Re{a0 exp(gamma*t) Phi_L(r,t) exp(i*m*theta)}.
%
% The weakly nonlinear prediction uses the Landau amplitude A(t) and the
% full-cylinder neutral interface mode. When available, it also includes
% the slaved O(A^2) second-azimuthal-harmonic and mean-interface fields,
%
%   zeta_WNL/h = Re{A*phi + A^2*q_AA + |A|^2*q_AbarA}.
%
% The cubic coefficient affects A(t). A third-order correction to the
% reconstructed spatial shape is not included because it is not needed for
% the cubic solvability coefficient and is not solved by this module.

requiredLinear = {'timeStar', 'forcingPeriods', 'displacementOverH', ...
    'floquetOscillation', 'betaStar'};
for fieldIndex = 1:numel(requiredLinear)
    assert(isfield(linearResult, requiredLinear{fieldIndex}), ...
        'linearResult.%s is required.', requiredLinear{fieldIndex});
end
assert(isfield(metadata, 'layout') && ...
    isfield(metadata.layout, 'zeta'), ...
    'metadata.layout.zeta is required.');
assert(isfield(metadata, 'discretization') && ...
    isfield(metadata.discretization, 'r'), ...
    'metadata.discretization.r is required.');
assert(isfield(parameters, 'omegaStar') && isfield(parameters, 'R0'), ...
    'parameters.omegaStar and parameters.R0 are required.');
validateattributes(amplitudeScale, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'});

[mode, selfResult] = first_mode_and_self(wnlResult);
timeStar = linearResult.timeStar(:).';
forcingPeriods = linearResult.forcingPeriods(:).';
wnlAmplitudeOverH = wnlAmplitudeOverH(:).';
assert(numel(wnlAmplitudeOverH) == numel(timeStar), ...
    'The WNL amplitude and linear time vectors must have the same length.');

radialGrid = metadata.discretization.r(:);
[radialGrid, radialOrder] = sort(radialGrid, 'ascend');
zetaRows = metadata.layout.zeta;

primaryCarrier = interface_carrier(mode.field, zetaRows, ...
    parameters.omegaStar, timeStar);
primaryCarrier = primaryCarrier(radialOrder, :);

linearRadialShape = besselj(mode.spec.m, ...
    linearResult.betaStar*radialGrid);
linearRadialScale = max(abs(linearRadialShape));
if linearRadialScale <= eps
    error('vi_compare_interface_dynamics:ZeroLinearRadialMode', ...
        'The reduced linear radial interface mode is zero on the grid.');
end
linearRadialShape = linearRadialShape/linearRadialScale;
linearModal = linearRadialShape*linearResult.displacementOverH(:).';

% A direct eigenvector has an arbitrary complex phase. Align the WNL
% critical carrier with the reduced linear carrier over the first forcing
% period so both predictions represent the same initial interface state.
% The Landau equation is phase equivariant, so this is exactly equivalent
% to choosing the corresponding phase of its initial amplitude. The A^2
% slaved field then receives twice this phase automatically.
linearCarrierForAlignment = linearRadialShape * ...
    linearResult.floquetOscillation(:).';
firstPeriod = forcingPeriods <= forcingPeriods(1)+1+64*eps;
alignmentProduct = conj(primaryCarrier(:, firstPeriod)) .* ...
    linearCarrierForAlignment(:, firstPeriod);
alignmentOverlap = sum(alignmentProduct(:));
alignmentDenominator = norm(primaryCarrier(:, firstPeriod), 'fro') * ...
    norm(linearCarrierForAlignment(:, firstPeriod), 'fro');
if abs(alignmentOverlap) <= eps*max(alignmentDenominator, 1)
    alignmentPhase = 1.0;
    normalizedAlignmentOverlap = 0.0;
    warning('vi_compare_interface_dynamics:PhaseAlignmentFailed', ...
        ['The linear and WNL primary interface carriers have nearly zero ', ...
         'overlap. Their arbitrary eigenvector phases were not aligned.']);
else
    alignmentPhase = alignmentOverlap/abs(alignmentOverlap);
    normalizedAlignmentOverlap = ...
        abs(alignmentOverlap)/alignmentDenominator;
end
modalAmplitude = ...
    (wnlAmplitudeOverH/amplitudeScale)*alignmentPhase;
wnlPrimaryModal = primaryCarrier.*modalAmplitude;

includeSlaved = option(settings, 'includeSlavedHarmonics', true);
secondHarmonicModal = complex(zeros(size(wnlPrimaryModal)));
meanModal = complex(zeros(size(wnlPrimaryModal)));
secondHarmonicM = 2*mode.spec.m;
meanM = 0;
includedSecondHarmonic = false;
includedMean = false;
if includeSlaved && isstruct(selfResult)
    if isfield(selfResult, 'qAA') && isstruct(selfResult.qAA) && ...
            isfield(selfResult.qAA, 'field') && ...
            ~isempty(selfResult.qAA.field)
        carrierAA = interface_carrier(selfResult.qAA.field, zetaRows, ...
            parameters.omegaStar, timeStar);
        carrierAA = carrierAA(radialOrder, :);
        secondHarmonicModal = carrierAA.*(modalAmplitude.^2);
        secondHarmonicM = selfResult.qAA.field.spec.m;
        includedSecondHarmonic = true;
    end
    if isfield(selfResult, 'qAbarA') && ...
            isstruct(selfResult.qAbarA) && ...
            isfield(selfResult.qAbarA, 'field') && ...
            ~isempty(selfResult.qAbarA.field)
        carrierMean = interface_carrier( ...
            selfResult.qAbarA.field, zetaRows, ...
            parameters.omegaStar, timeStar);
        carrierMean = carrierMean(radialOrder, :);
        meanModal = carrierMean.*abs(modalAmplitude).^2;
        meanM = selfResult.qAbarA.field.spec.m;
        includedMean = true;
    end
end

probeTheta = option(settings, 'probeTheta', 0.0);
validateattributes(probeTheta, {'numeric'}, ...
    {'scalar', 'real', 'finite'});
azimuthPrimary = exp(1i*mode.spec.m*probeTheta);
azimuthSecond = exp(1i*secondHarmonicM*probeTheta);
azimuthMean = exp(1i*meanM*probeTheta);

linearComplexAtTheta = linearModal*azimuthPrimary;
wnlPrimaryComplexAtTheta = wnlPrimaryModal*azimuthPrimary;
wnlSecondComplexAtTheta = ...
    secondHarmonicModal*azimuthSecond + meanModal*azimuthMean;
wnlTotalComplexAtTheta = ...
    wnlPrimaryComplexAtTheta + wnlSecondComplexAtTheta;

linearPhysical = real(linearComplexAtTheta);
wnlPrimaryPhysical = real(wnlPrimaryComplexAtTheta);
wnlSecondPhysical = real(wnlSecondComplexAtTheta);
wnlTotalPhysical = real(wnlTotalComplexAtTheta);
difference = linearPhysical-wnlTotalPhysical;

probeRadiusOverR = option(settings, 'probeRadiusOverR', []);
if isempty(probeRadiusOverR)
    [~, probeIndex] = max(max(abs(primaryCarrier), [], 2));
    probeSelection = 'automatic maximum of the critical interface mode';
else
    validateattributes(probeRadiusOverR, {'numeric'}, ...
        {'scalar', 'real', '>=', 0, '<=', 1, 'finite'});
    [~, probeIndex] = min(abs( ...
        radialGrid/parameters.R0-probeRadiusOverR));
    probeSelection = 'nearest radial collocation point';
end

linearProbe = linearPhysical(probeIndex, :).';
wnlPrimaryProbe = wnlPrimaryPhysical(probeIndex, :).';
wnlSecondProbe = wnlSecondPhysical(probeIndex, :).';
wnlTotalProbe = wnlTotalPhysical(probeIndex, :).';
probeDifference = linearProbe-wnlTotalProbe;

fieldRmse = sqrt(mean(difference(:).^2));
fieldRelativeL2 = norm(difference(:)) / ...
    max(norm(linearPhysical(:)), 1.0e-12);
fieldSymmetricRelativeL2 = norm(difference(:)) / ...
    max([norm(linearPhysical(:)), norm(wnlTotalPhysical(:)), 1.0e-12]);
fieldMaximumDifference = max(abs(difference(:)));
probeRmse = sqrt(mean(probeDifference.^2));
probeRelativeL2 = norm(probeDifference) / ...
    max(norm(linearProbe), 1.0e-12);
probeMaximumDifference = max(abs(probeDifference));
secondOrderRelativeL2 = norm(wnlSecondPhysical(:)) / ...
    max(norm(wnlPrimaryPhysical(:)), 1.0e-12);

radialDifferenceNorm = sqrt(sum(difference.^2, 1));
radialReferenceNorm = max(sqrt(sum(linearPhysical.^2, 1)), ...
    sqrt(sum(wnlTotalPhysical.^2, 1)));
timeResolvedSymmetricError = radialDifferenceNorm ./ ...
    max(radialReferenceNorm, 1.0e-12);

snapshotForcingPeriod = option(settings, 'snapshotForcingPeriod', []);
if isempty(snapshotForcingPeriod)
    [~, snapshotIndex] = max(radialDifferenceNorm);
    snapshotSelection = 'maximum radial L2 difference';
else
    validateattributes(snapshotForcingPeriod, {'numeric'}, ...
        {'scalar', 'real', 'finite', 'nonnegative'});
    [~, snapshotIndex] = min(abs( ...
        forcingPeriods-snapshotForcingPeriod));
    snapshotSelection = 'nearest requested forcing period';
end

numberOfTheta = option(settings, 'snapshotAzimuthalSamples', 129);
validateattributes(numberOfTheta, {'numeric'}, ...
    {'scalar', 'integer', '>=', 33});
thetaSnapshot = linspace(0, 2*pi, numberOfTheta);
linearSnapshot = real(linearModal(:, snapshotIndex) * ...
    exp(1i*mode.spec.m*thetaSnapshot));
wnlPrimarySnapshot = real(wnlPrimaryModal(:, snapshotIndex) * ...
    exp(1i*mode.spec.m*thetaSnapshot));
wnlSecondSnapshot = real( ...
    secondHarmonicModal(:, snapshotIndex) * ...
    exp(1i*secondHarmonicM*thetaSnapshot) + ...
    meanModal(:, snapshotIndex) * exp(1i*meanM*thetaSnapshot));
wnlTotalSnapshot = wnlPrimarySnapshot+wnlSecondSnapshot;
snapshotDifference = linearSnapshot-wnlTotalSnapshot;

result = struct();
result.nonlinearModelLabel = option( ...
    settings, 'nonlinearModelLabel', 'WNL');
result.description = [ ...
    'Real interface reconstructed with one common peak-amplitude ', ...
    'normalization; WNL contains its Landau-modulated critical mode ', ...
    'and available O(A^2) slaved interface harmonics.'];
result.realFieldConvention = ...
    'zeta/h = real{complex modal representation}';
result.timeStar = timeStar(:);
result.forcingPeriods = forcingPeriods(:);
result.r = radialGrid;
result.rOverR = radialGrid/parameters.R0;
result.probeTheta = probeTheta;
result.probeIndex = probeIndex;
result.probeRadius = radialGrid(probeIndex);
result.probeRadiusOverR = radialGrid(probeIndex)/parameters.R0;
result.probeSelection = probeSelection;
result.phaseAlignment = alignmentPhase;
result.normalizedPrimaryCarrierOverlap = normalizedAlignmentOverlap;
result.linearField = linearPhysical;
result.wnlPrimaryField = wnlPrimaryPhysical;
result.wnlSecondOrderField = wnlSecondPhysical;
result.wnlTotalField = wnlTotalPhysical;
result.differenceLinearMinusWnl = difference;
result.linearProbeSignal = linearProbe;
result.wnlPrimaryProbeSignal = wnlPrimaryProbe;
result.wnlSecondOrderProbeSignal = wnlSecondProbe;
result.wnlTotalProbeSignal = wnlTotalProbe;
result.probeDifferenceLinearMinusWnl = probeDifference;
result.fieldRmse = fieldRmse;
result.fieldRelativeL2 = fieldRelativeL2;
result.fieldSymmetricRelativeL2 = fieldSymmetricRelativeL2;
result.fieldMaximumAbsoluteDifference = fieldMaximumDifference;
result.probeRmse = probeRmse;
result.probeRelativeL2 = probeRelativeL2;
result.probeMaximumAbsoluteDifference = probeMaximumDifference;
result.secondOrderRelativeL2 = secondOrderRelativeL2;
result.timeResolvedSymmetricRelativeL2 = ...
    timeResolvedSymmetricError(:);
result.includedSecondAzimuthalHarmonic = includedSecondHarmonic;
result.includedMeanInterfaceCorrection = includedMean;
result.primaryAzimuthalWavenumber = mode.spec.m;
result.secondAzimuthalWavenumber = secondHarmonicM;
result.meanAzimuthalWavenumber = meanM;
result.snapshotIndex = snapshotIndex;
result.snapshotForcingPeriod = forcingPeriods(snapshotIndex);
result.snapshotSelection = snapshotSelection;
result.snapshotTheta = thetaSnapshot;
result.linearSnapshot = linearSnapshot;
result.wnlPrimarySnapshot = wnlPrimarySnapshot;
result.wnlSecondOrderSnapshot = wnlSecondSnapshot;
result.wnlTotalSnapshot = wnlTotalSnapshot;
result.snapshotDifferenceLinearMinusWnl = snapshotDifference;

fprintf('\nReconstructed interfacial-dynamics comparison\n');
fprintf('  probe r/R, theta            = %.6g, %.6g rad\n', ...
    result.probeRadiusOverR, probeTheta);
fprintf('  included WNL m=%d harmonic  = %s\n', ...
    secondHarmonicM, yes_no(includedSecondHarmonic));
fprintf('  included WNL mean field     = %s\n', ...
    yes_no(includedMean));
fprintf('  aligned carrier overlap     = %.6e\n', ...
    normalizedAlignmentOverlap);
fprintf('  field RMSE                  = %.6e [zeta/h]\n', ...
    fieldRmse);
fprintf('  field relative L2 error     = %.6e\n', fieldRelativeL2);
fprintf('  field symmetric relative L2 = %.6e\n', ...
    fieldSymmetricRelativeL2);
fprintf('  max |field linear - WNL|    = %.6e [zeta/h]\n', ...
    fieldMaximumDifference);
fprintf('  probe RMSE                  = %.6e [zeta/h]\n', ...
    probeRmse);
fprintf('  O(A^2)/primary relative L2  = %.6e\n', ...
    secondOrderRelativeL2);

if option(settings, 'plotInterfaceDynamics', true)
    plot_radial_time_comparison(result);
    plot_interface_snapshot(result, parameters.R0);
end
end

function carrier = interface_carrier(field, zetaRows, omegaStar, timeStar)
frequency = field.spec.n(:)+field.spec.s;
phase = exp(1i*frequency*(omegaStar*timeStar));
carrier = field.coeff(zetaRows, :)*phase;
end

function [mode, selfResult] = first_mode_and_self(wnlResult)
if isfield(wnlResult, 'mode')
    mode = wnlResult.mode;
    selfResult = wnlResult.self;
elseif isfield(wnlResult, 'modes') && ~isempty(wnlResult.modes)
    mode = wnlResult.modes{1};
    selfResult = wnlResult.self{1};
else
    error('vi_compare_interface_dynamics:MissingMode', ...
        'The WNL result does not contain a retained direct mode.');
end
end

function value = option(settings, name, defaultValue)
if isfield(settings, name)
    value = settings.(name);
else
    value = defaultValue;
end
end

function value = yes_no(flag)
if flag
    value = 'yes';
else
    value = 'no';
end
end

function plot_radial_time_comparison(result)
periods = result.forcingPeriods;
rOverR = result.rOverR;
modelLabel = result.nonlinearModelLabel;
fieldLimit = max(abs([result.linearField(:); ...
    result.wnlTotalField(:)]));
fieldLimit = max(fieldLimit, 1.0e-12);
differenceLimit = max(abs(result.differenceLinearMinusWnl(:)));
differenceLimit = max(differenceLimit, 1.0e-12);

figure('Name',sprintf('Linear versus %s interfacial dynamics',modelLabel));
tiledlayout(3, 2, 'TileSpacing', 'compact');

nexttile;
plot(periods, result.linearProbeSignal, 'LineWidth', 1.2);
hold on;
plot(periods, result.wnlPrimaryProbeSignal, '--', 'LineWidth', 1.2);
plot(periods, result.wnlTotalProbeSignal, 'LineWidth', 1.4);
grid on;
xlabel('forcing periods, t/T_f');
ylabel('\zeta/h');
legend('linear',sprintf('%s primary',modelLabel), ...
    sprintf('%s primary + slaved O(A^2)',modelLabel), ...
    'Location', 'best');
title(sprintf('probe r/R=%.3g, \theta=%.3g', ...
    result.probeRadiusOverR, result.probeTheta));

nexttile;
plot(periods, result.probeDifferenceLinearMinusWnl, ...
    'LineWidth', 1.2);
hold on;
yline(0, ':');
grid on;
xlabel('forcing periods, t/T_f');
ylabel(sprintf('(linear-%s) \\zeta/h',modelLabel));
title(sprintf('probe difference, RMSE=%.3g', result.probeRmse));

nexttile;
imagesc(periods, rOverR, result.linearField);
axis xy;
caxis([-fieldLimit, fieldLimit]);
colorbar;
xlabel('forcing periods, t/T_f');
ylabel('r/R');
title('linear interface');

nexttile;
imagesc(periods, rOverR, result.wnlTotalField);
axis xy;
caxis([-fieldLimit, fieldLimit]);
colorbar;
xlabel('forcing periods, t/T_f');
ylabel('r/R');
title(sprintf('%s interface: primary + slaved O(A^2)',modelLabel));

nexttile;
imagesc(periods, rOverR, result.wnlSecondOrderField);
axis xy;
colorbar;
xlabel('forcing periods, t/T_f');
ylabel('r/R');
title(sprintf('%s O(A^2) correction, relative L_2=%.3g', ...
    modelLabel,result.secondOrderRelativeL2));

nexttile;
imagesc(periods, rOverR, result.differenceLinearMinusWnl);
axis xy;
caxis([-differenceLimit, differenceLimit]);
colorbar;
xlabel('forcing periods, t/T_f');
ylabel('r/R');
title(sprintf('linear-%s, relative L_2=%.3g', ...
    modelLabel,result.fieldRelativeL2));
end

function plot_interface_snapshot(result, cylinderRadius)
modelLabel = result.nonlinearModelLabel;
theta = result.snapshotTheta;
r = result.r;
[thetaGrid, radialGrid] = meshgrid(theta, r);
x = radialGrid.*cos(thetaGrid)/cylinderRadius;
y = radialGrid.*sin(thetaGrid)/cylinderRadius;
fieldLimit = max(abs([result.linearSnapshot(:); ...
    result.wnlTotalSnapshot(:)]));
fieldLimit = max(fieldLimit, 1.0e-12);
differenceLimit = max(abs( ...
    result.snapshotDifferenceLinearMinusWnl(:)));
differenceLimit = max(differenceLimit, 1.0e-12);

figure('Name',sprintf('Linear versus %s interface snapshot',modelLabel));
tiledlayout(1, 3, 'TileSpacing', 'compact');

nexttile;
surf(x, y, result.linearSnapshot, 'EdgeColor', 'none');
view(2);
axis equal tight;
caxis([-fieldLimit, fieldLimit]);
colorbar;
xlabel('x/R');
ylabel('y/R');
title('linear \zeta/h');

nexttile;
surf(x, y, result.wnlTotalSnapshot, 'EdgeColor', 'none');
view(2);
axis equal tight;
caxis([-fieldLimit, fieldLimit]);
colorbar;
xlabel('x/R');
ylabel('y/R');
title(sprintf('%s \\zeta/h',modelLabel));

nexttile;
surf(x, y, result.snapshotDifferenceLinearMinusWnl, ...
    'EdgeColor', 'none');
view(2);
axis equal tight;
caxis([-differenceLimit, differenceLimit]);
colorbar;
xlabel('x/R');
ylabel('y/R');
title(sprintf('linear-%s',modelLabel));
sgtitle(sprintf('interface at t/T_f=%.4g (%s)', ...
    result.snapshotForcingPeriod, result.snapshotSelection));
end
