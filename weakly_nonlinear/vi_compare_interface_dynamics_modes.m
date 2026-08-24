function result = vi_compare_interface_dynamics_modes(linearResult, ...
        modeAmplitudesOverH, wnlResult, metadata, parameters, ...
        settings, amplitudeScale, linearModeAmplitudesOverH)
%VI_COMPARE_INTERFACE_DYNAMICS_MODES One/two-mode interface reconstruction.
%
% For one mode this calls vi_compare_interface_dynamics. For two modes it
% reconstructs both primary modes, every available self-generated mean and
% second harmonic, and the sum/difference fields q_AB and q_AbarB. An
% optional eighth argument supplies the exact linear amplitudes of both
% modes; this is used by the operating-point driver.

if nargin < 8
    linearModeAmplitudesOverH = [];
end

if isfield(wnlResult, 'mode')
    numberOfModes = 1;
else
    assert(isfield(wnlResult, 'modes') && ...
        numel(wnlResult.modes) == 2, ...
        'A two-mode WNL result must contain exactly two modes.');
    numberOfModes = 2;
end
numberOfTimes = numel(linearResult.timeStar);
amplitudes = orient_amplitudes( ...
    modeAmplitudesOverH, numberOfTimes, numberOfModes);
if numberOfModes == 1
    result = vi_compare_interface_dynamics(linearResult, ...
        amplitudes(:, 1), wnlResult, metadata, parameters, ...
        settings, amplitudeScale);
    result.numberOfModes = 1;
    result.modeAmplitudesOverH = amplitudes;
    return;
end

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
modes = wnlResult.modes;
selfResults = wnlResult.self;
crossResult = wnlResult.cross{1, 2};
timeStar = linearResult.timeStar(:).';
forcingPeriods = linearResult.forcingPeriods(:).';
radialGrid = metadata.discretization.r(:);
[radialGrid, radialOrder] = sort(radialGrid, 'ascend');
zetaRows = metadata.layout.zeta;

carriers = cell(2, 1);
for modeIndex = 1:2
    carriers{modeIndex} = interface_carrier( ...
        modes{modeIndex}.field, zetaRows, ...
        parameters.omegaStar, timeStar);
    carriers{modeIndex} = carriers{modeIndex}(radialOrder, :);
end

linearRadialShape = besselj(modes{1}.spec.m, ...
    linearResult.betaStar*radialGrid);
linearRadialScale = max(abs(linearRadialShape));
if linearRadialScale <= eps
    error('vi_compare_interface_dynamics_modes:ZeroLinearRadialMode', ...
        'The reduced first-mode radial interface shape is zero on the grid.');
end
linearRadialShape = linearRadialShape/linearRadialScale;
linearModal = linearRadialShape*linearResult.displacementOverH(:).';
linearCarrier = linearRadialShape* ...
    linearResult.floquetOscillation(:).';
firstPeriod = forcingPeriods <= forcingPeriods(1)+1+64*eps;
alignmentProduct = conj(carriers{1}(:, firstPeriod)).* ...
    linearCarrier(:, firstPeriod);
overlap = sum(alignmentProduct(:));
overlapScale = norm(carriers{1}(:, firstPeriod), 'fro')* ...
    norm(linearCarrier(:, firstPeriod), 'fro');
if abs(overlap) <= eps*max(overlapScale, 1)
    alignment = 1;
    normalizedOverlap = 0;
else
    alignment = overlap/abs(overlap);
    normalizedOverlap = abs(overlap)/overlapScale;
end

modalAmplitudes = amplitudes/amplitudeScale;
modalAmplitudes(:, 1) = modalAmplitudes(:, 1)*alignment;
linearModeAmplitudes = [];
linearTwoModeAmplitudes = [];
if ~isempty(linearModeAmplitudesOverH)
    linearModeAmplitudes = orient_amplitudes( ...
        linearModeAmplitudesOverH, numberOfTimes, numberOfModes);
    linearTwoModeAmplitudes = linearModeAmplitudes/amplitudeScale;
    linearTwoModeAmplitudes(:, 1) = ...
        linearTwoModeAmplitudes(:, 1)*alignment;
end
primaryModal = cell(2, 1);
for modeIndex = 1:2
    primaryModal{modeIndex} = carriers{modeIndex}.* ...
        modalAmplitudes(:, modeIndex).';
end

terms = struct('modal', {}, 'm', {}, 'label', {});
if get_option(settings, 'includeSlavedHarmonics', true)
    for modeIndex = 1:2
        selfValue = selfResults{modeIndex};
        if has_forced_field(selfValue, 'qAA')
            carrier = forced_carrier(selfValue.qAA, zetaRows, ...
                parameters.omegaStar, timeStar, radialOrder);
            terms(end+1) = make_term(carrier.* ... %#ok<AGROW>
                (modalAmplitudes(:, modeIndex).'.^2), ...
                selfValue.qAA.field.spec.m, ...
                sprintf('mode %d self second harmonic', modeIndex));
        end
        if has_forced_field(selfValue, 'qAbarA')
            carrier = forced_carrier(selfValue.qAbarA, zetaRows, ...
                parameters.omegaStar, timeStar, radialOrder);
            terms(end+1) = make_term(carrier.* ... %#ok<AGROW>
                (abs(modalAmplitudes(:, modeIndex).').^2), ...
                selfValue.qAbarA.field.spec.m, ...
                sprintf('mode %d mean', modeIndex));
        end
    end
    if isstruct(crossResult) && has_forced_field(crossResult, 'qAB')
        carrier = forced_carrier(crossResult.qAB, zetaRows, ...
            parameters.omegaStar, timeStar, radialOrder);
        terms(end+1) = make_term(carrier.* ... %#ok<AGROW>
            (modalAmplitudes(:, 1).*modalAmplitudes(:, 2)).', ...
            crossResult.qAB.field.spec.m, 'mode sum interaction');
    end
    if isstruct(crossResult) && has_forced_field(crossResult, 'qAbarB')
        carrier = forced_carrier(crossResult.qAbarB, zetaRows, ...
            parameters.omegaStar, timeStar, radialOrder);
        terms(end+1) = make_term(carrier.* ... %#ok<AGROW>
            (modalAmplitudes(:, 1).*conj(modalAmplitudes(:, 2))).', ...
            crossResult.qAbarB.field.spec.m, 'mode difference interaction');
    end
end

probeTheta = get_option(settings, 'probeTheta', 0);
validateattributes(probeTheta, {'numeric'}, ...
    {'scalar', 'real', 'finite'});
if isempty(linearTwoModeAmplitudes)
    linearComplex = linearModal*exp(1i*modes{1}.spec.m*probeTheta);
else
    linearComplex = complex(zeros(size(linearModal)));
    for modeIndex = 1:2
        linearComplex = linearComplex + carriers{modeIndex}.* ...
            linearTwoModeAmplitudes(:, modeIndex).' * ...
            exp(1i*modes{modeIndex}.spec.m*probeTheta);
    end
end
primaryComplex = complex(zeros(size(linearComplex)));
for modeIndex = 1:2
    primaryComplex = primaryComplex + primaryModal{modeIndex}* ...
        exp(1i*modes{modeIndex}.spec.m*probeTheta);
end
secondComplex = complex(zeros(size(linearComplex)));
for termIndex = 1:numel(terms)
    secondComplex = secondComplex + terms(termIndex).modal* ...
        exp(1i*terms(termIndex).m*probeTheta);
end
totalComplex = primaryComplex+secondComplex;
linearPhysical = real(linearComplex);
primaryPhysical = real(primaryComplex);
secondPhysical = real(secondComplex);
totalPhysical = real(totalComplex);
difference = linearPhysical-totalPhysical;

probeRadiusOverR = get_option(settings, 'probeRadiusOverR', []);
if isempty(probeRadiusOverR)
    [~, probeIndex] = max(max(abs(carriers{1}), [], 2));
    probeSelection = 'automatic maximum of first critical mode';
else
    validateattributes(probeRadiusOverR, {'numeric'}, ...
        {'scalar', 'real', '>=', 0, '<=', 1, 'finite'});
    [~, probeIndex] = min(abs( ...
        radialGrid/parameters.R0-probeRadiusOverR));
    probeSelection = 'nearest radial collocation point';
end
wnlPrimaryModeProbeSignals = zeros(numberOfTimes, numberOfModes);
linearModeProbeSignals = [];
if ~isempty(linearTwoModeAmplitudes)
    linearModeProbeSignals = zeros(numberOfTimes, numberOfModes);
end
for modeIndex = 1:numberOfModes
    azimuthalPhase = exp(1i*modes{modeIndex}.spec.m*probeTheta);
    wnlPrimaryModeProbeSignals(:,modeIndex) = real( ...
        primaryModal{modeIndex}(probeIndex,:)*azimuthalPhase).';
    if ~isempty(linearTwoModeAmplitudes)
        linearModeProbeSignals(:,modeIndex) = real( ...
            carriers{modeIndex}(probeIndex,:).* ...
            linearTwoModeAmplitudes(:,modeIndex).' * ...
            azimuthalPhase).';
    end
end
linearProbe = linearPhysical(probeIndex, :).';
primaryProbe = primaryPhysical(probeIndex, :).';
secondProbe = secondPhysical(probeIndex, :).';
totalProbe = totalPhysical(probeIndex, :).';
probeDifference = linearProbe-totalProbe;

fieldRmse = sqrt(mean(difference(:).^2));
fieldRelativeL2 = norm(difference(:))/ ...
    max(norm(linearPhysical(:)), 1e-12);
fieldSymmetricRelativeL2 = norm(difference(:))/max( ...
    [norm(linearPhysical(:)), norm(totalPhysical(:)), 1e-12]);
radialDifferenceNorm = sqrt(sum(difference.^2, 1));
radialReferenceNorm = max(sqrt(sum(linearPhysical.^2, 1)), ...
    sqrt(sum(totalPhysical.^2, 1)));
timeResolvedError = radialDifferenceNorm./ ...
    max(radialReferenceNorm, 1e-12);

snapshotPeriod = get_option(settings, 'snapshotForcingPeriod', []);
if isempty(snapshotPeriod)
    [~, snapshotIndex] = max(radialDifferenceNorm);
    snapshotSelection = 'maximum radial L2 difference';
else
    validateattributes(snapshotPeriod, {'numeric'}, ...
        {'scalar', 'real', 'finite'});
    [~, snapshotIndex] = min(abs(forcingPeriods-snapshotPeriod));
    snapshotSelection = 'nearest requested forcing period';
end
numberOfTheta = get_option(settings, 'snapshotAzimuthalSamples', 129);
validateattributes(numberOfTheta, {'numeric'}, ...
    {'scalar', 'integer', '>=', 3});
theta = linspace(0, 2*pi, numberOfTheta);
if isempty(linearTwoModeAmplitudes)
    linearSnapshot = real(linearModal(:, snapshotIndex)* ...
        exp(1i*modes{1}.spec.m*theta));
else
    linearSnapshot = complex(zeros(numel(radialGrid), numberOfTheta));
    for modeIndex = 1:2
        linearSnapshot = linearSnapshot + ...
            carriers{modeIndex}(:, snapshotIndex) * ...
            linearTwoModeAmplitudes(snapshotIndex, modeIndex) * ...
            exp(1i*modes{modeIndex}.spec.m*theta);
    end
    linearSnapshot = real(linearSnapshot);
end
primarySnapshot = complex(zeros(numel(radialGrid), numberOfTheta));
for modeIndex = 1:2
    primarySnapshot = primarySnapshot + ...
        primaryModal{modeIndex}(:, snapshotIndex)* ...
        exp(1i*modes{modeIndex}.spec.m*theta);
end
secondSnapshot = complex(zeros(size(primarySnapshot)));
for termIndex = 1:numel(terms)
    secondSnapshot = secondSnapshot + ...
        terms(termIndex).modal(:, snapshotIndex)* ...
        exp(1i*terms(termIndex).m*theta);
end
primarySnapshot = real(primarySnapshot);
secondSnapshot = real(secondSnapshot);
totalSnapshot = primarySnapshot+secondSnapshot;

result = struct();
result.nonlinearModelLabel = get_option( ...
    settings, 'nonlinearModelLabel', 'WNL');
result.description = ['Two-mode real interface: both primary modes plus ', ...
    'available self and cross O(A^2) slaved fields.'];
result.realFieldConvention = ...
    'zeta/h = real{complex modal representation}';
result.numberOfModes = 2;
result.modeAmplitudesOverH = amplitudes;
result.linearModeAmplitudesOverH = linearModeAmplitudes;
result.realModeAmplitudesWnl = real(amplitudes);
result.realModeAmplitudesLinear = real(linearModeAmplitudes);
if isempty(linearModeAmplitudes)
    result.realModeAmplitudeDifferenceLinearMinusWnl = [];
else
    result.realModeAmplitudeDifferenceLinearMinusWnl = ...
        real(linearModeAmplitudes)-real(amplitudes);
end
result.linearReferenceIncludesBothModes = ...
    ~isempty(linearTwoModeAmplitudes);
result.modeLabels = {modes{1}.spec.label; modes{2}.spec.label};
result.primaryAzimuthalWavenumbers = ...
    [modes{1}.spec.m; modes{2}.spec.m];
result.secondOrderAzimuthalWavenumbers = [terms.m].';
result.secondOrderLabels = {terms.label}.';
result.timeStar = timeStar(:);
result.forcingPeriods = forcingPeriods(:);
result.r = radialGrid;
result.rOverR = radialGrid/parameters.R0;
result.probeTheta = probeTheta;
result.probeIndex = probeIndex;
result.probeRadius = radialGrid(probeIndex);
result.probeRadiusOverR = radialGrid(probeIndex)/parameters.R0;
result.probeSelection = probeSelection;
result.phaseAlignment = [alignment; 1];
result.normalizedPrimaryCarrierOverlap = normalizedOverlap;
result.linearField = linearPhysical;
result.wnlPrimaryField = primaryPhysical;
result.wnlSecondOrderField = secondPhysical;
result.wnlTotalField = totalPhysical;
result.differenceLinearMinusWnl = difference;
result.linearProbeSignal = linearProbe;
result.linearModeProbeSignals = linearModeProbeSignals;
result.wnlPrimaryProbeSignal = primaryProbe;
result.wnlPrimaryModeProbeSignals = wnlPrimaryModeProbeSignals;
result.wnlSecondOrderProbeSignal = secondProbe;
result.wnlTotalProbeSignal = totalProbe;
result.probeDifferenceLinearMinusWnl = probeDifference;
result.fieldRmse = fieldRmse;
result.fieldRelativeL2 = fieldRelativeL2;
result.fieldSymmetricRelativeL2 = fieldSymmetricRelativeL2;
result.fieldMaximumAbsoluteDifference = max(abs(difference(:)));
result.probeRmse = sqrt(mean(probeDifference.^2));
result.probeRelativeL2 = norm(probeDifference)/ ...
    max(norm(linearProbe), 1e-12);
result.probeMaximumAbsoluteDifference = max(abs(probeDifference));
result.secondOrderRelativeL2 = norm(secondPhysical(:))/ ...
    max(norm(primaryPhysical(:)), 1e-12);
result.timeResolvedSymmetricRelativeL2 = timeResolvedError(:);
result.includedSecondAzimuthalHarmonic = ~isempty(terms);
result.includedMeanInterfaceCorrection = ...
    any([terms.m] == 0);
result.primaryAzimuthalWavenumber = modes{1}.spec.m;
result.secondAzimuthalWavenumber = 2*modes{1}.spec.m;
result.meanAzimuthalWavenumber = 0;
result.snapshotIndex = snapshotIndex;
result.snapshotForcingPeriod = forcingPeriods(snapshotIndex);
result.snapshotSelection = snapshotSelection;
result.snapshotTheta = theta;
result.linearSnapshot = linearSnapshot;
result.wnlPrimarySnapshot = primarySnapshot;
result.wnlSecondOrderSnapshot = secondSnapshot;
result.wnlTotalSnapshot = totalSnapshot;
result.snapshotDifferenceLinearMinusWnl = ...
    linearSnapshot-totalSnapshot;

if get_option(settings, 'verbose', true)
    fprintf('\nTwo-mode reconstructed interfacial dynamics\n');
    fprintf('  primary azimuthal modes      = m=%d, m=%d\n', ...
        modes{1}.spec.m, modes{2}.spec.m);
    fprintf('  retained O(A^2) fields       = %d\n', numel(terms));
    fprintf('  field relative L2 error      = %.6e\n', fieldRelativeL2);
    fprintf('  O(A^2)/primary relative L2   = %.6e\n', ...
        result.secondOrderRelativeL2);
end

if get_option(settings, 'plotInterfaceDynamics', true)
    plot_two_mode_fields(result, parameters.R0);
end
end

function amplitudes = orient_amplitudes(value, numberOfTimes, numberOfModes)
if isvector(value) && numberOfModes == 1
    amplitudes = value(:);
elseif isequal(size(value), [numberOfTimes, numberOfModes])
    amplitudes = value;
elseif isequal(size(value), [numberOfModes, numberOfTimes])
    amplitudes = value.';
else
    error('vi_compare_interface_dynamics_modes:AmplitudeSize', ...
        'Mode amplitudes must have size [%d,%d].', ...
        numberOfTimes, numberOfModes);
end
end

function carrier = interface_carrier(field, zetaRows, omegaStar, timeStar)
frequency = field.spec.n(:)+field.spec.s;
carrier = field.coeff(zetaRows, :)* ...
    exp(1i*frequency*(omegaStar*timeStar));
end

function carrier = forced_carrier(value, zetaRows, omegaStar, ...
        timeStar, radialOrder)
carrier = interface_carrier(value.field, zetaRows, omegaStar, timeStar);
carrier = carrier(radialOrder, :);
end

function flag = has_forced_field(value, name)
flag = isstruct(value) && isfield(value, name) && ...
    isstruct(value.(name)) && isfield(value.(name), 'field') && ...
    ~isempty(value.(name).field);
end

function term = make_term(modal, m, label)
term = struct('modal', modal, 'm', m, 'label', label);
end

function value = get_option(settings, name, defaultValue)
if isfield(settings, name)
    value = settings.(name);
else
    value = defaultValue;
end
end

function plot_two_mode_fields(result, cylinderRadius)
periods = result.forcingPeriods;
modelLabel = result.nonlinearModelLabel;
figure('Name',sprintf('Two-mode %s interfacial dynamics',modelLabel));
tiledlayout(2, 2, 'TileSpacing', 'compact');
nexttile;
hold on;
numberOfModes = result.numberOfModes;
modeColors = lines(numberOfModes);
if result.linearReferenceIncludesBothModes
    amplitudeHandles = gobjects(2*numberOfModes,1);
    amplitudeLabels = cell(2*numberOfModes,1);
    for modeIndex = 1:numberOfModes
        amplitudeHandles(2*modeIndex-1) = plot(periods, ...
            result.realModeAmplitudesLinear(:,modeIndex),'-', ...
            'Color',modeColors(modeIndex,:),'LineWidth',1.3);
        amplitudeHandles(2*modeIndex) = plot(periods, ...
            result.realModeAmplitudesWnl(:,modeIndex),'--', ...
            'Color',modeColors(modeIndex,:),'LineWidth',1.5);
        amplitudeLabels{2*modeIndex-1} = sprintf( ...
            'linear Re(A_%d): %s',modeIndex,result.modeLabels{modeIndex});
        amplitudeLabels{2*modeIndex} = sprintf( ...
            '%s Re(A_%d): %s',modelLabel,modeIndex, ...
            result.modeLabels{modeIndex});
    end
else
    amplitudeHandles = gobjects(numberOfModes,1);
    amplitudeLabels = cell(numberOfModes,1);
    for modeIndex = 1:numberOfModes
        amplitudeHandles(modeIndex) = plot(periods, ...
            result.realModeAmplitudesWnl(:,modeIndex),'-', ...
            'Color',modeColors(modeIndex,:),'LineWidth',1.3);
        amplitudeLabels{modeIndex} = sprintf( ...
            '%s Re(A_%d): %s',modelLabel,modeIndex, ...
            result.modeLabels{modeIndex});
    end
end
zeroHandle = yline(0,':');
zeroHandle.HandleVisibility = 'off';
grid on;
xlabel('forcing periods, t/T_f');
ylabel('Re(A_j), \zeta/h');
legend(amplitudeHandles,amplitudeLabels, ...
    'Interpreter','none','Location','best');
title(sprintf('linear and %s modal trajectories',modelLabel));
nexttile;
hold on;
if result.linearReferenceIncludesBothModes
    probeHandles = gobjects(2*numberOfModes,1);
    probeLabels = cell(2*numberOfModes,1);
    for modeIndex = 1:numberOfModes
        probeHandles(2*modeIndex-1) = plot(periods, ...
            result.linearModeProbeSignals(:,modeIndex),'-', ...
            'Color',modeColors(modeIndex,:),'LineWidth',1.2);
        probeHandles(2*modeIndex) = plot(periods, ...
            result.wnlPrimaryModeProbeSignals(:,modeIndex),'--', ...
            'Color',modeColors(modeIndex,:),'LineWidth',1.4);
        probeLabels{2*modeIndex-1} = sprintf( ...
            'linear mode %d: %s',modeIndex,result.modeLabels{modeIndex});
        probeLabels{2*modeIndex} = sprintf( ...
            '%s mode %d: %s',modelLabel,modeIndex, ...
            result.modeLabels{modeIndex});
    end
else
    probeHandles = plot(periods,result.wnlPrimaryModeProbeSignals, ...
        'LineWidth',1.3);
    probeLabels = cell(numberOfModes,1);
    for modeIndex = 1:numberOfModes
        probeLabels{modeIndex} = sprintf( ...
            '%s mode %d: %s',modelLabel,modeIndex, ...
            result.modeLabels{modeIndex});
    end
end
zeroHandle = yline(0,':');
zeroHandle.HandleVisibility = 'off';
grid on;
xlabel('forcing periods, t/T_f');
ylabel('\zeta/h at probe');
if result.linearReferenceIncludesBothModes
    differenceTitle = sprintf('two-mode linear - %s',modelLabel);
    snapshotLinearTitle = 'linear two-mode';
else
    differenceTitle = sprintf('first-mode linear - %s',modelLabel);
    snapshotLinearTitle = 'linear first mode';
end
legend(probeHandles,probeLabels,'Interpreter','none','Location','best');
title('mode-resolved primary probe contributions');
nexttile;
imagesc(periods, result.rOverR, result.wnlTotalField);
axis xy; colorbar;
xlabel('forcing periods, t/T_f'); ylabel('r/R');
title(sprintf('%s interface',modelLabel));
nexttile;
imagesc(periods, result.rOverR, result.differenceLinearMinusWnl);
axis xy; colorbar;
xlabel('forcing periods, t/T_f'); ylabel('r/R');
title(differenceTitle);

theta = result.snapshotTheta;
r = result.r;
[thetaGrid, radialGrid] = meshgrid(theta, r);
x = radialGrid.*cos(thetaGrid)/cylinderRadius;
y = radialGrid.*sin(thetaGrid)/cylinderRadius;
figure('Name',sprintf('Two-mode %s interface snapshot',modelLabel));
tiledlayout(1, 3, 'TileSpacing', 'compact');
fields = {result.linearSnapshot, result.wnlTotalSnapshot, ...
    result.snapshotDifferenceLinearMinusWnl};
titles = {snapshotLinearTitle,modelLabel,sprintf('linear - %s',modelLabel)};
for panel = 1:3
    nexttile;
    surf(x, y, fields{panel}, 'EdgeColor', 'none');
    view(2); axis equal tight; colorbar;
    xlabel('x/R'); ylabel('y/R'); title(titles{panel});
end
sgtitle(sprintf('interface at t/T_f=%.4g', ...
    result.snapshotForcingPeriod));
end
