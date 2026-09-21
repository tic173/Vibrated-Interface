function [figureHandle,summary] = ...
        vi_plot_bessel_corrected_structure(analysis,userSettings)
%VI_PLOT_BESSEL_CORRECTED_STRUCTURE Compare a Bessel carrier with c_n(r).
%
% The planforms use Re[zeta(r,tau)*exp(i*m*theta)] at the forcing phase
% where the recovered interface has its largest weighted radial norm. All
% planform and radial-profile panels share one normalization.

if nargin < 2
    userSettings = struct();
end
required = {'m','betaStar','rOverR','quadratureWeightsRdr', ...
    'relativeCorrectionL2','correctionEnergyFraction', ...
    'temporalPhaseSamples','fullCarrierSamples', ...
    'besselCarrierSamples','correctionCarrierSamples', ...
    'timeResolvedRelativeCorrectionL2'};
for fieldIndex = 1:numel(required)
    assert(isfield(analysis,required{fieldIndex}), ...
        'Radial-correction analysis is missing %s.',required{fieldIndex});
end

numberOfRadialPoints = setting(userSettings,'numberOfRadialPoints',240);
numberOfAzimuthalPoints = setting( ...
    userSettings,'numberOfAzimuthalPoints',360);
validateattributes(numberOfRadialPoints,{'numeric'}, ...
    {'scalar','integer','>=',64,'finite'});
validateattributes(numberOfAzimuthalPoints,{'numeric'}, ...
    {'scalar','integer','>=',120,'finite'});

weights = analysis.quadratureWeightsRdr(:);
fullSamples = analysis.fullCarrierSamples;
weightedNormSquared = real(sum(conj(fullSamples).* ...
    (weights.*fullSamples),1));
if isfield(userSettings,'forcingPhaseFraction') && ...
        ~isempty(userSettings.forcingPhaseFraction)
    phaseFraction = mod(userSettings.forcingPhaseFraction,1);
    availableFractions = analysis.temporalPhaseSamples(:).'/(2*pi);
    phaseDistance = abs(mod(availableFractions-phaseFraction+0.5,1)-0.5);
    [~,phaseIndex] = min(phaseDistance);
else
    [~,phaseIndex] = max(weightedNormSquared);
end

fullRadial = fullSamples(:,phaseIndex);
besselRadial = analysis.besselCarrierSamples(:,phaseIndex);
correctionRadial = analysis.correctionCarrierSamples(:,phaseIndex);
phaseRotation = orientation_rotation(besselRadial,fullRadial);
fullRadial = phaseRotation*fullRadial;
besselRadial = phaseRotation*besselRadial;
correctionRadial = phaseRotation*correctionRadial;

[radialNodes,order] = sort(real(analysis.rOverR(:)));
[radialNodes,uniquePosition] = unique(radialNodes,'stable');
order = order(uniquePosition);
fineRadius = linspace(max(0,min(radialNodes)),max(radialNodes), ...
    numberOfRadialPoints).';
fullFine = interpolate_complex(radialNodes,fullRadial(order),fineRadius);
besselFine = interpolate_complex( ...
    radialNodes,besselRadial(order),fineRadius);
correctionFine = interpolate_complex( ...
    radialNodes,correctionRadial(order),fineRadius);

theta = linspace(0,2*pi,numberOfAzimuthalPoints+1);
azimuthalPhase = exp(1i*analysis.m*theta);
fullPlanform = real(fullFine*azimuthalPhase);
besselPlanform = real(besselFine*azimuthalPhase);
correctionPlanform = real(correctionFine*azimuthalPhase);
[thetaGrid,radiusGrid] = meshgrid(theta,fineRadius);
xGrid = radiusGrid.*cos(thetaGrid);
yGrid = radiusGrid.*sin(thetaGrid);
normalization = max(abs([fullPlanform(:);besselPlanform(:); ...
    correctionPlanform(:)]));
normalization = max(normalization,eps);
fullPlanform = fullPlanform/normalization;
besselPlanform = besselPlanform/normalization;
correctionPlanform = correctionPlanform/normalization;

visibility = char(string(setting(userSettings,'visibility','on')));
figureHandle = figure('Name','Bessel and corrected Floquet structure', ...
    'Color','w','Visible',visibility,'Position',[100,100,1500,900]);
layout = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
plot_planform(nexttile(layout,1),xGrid,yGrid,besselPlanform, ...
    'Bessel projection, b_n J_m');
plot_planform(nexttile(layout,2),xGrid,yGrid,correctionPlanform, ...
    'nonseparable correction, c_n');
fullAxis = nexttile(layout,3);
plot_planform(fullAxis,xGrid,yGrid,fullPlanform, ...
    'full corrected structure');
colorbar(fullAxis);
colormap(figureHandle,blue_red_map(256));

profileAxis = nexttile(layout,4,[1,2]);
hold(profileAxis,'on');
plot(profileAxis,fineRadius,real(fullFine)/normalization,'k-', ...
    'LineWidth',2.0,'DisplayName','full corrected');
plot(profileAxis,fineRadius,real(besselFine)/normalization,'--', ...
    'LineWidth',1.8,'DisplayName','Bessel projection');
plot(profileAxis,fineRadius,real(correctionFine)/normalization,':', ...
    'LineWidth',2.0,'DisplayName','correction c_n');
yline(profileAxis,0,':','HandleVisibility','off');
xlabel(profileAxis,'r/R');
ylabel(profileAxis,'normalized surface displacement at \theta=0');
title(profileAxis,'radial decomposition at the displayed phase');
legend(profileAxis,'Location','best');
grid(profileAxis,'on');
xlim(profileAxis,[0,1]);

phaseAxis = nexttile(layout,6);
phaseFractions = analysis.temporalPhaseSamples/(2*pi);
plot(phaseAxis,phaseFractions, ...
    analysis.timeResolvedRelativeCorrectionL2,'LineWidth',1.8);
hold(phaseAxis,'on');
yline(phaseAxis,analysis.relativeCorrectionL2,'--', ...
    'DisplayName','harmonic aggregate');
xline(phaseAxis,phaseFractions(phaseIndex),':', ...
    'DisplayName','displayed phase');
xlabel(phaseAxis,'forcing phase, \tau/(2\pi)');
ylabel(phaseAxis,'||c||/||\zeta||');
title(phaseAxis,sprintf('correction norm: %.1f%%; energy: %.1f%%', ...
    100*analysis.relativeCorrectionL2, ...
    100*analysis.correctionEnergyFraction));
grid(phaseAxis,'on');
xlim(phaseAxis,[0,1]);

label = setting(userSettings,'label',analysis.label);
sgtitle(layout,sprintf(['%s: m=%d, \\beta^*=%.5g, ', ...
    'displayed at \\tau/(2\\pi)=%.3f'],label,analysis.m, ...
    analysis.betaStar,phaseFractions(phaseIndex)), ...
    'Interpreter','tex');

summary = struct();
summary.phaseIndex = phaseIndex;
summary.forcingPhaseFraction = phaseFractions(phaseIndex);
summary.normalization = normalization;
summary.relativeCorrectionL2 = analysis.relativeCorrectionL2;
summary.correctionEnergyFraction = analysis.correctionEnergyFraction;
summary.maximumTimeResolvedRelativeCorrectionL2 = ...
    analysis.maximumTimeResolvedRelativeCorrectionL2;
summary.fullRadialAtSelectedPhase = fullRadial;
summary.besselRadialAtSelectedPhase = besselRadial;
summary.correctionRadialAtSelectedPhase = correctionRadial;
summary.fineRadius = fineRadius;
summary.fullFine = fullFine;
summary.besselFine = besselFine;
summary.correctionFine = correctionFine;
end

function rotation = orientation_rotation(besselRadial,fullRadial)
[besselPeak,position] = max(abs(besselRadial));
if besselPeak > 64*eps
    reference = besselRadial(position);
else
    [~,position] = max(abs(fullRadial));
    reference = fullRadial(position);
end
if abs(reference) > 64*eps
    rotation = exp(-1i*angle(reference));
else
    rotation = 1;
end
end

function interpolated = interpolate_complex(nodes,values,query)
interpolated = interp1(nodes,real(values),query,'pchip') + ...
    1i*interp1(nodes,imag(values),query,'pchip');
end

function plot_planform(axisHandle,x,y,value,titleText)
surface(axisHandle,x,y,zeros(size(value)),value, ...
    'EdgeColor','none','FaceColor','interp');
view(axisHandle,2);
axis(axisHandle,'equal','tight');
clim(axisHandle,[-1,1]);
xlabel(axisHandle,'x/R');
ylabel(axisHandle,'y/R');
title(axisHandle,titleText);
end

function map = blue_red_map(numberOfColors)
anchors = [0.12,0.27,0.74; 0.28,0.63,0.91; 0.94,0.97,0.98; ...
    0.96,0.57,0.42; 0.70,0.05,0.12];
anchorPosition = linspace(0,1,size(anchors,1));
query = linspace(0,1,numberOfColors);
map = interp1(anchorPosition,anchors,query,'pchip');
map = min(max(map,0),1);
end

function value = setting(settings,name,defaultValue)
if isfield(settings,name) && ~isempty(settings.(name))
    value = settings.(name);
else
    value = defaultValue;
end
end
