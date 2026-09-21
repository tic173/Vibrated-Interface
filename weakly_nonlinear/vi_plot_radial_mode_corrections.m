function vi_plot_radial_mode_corrections(analyses)
%VI_PLOT_RADIAL_MODE_CORRECTIONS Display Bessel and nonseparable content.

if ~iscell(analyses)
    analyses = {analyses};
end
analyses = analyses(:);
if isempty(analyses)
    return;
end

figure('Name','Floquet radial mode corrections');
layout = tiledlayout(numel(analyses),2, ...
    'TileSpacing','compact','Padding','compact');
for modeIndex = 1:numel(analyses)
    value = analyses{modeIndex};
    dominant = value.dominantHarmonicPosition;
    phaseRotation = profile_phase(value,dominant);

    nexttile;
    hold on;
    plot(value.rOverR,real(phaseRotation* ...
        value.fullCoefficients(:,dominant)),'k-', ...
        'LineWidth',1.5,'DisplayName','full \zeta_n');
    plot(value.rOverR,real(phaseRotation* ...
        value.separableCoefficients(:,dominant)),'--', ...
        'LineWidth',1.35,'DisplayName','b_n J_m');
    plot(value.rOverR,real(phaseRotation* ...
        value.correctionCoefficients(:,dominant)),':', ...
        'LineWidth',1.5,'DisplayName','c_n');
    yline(0,':','HandleVisibility','off');
    grid on;
    xlabel('r/R');
    ylabel('phase-aligned \zeta_n');
    title(sprintf('%s, n+s=%g',value.label, ...
        value.frequencies(dominant)),'Interpreter','none');
    legend('Location','best');

    nexttile;
    imagesc(value.temporalPhaseSamples/(2*pi),value.rOverR, ...
        real(phaseRotation*value.correctionCarrierSamples));
    axis xy;
    colorbar;
    xlabel('forcing phase, \tau/(2\pi)');
    ylabel('r/R');
    title(sprintf('c_n contribution: L2 %.3g, included %s', ...
        value.relativeCorrectionL2,yes_no( ...
        value.includeInSurfacePattern)));
end
sgtitle(layout,['Radial nonseparability: ', ...
    '\zeta_n=b_nJ_m+c_n']);
end

function phaseRotation = profile_phase(value,dominant)
coefficient = value.besselTemporalCoefficients(dominant);
if abs(coefficient) > 64*eps
    phaseRotation = exp(-1i*angle(coefficient));
    return;
end
[~,peakIndex] = max(abs(value.fullCoefficients(:,dominant)));
peak = value.fullCoefficients(peakIndex,dominant);
if abs(peak) > 64*eps
    phaseRotation = exp(-1i*angle(peak));
else
    phaseRotation = 1;
end
end

function value = yes_no(flag)
if flag
    value = 'yes';
else
    value = 'no';
end
end
