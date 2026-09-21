function result = vi_linear_dynamics(result)
%VI_LINEAR_DYNAMICS Reconstruct previously computed Floquet modes.
% Input comes from vi_linear_growth. No root search is repeated.
cfg=result.config; p=result.parameters; modes=[result.modes.input];
t=cfg.sampling.timeSeconds; ts=t/p.timeScale;
for i=1:numel(result.modes)
    mode=result.modes(i); md=mode.input;
    % Normalize the periodic factor at the actual initial forcing phase.
    % p(0)=1, so amplitude_m is the initial complex modal displacement.
    z=mode.harmonicCoefficients; p0=sum(z);
    if abs(p0)<1e-8*norm(z)
        error('vi_linear:InitialNode', ...
            'Mode %s has p(0) near zero; change forcing.phaseRad to normalize at a non-node.',md.label);
    end
    z=z/p0;
    periodic=z.'*exp(1i*mode.harmonicIndices*(p.omegaStar*ts));
    amplitude=md.amplitude_m*exp(1i*md.phaseRad);
    temporal=amplitude*exp(mode.exponentStar*ts).*periodic;
    if any(~isfinite(temporal))
        error('vi_linear:Overflow','Shorten sampling.timeSeconds: the exponential response overflowed.');
    end
    result.modes(i).harmonicCoefficients=z;
    result.modes(i).periodicFactor=periodic;
    result.modes(i).temporalDisplacement_m=temporal;
end
result.timeSeconds=t;
result.modalDisplacement_m=vertcat(result.modes.temporalDisplacement_m);
result.maximumGrowthRatePerSecond=max([result.modes.growthRatePerSecond]);
result.stabilityScope='Largest accepted root for each requested spatial mode; multistart search is not a complete spectrum.';
result.grid=[]; result.interface_m=[];
if cfg.output.reconstruct
    [grid,basis]=vi_linear_spatial(cfg,modes);
    interface=zeros([size(grid.X) numel(t)]);
    for i=1:numel(modes)
        interface=interface+real(reshape(basis{i},[size(grid.X) 1]).* ...
            reshape(result.modalDisplacement_m(i,:),1,1,[]));
    end
    result.grid=grid; result.interface_m=interface;
end
result.outputDirectory='';
if cfg.output.save, result=vi_linear_save(result); end
if cfg.output.plot, vi_linear_plot(result); end
end
