function result = vi_linear_floquet(cfg)
%VI_LINEAR_FLOQUET Viscous two-fluid linear growth and interface dynamics.
% cfg=vi_linear_defaults('cylindrical3d'); result=vi_linear_floquet(cfg);
% All user inputs are SI, except explicitly named dimensionless options.
% Requires Optimization Toolbox (fsolve). See linear_floquet/README.md.
rootFolder=fileparts(fileparts(mfilename('fullpath')));
if exist(fullfile(rootFolder,'vi_reduced_cylinder_coefficients.m'),'file')
    addpath(rootFolder);
end
if exist('fsolve','file')~=2
    error('vi_linear:Toolbox','The root search requires Optimization Toolbox (fsolve).');
end
[cfg,p,modes]=vi_linear_inputs(cfg);
result=struct('config',cfg,'parameters',p,'modes',struct([]));
t=cfg.sampling.timeSeconds; ts=t/p.timeScale;
for i=1:numel(modes)
    md=modes(i);
    fprintf('Linear Floquet: %s, kh=%.6g\n',md.label,md.kStar);
    root=vi_linear_root(p,md.kStar,cfg.numerics,cfg.numerics.harmonics);
    change=NaN;
    if cfg.numerics.checkConvergence
        refined=vi_linear_root(p,md.kStar,cfg.numerics, ...
            cfg.numerics.harmonics+cfg.numerics.harmonicIncrement, ...
            [root.candidates.exponent]);
        % A conjugate pair describes the same real growth/oscillation rate.
        delta=[real(refined.exponent)-real(root.exponent), ...
            abs(imag(refined.exponent))-abs(imag(root.exponent))];
        change=norm(delta)/max([1,abs(root.exponent),p.omegaStar]);
        if change>cfg.numerics.convergenceTolerance
            error('vi_linear:HarmonicConvergence', ...
                ['Mode %s changed by %.3g between N=%d and N=%d. Increase ', ...
                 'harmonics or broaden initialGuesses before using the result.'], ...
                md.label,change,root.harmonics,refined.harmonics);
        end
        root=refined;
    end
    % Normalize the periodic factor at the actual initial forcing phase.
    % p(0)=1, so amplitude_m is the initial complex modal displacement.
    z=root.zeta; p0=sum(z);
    if abs(p0)<1e-8*norm(z)
        error('vi_linear:InitialNode', ...
            'Mode %s has p(0) near zero; change forcing.phaseRad to normalize at a non-node.',md.label);
    end
    z=z/p0;
    periodic=z.'*exp(1i*root.indices*(p.omegaStar*ts));
    amplitude=md.amplitude_m*exp(1i*md.phaseRad);
    temporal=amplitude*exp(root.exponent*ts).*periodic;
    if any(~isfinite(temporal))
        error('vi_linear:Overflow','Shorten sampling.timeSeconds: the exponential response overflowed.');
    end
    mode=struct('input',md,'exponentStar',root.exponent, ...
        'exponentPerSecond',root.exponent/p.timeScale, ...
        'growthRatePerSecond',real(root.exponent)/p.timeScale, ...
        'multiplier',exp(root.exponent*2*pi/p.omegaStar), ...
        'harmonicIndices',root.indices,'harmonicCoefficients',z, ...
        'periodicFactor',periodic,'temporalDisplacement_m',temporal, ...
        'harmonicConvergenceError',change,'rootSearch',root);
    mode.quasiFrequencyHz=imag(mode.exponentPerSecond)/(2*pi);
    if i==1, result.modes=mode; else, result.modes(i)=mode; end
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
