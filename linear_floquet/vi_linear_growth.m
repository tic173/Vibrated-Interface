function result = vi_linear_growth(cfg)
%VI_LINEAR_GROWTH Growth-only analysis of every cfg.modes entry.
% No time histories, spatial arrays, plots, or files are created. Returned
% harmonic vectors are unnormalized; vi_linear_dynamics normalizes p(0)=1.
rootFolder=fileparts(fileparts(mfilename('fullpath')));
if exist(fullfile(rootFolder,'vi_reduced_cylinder_coefficients.m'),'file')
    addpath(rootFolder);
end
if exist('fsolve','file')~=2
    error('vi_linear:Toolbox','The root search requires Optimization Toolbox (fsolve).');
end
[cfg,p,modes]=vi_linear_inputs(cfg);
result=struct('config',cfg,'parameters',p,'modes',struct([]));
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
    mode=struct('input',md,'exponentStar',root.exponent, ...
        'exponentPerSecond',root.exponent/p.timeScale, ...
        'growthRatePerSecond',real(root.exponent)/p.timeScale, ...
        'multiplier',exp(root.exponent*2*pi/p.omegaStar), ...
        'harmonicIndices',root.indices,'harmonicCoefficients',root.zeta, ...
        'harmonicConvergenceError',change,'rootSearch',root);
    mode.quasiFrequencyHz=imag(mode.exponentPerSecond)/(2*pi);
    if i==1, result.modes=mode; else, result.modes(i)=mode; end
end
result.maximumGrowthRatePerSecond=max([result.modes.growthRatePerSecond]);
result.stabilityScope='Largest accepted root for each requested spatial mode; multistart search is not a complete spectrum.';
end
