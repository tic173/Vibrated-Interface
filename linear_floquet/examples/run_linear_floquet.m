%% Linear Floquet wavenumber scan, then dynamics of the fastest mode
% Run from any folder with this script's full path, or from the repo root:
% run('linear_floquet/examples/run_linear_floquet.m')
linearFolder=fileparts(fileparts(mfilename('fullpath')));
addpath(linearFolder);

%% User inputs
% Choose: 'cartesian2d', 'cartesian3d', or 'cylindrical3d'.
cfg=vi_linear_defaults('cartesian2d');
cfg.fluids.rhoLower=997;             % kg/m^3
cfg.fluids.rhoUpper=1.20;            % kg/m^3
cfg.fluids.muLower=1e-3;             % dynamic viscosity, Pa s
cfg.fluids.muUpper=1.81e-5;          % dynamic viscosity, Pa s
cfg.fluids.surfaceTension=0.072;     % N/m
cfg.geometry.lowerDepth=22e-3;      % m, reference length h
cfg.geometry.upperDepth=22e-3;      % m (may differ from lowerDepth)
cfg.geometry.radius=35e-3;          % cylinder radius, m
cfg.geometry.Lx=70e-3;              % Cartesian period / plotting window, m
cfg.geometry.Ly=70e-3;              % 3D Cartesian period / plotting window, m
cfg.gravity.magnitude=9.81;         % positive gravity scale, m/s^2
cfg.gravity.sign=1;                 % +1 usual gravity; -1 reversed gravity

cfg.forcing.kind='acceleration_g';  % or 'acceleration_m_s2', 'displacement_m'
cfg.forcing.amplitude=0.3;          % units specified by kind
cfg.forcing.frequencyHz=15;         % Hz, not rad/s
cfg.forcing.phaseRad=0;             % phase of ACCELERATION cosine at t=0

% Scan controls: effective nondimensional wavenumber k_eff*h.
% Cartesian modes use continuous Fourier wavenumbers ('unbounded'); Lx/Ly
% define the plotted window. In 3D, directionRad sets atan2(ky,kx).
kStarRange=linspace(0.2,20,81);
directionRad=0;
% The cylinder must instead use discrete Bessel modes at the fixed radius.
azimuthalOrders=0;
radialIndices=1:20;
initialAmplitude_m=1e-5;
initialDisturbancePhaseRad=0;

seed=cfg.modes(1);
seed.amplitude_m=initialAmplitude_m;
seed.phaseRad=initialDisturbancePhaseRad;
if strcmp(cfg.geometry.type,'cylindrical3d')
    cfg.modes=repmat(seed,1,numel(azimuthalOrders)*numel(radialIndices));
    index=0;
    for m=azimuthalOrders
        for j=radialIndices
            index=index+1;
            cfg.modes(index).m=m; cfg.modes(index).radialIndex=j;
            cfg.modes(index).label=sprintf('m=%d, j=%d',m,j);
        end
    end
else
    validateattributes(kStarRange,{'numeric'},{'vector','positive','real','finite','nonempty'});
    validateattributes(directionRad,{'numeric'},{'scalar','real','finite'});
    cfg.geometry.horizontalBoundary='unbounded';
    cfg.modes=repmat(seed,1,numel(kStarRange));
    for index=1:numel(kStarRange)
        k=kStarRange(index)/cfg.geometry.lowerDepth;
        cfg.modes(index).kx=k;
        if strcmp(cfg.geometry.type,'cartesian3d')
            cfg.modes(index).kx=k*cos(directionRad);
            cfg.modes(index).ky=k*sin(directionRad);
        end
        cfg.modes(index).label=sprintf('kh=%.4g',kStarRange(index));
    end
end

cfg.numerics.harmonics=10;          % temporal indices -N:N
cfg.numerics.checkConvergence=true; % compare N against N+harmonicIncrement
cfg.sampling.timeSeconds=linspace(0,10/cfg.forcing.frequencyHz,401);
cfg.output.plot=true;
cfg.output.reconstruct=true;
cfg.output.save=true;               % CSV + MAT in a unique results/run_* folder

%% Run
result=vi_linear_most_unstable(cfg);
growthRates=result.sweep.table; % Every sampled wavenumber and its growth rate
for i=1:numel(result.modes)
    fprintf('%s: growth=%+.6g 1/s, quasi-frequency=%+.6g Hz\n', ...
        result.modes(i).input.label,result.modes(i).growthRatePerSecond, ...
        result.modes(i).quasiFrequencyHz);
end
% result.modes contains ONLY the selected mode; result.sweep stores the scan.
% result.interface_m(:,:,it) is the interface at result.timeSeconds(it).
% For 2D Cartesian, size(interface_m) is [1 Nx Nt].
% For 3D Cartesian, size(interface_m) is [Ny Nx Nt].
% For cylindrical, size(interface_m) is [Ntheta Nr Nt].
