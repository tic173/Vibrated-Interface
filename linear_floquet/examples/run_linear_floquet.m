%% Linear Floquet stability and interfacial dynamics -- editable SI inputs
% Run from any folder with this script's full path, or from the repo root:
% run('linear_floquet/examples/run_linear_floquet.m')
linearFolder=fileparts(fileparts(mfilename('fullpath')));
addpath(linearFolder);

%% User inputs
% Choose: 'cartesian2d', 'cartesian3d', or 'cylindrical3d'.
cfg=vi_linear_defaults('cylindrical3d');
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

% Spatial mode selection. Edit the fields for the selected geometry.
% More modes can be added by copying cfg.modes(1) to cfg.modes(2), etc.
switch cfg.geometry.type
    case 'cylindrical3d'
        cfg.modes(1).m=0; cfg.modes(1).radialIndex=1;
    case 'cartesian2d'
        cfg.modes(1).kx=2*pi/cfg.geometry.Lx; % rad/m
    case 'cartesian3d'
        cfg.modes(1).kx=2*pi/cfg.geometry.Lx;
        cfg.modes(1).ky=2*pi/cfg.geometry.Ly;
end
cfg.modes(1).label='mode 1';
cfg.modes(1).amplitude_m=1e-5;       % initial complex modal displacement scale
cfg.modes(1).phaseRad=0;             % disturbance phase, distinct from forcing

cfg.numerics.harmonics=10;          % temporal indices -N:N
cfg.numerics.checkConvergence=true; % compare N against N+harmonicIncrement
cfg.sampling.timeSeconds=linspace(0,10/cfg.forcing.frequencyHz,401);
cfg.output.plot=true;
cfg.output.reconstruct=true;
cfg.output.save=true;               % CSV + MAT in a unique results/run_* folder

%% Run
result=vi_linear_floquet(cfg);
for i=1:numel(result.modes)
    fprintf('%s: growth=%+.6g 1/s, quasi-frequency=%+.6g Hz\n', ...
        result.modes(i).input.label,result.modes(i).growthRatePerSecond, ...
        result.modes(i).quasiFrequencyHz);
end
% result.interface_m(:,:,it) is the interface at result.timeSeconds(it).
% For 2D Cartesian, size(interface_m) is [1 Nx Nt].
% For 3D Cartesian, size(interface_m) is [Ny Nx Nt].
% For cylindrical, size(interface_m) is [Ntheta Nr Nt].
