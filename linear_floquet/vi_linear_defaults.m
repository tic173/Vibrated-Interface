function cfg = vi_linear_defaults(geometry)
%VI_LINEAR_DEFAULTS Editable SI inputs for linear Floquet analysis.
% geometry: 'cartesian2d', 'cartesian3d', or 'cylindrical3d'.
if nargin == 0, geometry = 'cylindrical3d'; end
geometry = validatestring(geometry,{'cartesian2d','cartesian3d','cylindrical3d'});
cfg.fluids = struct('rhoLower',997,'rhoUpper',1.20, ...
    'muLower',1e-3,'muUpper',1.81e-5,'surfaceTension',0.072);
cfg.geometry = struct('type',geometry,'lowerDepth',0.022,'upperDepth',0.022, ...
    'Lx',0.07,'Ly',0.07,'radius',0.035,'horizontalBoundary','periodic');
cfg.gravity = struct('magnitude',9.81,'sign',1);
% Acceleration convention: g_eff(t) = g*sign + a*cos(2*pi*f*t+phase).
% If kind='displacement_m', amplitude is the amplitude of a/omega^2;
% phase still refers to acceleration (base displacement has phase+pi).
cfg.forcing = struct('amplitude',0.3,'kind','acceleration_g', ...
    'frequencyHz',15,'phaseRad',0);
cfg.numerics = struct('harmonics',10,'checkConvergence',true,'harmonicIncrement',4, ...
    'rootTolerance',1e-8,'convergenceTolerance',1e-4, ...
    'edgeTolerance',1e-6,'initialGuesses',[],'maxIterations',100);
cfg.sampling = struct('timeSeconds',linspace(0,10/cfg.forcing.frequencyHz,401), ...
    'Nx',81,'Ny',61,'Nr',61,'Ntheta',97);
cfg.output = struct('plot',true,'save',false,'reconstruct',true, ...
    'directory',fullfile(fileparts(mfilename('fullpath')),'results'));
mode = struct('label','mode 1','amplitude_m',1e-5,'phaseRad',0);
switch geometry
    case 'cartesian2d'
        mode.kx = 2*pi/cfg.geometry.Lx;
    case 'cartesian3d'
        mode.kx = 2*pi/cfg.geometry.Lx; mode.ky = 2*pi/cfg.geometry.Ly;
    case 'cylindrical3d'
        cfg.geometry.horizontalBoundary = 'bessel-neumann';
        mode.m = 0; mode.radialIndex = 1;
end
cfg.modes = mode;
end
