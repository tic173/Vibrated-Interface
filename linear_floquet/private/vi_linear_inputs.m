function [cfg,p,modes] = vi_linear_inputs(cfg)
% Validate the complete, editable configuration before any expensive solve.
if ~isstruct(cfg) || ~isscalar(cfg)
    error('vi_linear:Input','cfg must be a scalar struct from vi_linear_defaults.');
end
required={'fluids','geometry','gravity','forcing','numerics','sampling','output','modes'};
for i=1:numel(required)
    if ~isfield(cfg,required{i}), error('vi_linear:MissingInput','Missing cfg.%s.',required{i}); end
end
cfg.geometry.type=validatestring(cfg.geometry.type, ...
    {'cartesian2d','cartesian3d','cylindrical3d'});
default=vi_linear_defaults(cfg.geometry.type);
unknown=setdiff(fieldnames(cfg),fieldnames(default));
if ~isempty(unknown), error('vi_linear:UnknownInput','Unknown cfg.%s.',unknown{1}); end
for s=required(1:end-1)
    key=s{1};
    if ~isstruct(cfg.(key)) || ~isscalar(cfg.(key))
        error('vi_linear:Input','cfg.%s must be a scalar struct.',key);
    end
    absent=setdiff(fieldnames(default.(key)),fieldnames(cfg.(key)));
    extra=setdiff(fieldnames(cfg.(key)),fieldnames(default.(key)));
    if ~isempty(absent), error('vi_linear:MissingInput','Missing cfg.%s.%s.',key,absent{1}); end
    if ~isempty(extra), error('vi_linear:UnknownInput','Unknown cfg.%s.%s.',key,extra{1}); end
end
for key={'rhoLower','rhoUpper','muLower','muUpper','surfaceTension'}
    positive(cfg.fluids.(key{1}),['fluids.' key{1}]);
end
for key={'lowerDepth','upperDepth','Lx','Ly','radius'}
    positive(cfg.geometry.(key{1}),['geometry.' key{1}]);
end
positive(cfg.gravity.magnitude,'gravity.magnitude');
if ~ismember(cfg.gravity.sign,[-1 1]) || ~isscalar(cfg.gravity.sign)
    error('vi_linear:Gravity','gravity.sign must be +1 or -1.');
end
positive(cfg.forcing.frequencyHz,'forcing.frequencyHz');
validateattributes(cfg.forcing.amplitude,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(cfg.forcing.phaseRad,{'numeric'},{'scalar','real','finite'});
cfg.forcing.kind=validatestring(cfg.forcing.kind, ...
    {'acceleration_g','acceleration_m_s2','displacement_m'});
for key={'harmonics','harmonicIncrement','maxIterations'}
    validateattributes(cfg.numerics.(key{1}),{'numeric'},{'scalar','integer','positive','finite'});
end
for key={'rootTolerance','convergenceTolerance','edgeTolerance'}
    positive(cfg.numerics.(key{1}),['numerics.' key{1}]);
end
if ~isempty(cfg.numerics.initialGuesses)
    validateattributes(cfg.numerics.initialGuesses,{'numeric'},{'vector','finite'});
end
for key={'Nx','Ny','Nr','Ntheta'}
    validateattributes(cfg.sampling.(key{1}),{'numeric'},{'scalar','integer','>=',3,'finite'});
end
validateattributes(cfg.sampling.timeSeconds,{'numeric'},{'vector','real','finite','nonnegative','nonempty'});
cfg.sampling.timeSeconds=cfg.sampling.timeSeconds(:).';
if any(diff(cfg.sampling.timeSeconds)<=0)
    error('vi_linear:Time','timeSeconds must be strictly increasing.');
end
for key={'plot','save','reconstruct'}
    validateattributes(cfg.output.(key{1}),{'logical','numeric'},{'scalar','binary'});
end
validateattributes(cfg.numerics.checkConvergence,{'logical','numeric'},{'scalar','binary'});
if ~(ischar(cfg.output.directory) || (isstring(cfg.output.directory) && isscalar(cfg.output.directory)))
    error('vi_linear:Output','output.directory must be a folder path.');
end
if isempty(cfg.output.directory), error('vi_linear:Output','output.directory cannot be empty.'); end
f=cfg.fluids; g=cfg.gravity.magnitude; h=cfg.geometry.lowerDepth;
p=struct('h',h,'g',g,'timeScale',sqrt(h/g), ...
    'At',(f.rhoLower-f.rhoUpper)/(f.rhoLower+f.rhoUpper), ...
    'eta',f.muUpper/f.muLower,'C',(f.muLower/f.rhoLower)/sqrt(g*h^3), ...
    'Bd',f.rhoLower*g*h^2/f.surfaceTension,'gSign',cfg.gravity.sign, ...
    'phase',cfg.forcing.phaseRad,'depths',[1 cfg.geometry.upperDepth/h]);
if abs(p.At)<1e-8
    error('vi_linear:DensityContrast','This density-jump formulation requires abs(At)>=1e-8.');
end
p.omega=2*pi*cfg.forcing.frequencyHz; p.omegaStar=p.omega*p.timeScale;
switch cfg.forcing.kind
    case 'acceleration_g', p.Ac=cfg.forcing.amplitude;
    case 'acceleration_m_s2', p.Ac=cfg.forcing.amplitude/g;
    case 'displacement_m', p.Ac=cfg.forcing.amplitude*p.omega^2/g;
end
p.accelerationSI=p.Ac*g; p.displacementSI=p.accelerationSI/p.omega^2;
if strcmp(cfg.geometry.type,'cylindrical3d')
    if ~strcmp(cfg.geometry.horizontalBoundary,'bessel-neumann')
        error('vi_linear:Boundary','The cylinder requires bessel-neumann; pinned/no-slip sidewalls are not this model.');
    end
    extra={'m','radialIndex'};
else
    cfg.geometry.horizontalBoundary=validatestring(cfg.geometry.horizontalBoundary,{'periodic','unbounded'});
    extra={'kx'};
    if strcmp(cfg.geometry.type,'cartesian3d'), extra={'kx','ky'}; end
end
if ~isstruct(cfg.modes) || isempty(cfg.modes) || ~isvector(cfg.modes)
    error('vi_linear:Modes','modes must be a nonempty struct vector.');
end
fields=[{'label','amplitude_m','phaseRad'},extra];
if ~isempty(setxor(fieldnames(cfg.modes),fields))
    error('vi_linear:Modes','Mode fields for %s must be: %s.',cfg.geometry.type,strjoin(fields,', '));
end
modes=cfg.modes;
for i=1:numel(modes)
    md=modes(i);
    if ~(ischar(md.label) || (isstring(md.label)&&isscalar(md.label)))
        error('vi_linear:Modes','Each mode needs a text label.');
    end
    validateattributes(md.amplitude_m,{'numeric'},{'scalar','real','finite','nonnegative'});
    validateattributes(md.phaseRad,{'numeric'},{'scalar','real','finite'});
    if strcmp(cfg.geometry.type,'cylindrical3d')
        validateattributes(md.m,{'numeric'},{'scalar','integer','nonnegative','finite'});
        validateattributes(md.radialIndex,{'numeric'},{'scalar','integer','positive','finite'});
        roots=bessel_derivative_root(md.m,md.radialIndex);
        k=roots(end)/cfg.geometry.radius;
    else
        validateattributes(md.kx,{'numeric'},{'scalar','real','finite'});
        vector=md.kx; lengths=cfg.geometry.Lx;
        if strcmp(cfg.geometry.type,'cartesian3d')
            validateattributes(md.ky,{'numeric'},{'scalar','real','finite'});
            vector=[md.kx md.ky]; lengths=[cfg.geometry.Lx cfg.geometry.Ly];
        end
        k=norm(vector);
        if strcmp(cfg.geometry.horizontalBoundary,'periodic')
            index=vector.*lengths/(2*pi);
            if any(abs(index-round(index))>1e-9*max(1,abs(index)))
                error('vi_linear:Periodicity','For periodic modes, kx*Lx and ky*Ly must be integer multiples of 2*pi.');
            end
        end
    end
    if k==0, error('vi_linear:ZeroMode','The volume-changing spatially uniform mode is excluded.'); end
    modes(i).k=k; modes(i).kStar=k*h;
end
end

function positive(value,name)
validateattributes(value,{'numeric'},{'scalar','real','positive','finite'},'',name);
end
