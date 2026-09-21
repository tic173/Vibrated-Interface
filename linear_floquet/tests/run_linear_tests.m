function run_linear_tests
%RUN_LINEAR_TESTS Physical invariants and user-facing workflow regression.
linearFolder=fileparts(fileparts(mfilename('fullpath')));
addpath(linearFolder); addpath(fileparts(linearFolder));
oldVisible=get(groot,'defaultFigureVisible');
set(groot,'defaultFigureVisible','off');
restore=onCleanup(@() set(groot,'defaultFigureVisible',oldVisible));
cfg=vi_linear_defaults('cylindrical3d');
cfg.output.plot=false; cfg.numerics.harmonics=6;
cfg.sampling.timeSeconds=(0:2)/cfg.forcing.frequencyHz;
a=vi_linear_floquet(cfg); s=a.modes.exponentStar;
assert(a.modes.rootSearch.residual<1e-8);
assert(a.modes.harmonicConvergenceError<1e-4);
assert(abs(a.modalDisplacement_m(1)-cfg.modes.amplitude_m)<1e-14);
assert(abs(a.modalDisplacement_m(2)/a.modalDisplacement_m(1)-a.modes.multiplier)<1e-8);
assert(abs(max(abs(a.interface_m(:,:,1)),[],'all')-cfg.modes.amplitude_m)<1e-12);

% Same k and fluids must give the same temporal equation in all geometries.
k=a.modes.input.k;
for type={'cartesian2d','cartesian3d'}
    c=vi_linear_defaults(type{1}); c.output.plot=false;
    c.geometry.Lx=2*pi/k; c.modes.kx=k;
    if strcmp(type{1},'cartesian3d')
        c.modes.kx=.6*k; c.modes.ky=.8*k;
        c.geometry.Lx=2*pi/c.modes.kx; c.geometry.Ly=2*pi/c.modes.ky;
    end
    c.numerics.harmonics=10; c.numerics.checkConvergence=false;
    c.numerics.initialGuesses=[s conj(s)];
    c.sampling.timeSeconds=cfg.sampling.timeSeconds;
    b=vi_linear_floquet(c);
    assert(abs(b.modes.exponentStar-s)<1e-7);
    assert(max(abs(b.modalDisplacement_m-a.modalDisplacement_m))<1e-9);
    vi_linear_plot(b); close all;
end

% Acceleration phase is a time-origin change: exponent invariant, periodic
% coefficients rotated by exp(i*n*phase), with p(0)=1 in both cases.
c=cfg; c.numerics.harmonics=10; c.numerics.checkConvergence=false;
c.numerics.initialGuesses=[s conj(s)]; c.forcing.phaseRad=.73;
b=vi_linear_floquet(c);
expected=a.modes.harmonicCoefficients.*exp(1i*a.modes.harmonicIndices*.73);
expected=expected/sum(expected);
assert(abs(b.modes.exponentStar-s)<1e-7);
assert(norm(b.modes.harmonicCoefficients-expected)<1e-6);

% The three amplitude units describe exactly the same forcing.
for kind={'acceleration_m_s2','displacement_m'}
    c=cfg; c.numerics.harmonics=10; c.numerics.checkConvergence=false;
    c.numerics.initialGuesses=[s conj(s)]; c.forcing.kind=kind{1};
    c.forcing.amplitude=a.parameters.accelerationSI;
    if strcmp(kind{1},'displacement_m'), c.forcing.amplitude=a.parameters.displacementSI; end
    b=vi_linear_floquet(c); assert(abs(b.modes.exponentStar-s)<1e-7);
end

% Existing cylindrical public solver stays compatible at the same root.
p=a.parameters;
[legacy,z,~,~,~,~,d]=faradayFloquet_RT_GenEIG_cylindrical( ...
    p.Ac,p.omegaStar,cfg.geometry.radius/p.h,0,1,p.C,p.Bd,p.At,p.eta,10,'H',1,s,0);
assert(abs(legacy-s)<1e-7 && d.relativeSingularResidual<1e-8);
z=z/sum(z); assert(norm(z-a.modes.harmonicCoefficients)<1e-6);

% Low-viscosity, unforced, unequal-depth limit approaches the known
% capillary-gravity dispersion relation with finite-depth coth factors.
c=cfg; c.geometry.upperDepth=.013; c.fluids.muLower=1e-8;
c.fluids.muUpper=1.81e-10; c.forcing.amplitude=0;
f=c.fluids; h1=c.geometry.lowerDepth; h2=c.geometry.upperDepth;
w=sqrt(((f.rhoLower-f.rhoUpper)*c.gravity.magnitude*k+f.surfaceTension*k^3)/ ...
    (f.rhoLower*coth(k*h1)+f.rhoUpper*coth(k*h2)));
c.numerics.initialGuesses=[-.001+1i*w*sqrt(h1/c.gravity.magnitude), ...
    -.001-1i*w*sqrt(h1/c.gravity.magnitude)];
b=vi_linear_floquet(c);
assert(abs(abs(imag(b.modes.exponentPerSecond))-w)/w<.01);
assert(b.modes.growthRatePerSecond<0);

% Long-wave RT mode grows without forcing. This is not an H/SH onset test.
c=cfg; c.gravity.sign=-1; c.forcing.amplitude=0;
b=vi_linear_floquet(c); assert(b.modes.growthRatePerSecond>0);
rt=b.modes.exponentStar;
c.gravity.sign=1;
c.fluids.rhoLower=cfg.fluids.rhoUpper; c.fluids.rhoUpper=cfg.fluids.rhoLower;
c.fluids.muLower=cfg.fluids.muUpper; c.fluids.muUpper=cfg.fluids.muLower;
b=vi_linear_floquet(c); assert(abs(b.modes.exponentStar-rt)<1e-6);


% Multiple modes, disturbance phase, spatial arrays, plotting and saved I/O.
c=cfg; c.numerics.checkConvergence=false;
c.modes(2)=c.modes(1); c.modes(2).m=1; c.modes(2).label='m=1';
c.modes(2).phaseRad=.4; c.output.plot=true; c.output.save=true;
tempFolder=tempname; mkdir(tempFolder);
cleanup=onCleanup(@() rmdir(tempFolder,'s'));
c.output.directory=tempFolder;
b=vi_linear_floquet(c);
assert(size(b.modalDisplacement_m,1)==2);
assert(abs(b.modalDisplacement_m(2,1)-c.modes(2).amplitude_m*exp(.4i))<1e-12);
loaded=load(fullfile(b.outputDirectory,'linear_floquet.mat'),'result');
assert(isequaln(loaded.result.interface_m,b.interface_m));
assert(height(readtable(fullfile(b.outputDirectory,'growth_rates.csv')))==2);
close all;

% Sweep computes all growth rates but reconstructs only the sampled maximum.
c=vi_linear_defaults('cartesian2d'); c.geometry.horizontalBoundary='unbounded';
c.numerics.harmonics=6; c.numerics.checkConvergence=false;
c.output.plot=true; c.output.save=true; c.output.directory=tempFolder;
c.sampling.timeSeconds=(0:2)/c.forcing.frequencyHz;
for j=1:3
    c.modes(j)=c.modes(1); c.modes(j).kx=(2*j-1)/c.geometry.lowerDepth;
    c.modes(j).label=sprintf('sample %d',j);
end
growth=vi_linear_growth(c);
assert(~isfield(growth,'interface_m') && ~isfield(growth.modes,'temporalDisplacement_m'));
b=vi_linear_most_unstable(c);
[expected,index]=max([growth.modes.growthRatePerSecond]);
assert(b.sweep.selectedIndex==index && numel(b.modes)==1);
assert(abs(b.modes.growthRatePerSecond-expected)<1e-7);
assert(size(b.modalDisplacement_m,1)==1 && numel(b.config.modes)==1);
assert(height(readtable(fullfile(b.outputDirectory,'growth_rates.csv')))==3);
assert(~isfield(b.sweep.modes,'temporalDisplacement_m'));
assert(height(readtable(fullfile(b.outputDirectory,'modal_dynamics.csv')))==numel(c.sampling.timeSeconds));
close all;

% Fail clearly for unsupported/ambiguous physical inputs.
c=cfg; c.fluids.rhoUpper=c.fluids.rhoLower; must_fail(c,'vi_linear:DensityContrast');
c=cfg; c.geometry.horizontalBoundary='pinned'; must_fail(c,'vi_linear:Boundary');
c=vi_linear_defaults('cartesian2d'); c.modes.kx=0; must_fail(c,'vi_linear:ZeroMode');
c=vi_linear_defaults('cartesian2d'); c.modes.kx=1; must_fail(c,'vi_linear:Periodicity');
c=cfg; c.forcing.amplitdue=1; must_fail(c,'vi_linear:UnknownInput');
fprintf('PASS: geometry equivalence, phase, units, legacy compatibility, inviscid limit, RT, multimode I/O, growth-only sweep selection, and input validation.\n');
end

function must_fail(c,identifier)
try
    vi_linear_floquet(c);
catch err
    assert(strcmp(err.identifier,identifier),'Unexpected failure: %s',err.message);
    return
end
error('Expected error %s was not raised.',identifier);
end
