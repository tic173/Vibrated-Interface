function report = benchmark_radial_mode_correction( ...
        plotResult,viscousParameter)
%BENCHMARK_RADIAL_MODE_CORRECTION Small full-cylinder nonseparability test.

if nargin < 1
    plotResult = false;
end
p = benchmark_parameters();
if nargin >= 2 && ~isempty(viscousParameter)
    validateattributes(viscousParameter,{'numeric'}, ...
        {'scalar','real','positive','finite'});
    p.C = viscousParameter;
end
models = {'legacyFreeSlide','stressFree'};
report = repmat(struct(),numel(models),1);
for modelIndex = 1:numel(models)
    p.boundary.sidewallTangentialCondition = models{modelIndex};
    [operators,metadata] = cylinder_wnl_operators(p);
    spec = wnl_spec(p.modes.m,0,0,metadata.ndof, ...
        ['benchmark_',models{modelIndex}]);
    spec.betaStar = p.modes.betaStar;
    B = operators.B0(spec);
    L = operators.Lhat(spec,0);
    solveClock = tic;
    [vectors,values] = eig(full(L),full(B),'qz');
    eigenSeconds = toc(solveClock);
    [vector,lambda,overlap] = select_bessel_branch( ...
        vectors,diag(values),metadata,p);
    spec.lambda = lambda;
    field = metadata.normalizeDirect(wnl_make_field(spec,vector));
    mode.spec = spec;
    mode.field = field;
    settings.radialCorrectionPolicy = 'auto';
    settings.radialCorrectionRelativeL2Threshold = 1.0e-2;
    settings.radialCorrectionTemporalSamples = 64;
    correctionClock = tic;
    correction = vi_radial_mode_correction( ...
        mode,metadata,p,settings);
    correctionSeconds = toc(correctionClock);
    residual = norm((lambda*B-L)*wnl_field_vector(field)) / ...
        max(norm(L*wnl_field_vector(field)),eps);
    report(modelIndex).sidewallTangentialCondition = models{modelIndex};
    report(modelIndex).linearOperatorVersion = ...
        metadata.linearOperatorVersion;
    report(modelIndex).lambda = lambda;
    report(modelIndex).eigenResidual = residual;
    report(modelIndex).besselOverlap = overlap;
    report(modelIndex).relativeCorrectionL2 = ...
        correction.relativeCorrectionL2;
    report(modelIndex).peakCorrectionFraction = ...
        correction.peakCorrectionFraction;
    report(modelIndex).includeInSurfacePattern = ...
        correction.includeInSurfacePattern;
    report(modelIndex).eigenSeconds = eigenSeconds;
    report(modelIndex).correctionSeconds = correctionSeconds;
    report(modelIndex).correction = correction;
    fprintf(['Radial-correction benchmark %-16s: lambda ', ...
        '= %.6g%+.6gi, eigen residual %.3e, Bessel overlap %.6f, ', ...
        '||c||/||zeta|| %.6e, auto included=%d, ', ...
        'eigen %.3f s, projection %.3g s\n'], ...
        models{modelIndex},real(lambda),imag(lambda),residual,overlap, ...
        correction.relativeCorrectionL2, ...
        correction.includeInSurfacePattern,eigenSeconds,correctionSeconds);
end
if plotResult
    vi_plot_radial_mode_corrections({report.correction});
end
end

function p = benchmark_parameters()
p.omegaStar = 8.0;
p.R0 = 35/22;
p.C = 9.82e-5;
p.Bd = 65.78;
p.At = 0.9976;
p.eta = 1.81e-2;
p.g_sgn = -1;
p.aCritical = 0;
p.phase = 0;
p.numerics.Nr = 7;
p.numerics.NzUpper = 7;
p.numerics.NzLower = 7;
p.numerics.verticalGrid.type = 'single';
p.modes.m = 2;
p.modes.radialIndex = 1;
roots = bessel_derivative_root( ...
    p.modes.m,p.modes.radialIndex);
p.modes.betaStar = roots(p.modes.radialIndex)/p.R0;
p.numerics.radialGrid.type = 'besselEnriched';
p.numerics.radialGrid.maximumProductOrder = 2;
p.numerics.radialGrid.maximumConditionNumber = 1.0e10;
p.numerics.radialGrid.independenceTolerance = 1.0e-10;
p.numerics.radialGrid.fallbackToChebyshev = false;
p.numerics.Ntheta = 12;
p.numerics.quadraticStep = 2.0e-4;
p.numerics.cubicStep = 2.0e-3;
p.boundary.contactLine = 'free';
end

function [selectedVector,selectedLambda,selectedOverlap] = ...
        select_bessel_branch(vectors,values,metadata,p)
r = metadata.discretization.r(:);
quadrature = metadata.discretization.radial.quadratureWeights(:).*r;
shape = besselj(abs(p.modes.m),p.modes.betaStar*r);
shapeNorm = sqrt(real(shape'*(quadrature.*shape)));
zetaRows = metadata.layout.zeta;
bestScore = -Inf;
selectedPosition = NaN;
selectedMagnitude = Inf;
for position = 1:numel(values)
    lambda = values(position);
    if ~isfinite(real(lambda)) || ~isfinite(imag(lambda)) || ...
            abs(lambda) > 1.0e3
        continue;
    end
    zeta = vectors(zetaRows,position);
    zetaNorm = sqrt(max(real(zeta'*(quadrature.*zeta)),0));
    if zetaNorm <= 1.0e-10*max(norm(vectors(:,position)),1)
        continue;
    end
    overlap = abs(shape'*(quadrature.*zeta)) / ...
        max(shapeNorm*zetaNorm,eps);
    sameShape = abs(overlap-bestScore) <= 1.0e-8;
    if overlap > bestScore+1.0e-8 || ...
            (sameShape && abs(lambda) < selectedMagnitude)
        bestScore = overlap;
        selectedPosition = position;
        selectedMagnitude = abs(lambda);
    end
end
if ~isfinite(selectedPosition)
    error('benchmark_radial_mode_correction:NoFiniteMode', ...
        'No finite interface eigenmode was found.');
end
selectedVector = vectors(:,selectedPosition);
selectedLambda = values(selectedPosition);
selectedOverlap = bestScore;
end
