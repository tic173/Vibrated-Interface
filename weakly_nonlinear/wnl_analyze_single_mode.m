function result = wnl_analyze_single_mode(model, spec, userOpts)
%WNL_ANALYZE_SINGLE_MODE Direct, adjoint, detuning, and self coefficient.

if nargin < 3
    userOpts = struct();
end
opts = wnl_options(userOpts);
analysisWallClock = tic;
modeWallClock = tic;
mode = wnl_compute_mode(model, spec, opts);
modeSeconds = toc(modeWallClock);
modeBar = wnl_conjugate_mode(model, mode, opts);
coordinateType = wnl_mode_coordinate_type( ...
    mode,opts.realModeEigenvalueTolerance, ...
    opts.realModeMinimumConjugateOverlap);
realCoordinateProjection = struct('applied',false);
if strcmp(coordinateType,'real')
    [mode,realCoordinateProjection] = ...
        wnl_enforce_real_mode(model,mode,modeBar,opts);
    modeBar = wnl_conjugate_mode(model,mode,opts);
    modeBar.coordinateType = 'real';
else
    mode.coordinateType = 'complex';
    modeBar.coordinateType = 'complex';
end
wnl_assert_mode_converged(mode, opts);
wnl_assert_mode_converged(modeBar, opts);
% Establish that the full mode still represents the reduced operating-point
% branch before applying the stricter coefficient-level residual gate. A
% large exponent shift is a discretization/model-consistency failure, not a
% reason to spend more iterations on the same grid.
referenceConsistency = ...
    wnl_assert_reference_eigenvalue_consistent(mode,opts);
coefficientOpts = opts;
coefficientOpts.modeResidualTolerance = ...
    opts.coefficientModeResidualTolerance;
wnl_assert_mode_converged(mode, coefficientOpts);
% The physical conjugate direct field enters the nonlinear actions and must
% meet the coefficient gate. Its mapped adjoint does not enter a
% nonresonant self-cubic projection; it has already passed the ordinary
% linear gate above. If a quadratic block matches the conjugate mode, the
% resonance test invalidates the nonresonant cubic scaling independently.
coefficientBarOpts = coefficientOpts;
coefficientBarOpts.checkAdjointModeResidual = false;
wnl_assert_mode_converged(modeBar, coefficientBarOpts);
neutralModes = wnl_unique_modes({mode, modeBar});
selfWallClock = tic;
self = wnl_self_coefficient(model, mode, modeBar, ...
    neutralModes, opts);
selfCoefficientSeconds = toc(selfWallClock);

mu = NaN;
muReality = struct();
detuningForcing = [];
if isfield(model, 'detuning') && isa(model.detuning, 'function_handle') ...
        && ~isempty(opts.detuning)
    detuningForcing = model.detuning(mode.field, ...
        opts.detuning, mode.spec);
    rawMu = (mode.left' * detuningForcing(:)) / ...
        mode.normalization;
    [mu,muReality] = wnl_project_modal_coefficient(rawMu,mode,opts);
end

result = struct();
result.mode = mode;
result.modeBar = modeBar;
result.linearCoefficient = wnl_spec_lambda(mode.spec);
if abs(result.linearCoefficient) > 1.0e-13
    result.referenceType = 'operating-point Floquet reduction';
else
    result.referenceType = 'neutral-point reduction';
end
result.mu = mu;
result.g = self.g;
result.gPhysicalPeak = wnl_physical_cubic_coefficients( ...
    self.g,{mode});
result.coordinateType = coordinateType;
result.amplitudeScaleToPeakZetaOverH = ...
    wnl_mode_amplitude_scale(mode);
result.realCoordinateProjection = realCoordinateProjection;
result.detuningCoefficientReality = muReality;
result.coefficientConventionVersion = ...
    'V2-real-self-conjugate-physical-peak';
result.self = self;
result.detuningForcing = detuningForcing;
result.validCubicScaling = self.validCubicScaling;
result.quadraticResonance = self.quadraticResonance;
result.forcedSolvesValid = self.forcedSolvesValid;
result.forcedSolvesExploratoryUsable = ...
    self.forcedSolvesExploratoryUsable;
result.referenceEigenvalueConsistency = referenceConsistency;
result.timing = struct('modeSeconds',modeSeconds, ...
    'selfCoefficientSeconds',selfCoefficientSeconds, ...
    'totalSeconds',toc(analysisWallClock));
if opts.verbose
    fprintf(['WNL single-mode timing: mode %.3f s, self coefficient ', ...
        '%.3f s, total %.3f s\n'],modeSeconds,selfCoefficientSeconds, ...
        result.timing.totalSeconds);
end
if strcmp(result.referenceType, 'operating-point Floquet reduction')
    result.amplitudeEquation = ...
        'dA/dt = lambda*A + g*abs(A)^2*A';
else
    result.amplitudeEquation = ...
        'dA/dT = mu*A + g*abs(A)^2*A';
end
end
