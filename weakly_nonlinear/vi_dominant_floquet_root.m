function result = vi_dominant_floquet_root(problem, userOptions)
%VI_DOMINANT_FLOQUET_ROOT Multistart search for one dominant Floquet root.
%
% Required problem fields:
%   Ac, omegaStar, R0, m, radialIndex, C, Bd, At, eta, N,
%   modeType, g_sgn, targetS; phase is optional and defaults to zero
%
% Options:
%   realGuesses          real parts of the initial guesses. Empty creates
%                        a symmetric set scaled by omegaStar.
%   imaginaryOffsets     offsets from targetS in units of omegaStar.
%   residualTolerance    accepted relative smallest singular value.
%   duplicateTolerance   relative tolerance for duplicate Floquet roots.
%   verbose              print accepted roots.
%   suppressSolverOutput hide rejected fsolve trial messages.
%
% The largest real part among distinct accepted roots is returned. The
% imaginary part is canonicalized only for root comparison; the raw root
% and its matching harmonic vector are retained for reconstruction.

required = {'Ac', 'omegaStar', 'R0', 'm', 'radialIndex', ...
    'C', 'Bd', 'At', 'eta', 'N', 'modeType', 'g_sgn', 'targetS'};
for j = 1:numel(required)
    if ~isfield(problem, required{j})
        error('vi_dominant_floquet_root:MissingProblemField', ...
            'problem.%s is required.', required{j});
    end
end

options = struct();
options.realGuesses = [];
options.imaginaryOffsets = [-0.05, 0, 0.05];
options.residualTolerance = 1.0e-7;
options.duplicateTolerance = 1.0e-6;
options.verbose = true;
options.suppressSolverOutput = false;
if nargin >= 2 && ~isempty(userOptions)
    names = fieldnames(userOptions);
    for j = 1:numel(names)
        options.(names{j}) = userOptions.(names{j});
    end
end

if ~isfield(problem, 'phase')
    problem.phase = 0;
end
validateattributes(problem.phase, {'numeric'}, ...
    {'scalar', 'real', 'finite'});

if isempty(options.realGuesses)
    realSpan = max(2.0, problem.omegaStar);
    options.realGuesses = realSpan * ...
        [-1, -0.5, -0.2, 0, 0.2, 0.5, 1];
end

validateattributes(options.realGuesses, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'});
validateattributes(options.imaginaryOffsets, {'numeric'}, ...
    {'vector', 'real', 'finite', 'nonempty'});
validateattributes(options.residualTolerance, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'});
validateattributes(options.duplicateTolerance, {'numeric'}, ...
    {'scalar', 'real', 'positive', 'finite'});
validateattributes(options.suppressSolverOutput, {'logical','numeric'}, ...
    {'scalar'});
options.suppressSolverOutput = logical(options.suppressSolverOutput);

candidate = struct('gamma', {}, 'gammaCanonical', {}, 'zeta', {}, ...
    'diagnostics', {}, 'accepted', {});
candidateIndex = 0;
for realGuess = options.realGuesses(:).'
    for imaginaryOffset = options.imaginaryOffsets(:).'
        guess = realGuess + 1i * ...
            (problem.targetS+imaginaryOffset)*problem.omegaStar;
        candidateIndex = candidateIndex + 1;
        try
            if options.suppressSolverOutput
                [~, gamma, zeta, unused3, unused4, unused5, ...
                    unused6, diagnostics] = evalc( ...
                    'faradayFloquet_RT_GenEIG_cylindrical(problem.Ac, problem.omegaStar, problem.R0, problem.m, problem.radialIndex, problem.C, problem.Bd, problem.At, problem.eta, problem.N, problem.modeType, problem.g_sgn, guess, problem.phase)');
                %#ok<NASGU>
            else
                [gamma, zeta, unused3, unused4, unused5, unused6, ...
                    diagnostics] = faradayFloquet_RT_GenEIG_cylindrical( ...
                    problem.Ac, problem.omegaStar, problem.R0, ...
                    problem.m, problem.radialIndex, problem.C, ...
                    problem.Bd, problem.At, problem.eta, problem.N, ...
                    problem.modeType, problem.g_sgn, guess, problem.phase);
                %#ok<NASGU>
            end
            canonicalS = wnl_wrap_quasifrequency( ...
                imag(gamma)/problem.omegaStar);
            gammaCanonical = real(gamma) + ...
                1i*canonicalS*problem.omegaStar;
            accepted = all(isfinite([real(gamma), imag(gamma)])) && ...
                diagnostics.relativeSingularResidual <= ...
                options.residualTolerance;
        catch solveError
            gamma = NaN;
            gammaCanonical = NaN;
            zeta = [];
            diagnostics = struct();
            diagnostics.initialGuess = guess;
            diagnostics.relativeSingularResidual = Inf;
            diagnostics.exitFlag = NaN;
            diagnostics.errorIdentifier = solveError.identifier;
            diagnostics.errorMessage = solveError.message;
            accepted = false;
        end
        candidate(candidateIndex).gamma = gamma;
        candidate(candidateIndex).gammaCanonical = gammaCanonical;
        candidate(candidateIndex).zeta = zeta;
        candidate(candidateIndex).diagnostics = diagnostics;
        candidate(candidateIndex).accepted = accepted;
    end
end

acceptedIndices = find([candidate.accepted]);
if isempty(acceptedIndices)
    residuals = arrayfun(@(x) ...
        x.diagnostics.relativeSingularResidual, candidate);
    [bestResidual, bestIndex] = min(residuals);
    error('vi_dominant_floquet_root:NoConvergedRoot', ...
        ['No root met residualTolerance=%.3e. Best residual %.3e ', ...
         'came from guess %.6g%+.6gi.'], ...
        options.residualTolerance, bestResidual, ...
        real(candidate(bestIndex).diagnostics.initialGuess), ...
        imag(candidate(bestIndex).diagnostics.initialGuess));
end

% Sort accepted candidates by residual, then retain only distinct roots.
acceptedResidual = arrayfun(@(idx) ...
    candidate(idx).diagnostics.relativeSingularResidual, acceptedIndices);
[~, residualOrder] = sort(acceptedResidual, 'ascend');
acceptedIndices = acceptedIndices(residualOrder);
uniqueIndices = zeros(0, 1);
rootScale = max(1.0, problem.omegaStar);
for idx = acceptedIndices
    isDuplicate = false;
    for kept = uniqueIndices(:).'
        difference = candidate(idx).gammaCanonical - ...
            candidate(kept).gammaCanonical;
        if abs(difference) <= options.duplicateTolerance*rootScale
            isDuplicate = true;
            break;
        end
    end
    if ~isDuplicate
        uniqueIndices(end+1, 1) = idx; %#ok<AGROW>
    end
end

uniqueGrowth = arrayfun(@(idx) ...
    real(candidate(idx).gammaCanonical), uniqueIndices);
[~, dominantLocation] = max(uniqueGrowth);
dominantIndex = uniqueIndices(dominantLocation);

result = struct();
result.gamma = candidate(dominantIndex).gamma;
result.gammaCanonical = candidate(dominantIndex).gammaCanonical;
result.zeta = candidate(dominantIndex).zeta;
result.diagnostics = candidate(dominantIndex).diagnostics;
result.candidates = candidate;
result.acceptedIndices = acceptedIndices;
result.uniqueIndices = uniqueIndices;
result.options = options;

if options.verbose
    fprintf('Distinct accepted Floquet roots at Ac=%.8g:\n', problem.Ac);
    for idx = uniqueIndices(:).'
        root = candidate(idx).gammaCanonical;
        residual = ...
            candidate(idx).diagnostics.relativeSingularResidual;
        marker = ' ';
        if idx == dominantIndex
            marker = '*';
        end
        fprintf(' %s %.10g%+.10gi, residual %.3e\n', ...
            marker, real(root), imag(root), residual);
    end
end
end
