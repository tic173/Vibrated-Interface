function result = vi_operating_point_floquet_mode(parameters, mode, ...
    temporalCutoff, rootOptions, precomputedRoot)
%VI_OPERATING_POINT_FLOQUET_MODE Exact reduced mode at a requested forcing.
%
% The returned periodic carrier and continuous exponent satisfy
%   q(t)=exp(lambda*t)*sum_n q_n exp(i*(n+s)*omega*t).
% Unlike a neutral-point seed, lambda need not be zero.  This permits two
% retained modes with different individual critical accelerations to be
% represented at one common physical operating amplitude.

if nargin < 4
    rootOptions = struct();
end
if nargin < 5
    precomputedRoot = [];
end
required = {'aAnalysis', 'omegaStar', 'R0', 'C', 'Bd', 'At', ...
    'eta', 'g_sgn', 'phase'};
for fieldIndex = 1:numel(required)
    if ~isfield(parameters, required{fieldIndex})
        error('vi_operating_point_floquet_mode:MissingParameter', ...
            'parameters.%s is required.', required{fieldIndex});
    end
end

targetS = wnl_wrap_quasifrequency(mode.s);
if abs(targetS) < 1.0e-12
    modeType = 'H';
elseif abs(targetS-0.5) < 1.0e-12
    modeType = 'SH';
else
    error('vi_operating_point_floquet_mode:UnsupportedGeneralS', ...
        ['The reduced cylinder root search currently supports the ', ...
         'harmonic (s=0) and subharmonic (s=1/2) branches.']);
end

problem = struct();
problem.Ac = parameters.aAnalysis;
problem.omegaStar = parameters.omegaStar;
problem.R0 = parameters.R0;
problem.m = mode.m;
problem.radialIndex = mode.radialIndex;
problem.C = parameters.C;
problem.Bd = parameters.Bd;
problem.At = parameters.At;
problem.eta = parameters.eta;
problem.N = temporalCutoff;
problem.modeType = modeType;
problem.g_sgn = parameters.g_sgn;
problem.targetS = targetS;
problem.phase = parameters.phase;

if isempty(precomputedRoot)
    root = vi_dominant_floquet_root(problem, rootOptions);
else
    if ~isstruct(precomputedRoot) || ...
            ~isfield(precomputedRoot, 'gammaCanonical')
        error('vi_operating_point_floquet_mode:BadPrecomputedRoot', ...
            'precomputedRoot must contain gammaCanonical.');
    end
    root = precomputedRoot;
end
gamma = root.gammaCanonical;
s = wnl_wrap_quasifrequency(imag(gamma)/parameters.omegaStar);
lambda = gamma-1i*s*parameters.omegaStar;
if abs(imag(lambda)) <= 1.0e-10*max(1, abs(gamma))
    lambda = real(lambda);
end

[zeta, diagnostics] = vi_reduced_cylinder_mode_at_exponent( ...
    parameters.aAnalysis, parameters.omegaStar, parameters.R0, ...
    mode.m, mode.radialIndex, parameters.C, parameters.Bd, ...
    parameters.At, parameters.eta, temporalCutoff, modeType, ...
    parameters.g_sgn, gamma, parameters.phase);
zeta = zeta(:, 1);
if strcmp(modeType, 'SH')
    harmonicIndices = -temporalCutoff-1:temporalCutoff;
else
    harmonicIndices = -temporalCutoff:temporalCutoff;
end

result = struct();
result.label = mode.label;
result.modeType = modeType;
result.analysisAmplitude = parameters.aAnalysis;
result.gamma = gamma;
result.lambda = lambda;
result.s = s;
result.zeta = zeta(:);
result.harmonicIndices = harmonicIndices;
result.diagnostics = diagnostics;
result.rootSearch = root;
result.problem = problem;
end
