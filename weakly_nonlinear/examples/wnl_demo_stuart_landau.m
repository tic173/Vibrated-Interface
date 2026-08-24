function result = wnl_demo_stuart_landau()
%WNL_DEMO_STUART_LANDAU Verify the coefficient engine analytically.
%
% The real two-component system is
%   q_t = J*q + mu*q + g0*(q.'*q)*q,
%   J = [0 -1; 1 0].
%
% With phi=[1;-i] at temporal harmonic n=1, the chosen amplitude
% normalization gives
%   dA/dT = mu*A + 4*g0*A*|A|^2.

g0 = -0.75;
N = 3;
J = [0, -1; 1, 0];

config = struct();
config.omega = 1.0;
config.N = N;
config.ndof = 2;
config.mass = speye(2);
config.linearFourier = @(spec, k) linear_fourier(J, k); %#ok<INUSD>
config.quadraticLocal = @(a, b, sa, sb, so) zeros(2, 1); %#ok<INUSD>
config.cubicLocal = @(a, b, c, sa, sb, sc, so) ...
    cubic_local(g0, a, b, c); %#ok<INUSD>
config.detuningBlock = @(spec, nOut, nIn, direction) ...
    detuning_block(nOut, nIn, direction); %#ok<INUSD>
model = wnl_fourier_model(config);

spec = model.makeSpec(0, 0, 'Stuart_Landau');
directCoeff = complex(zeros(spec.ndof, numel(spec.n)));
idx = find(spec.n == 1, 1);
directCoeff(:, idx) = [1; -1i];
spec.direct = directCoeff(:);
spec.left = directCoeff(:);

opts = struct();
opts.verbose = false;
opts.detuning = 1.0;
result = wnl_analyze_single_mode(model, spec, opts);

expectedMu = 1.0;
expectedG = 4.0 * g0;
assert(abs(result.mu - expectedMu) < 1.0e-10, ...
    'The detuning coefficient test failed.');
assert(abs(result.g - expectedG) < 1.0e-10, ...
    'The cubic coefficient test failed.');
assert(result.validCubicScaling, ...
    'The nonresonant demo was incorrectly marked resonant.');

fprintf('Stuart-Landau verification passed.\n');
fprintf('  mu: computed %.12g, expected %.12g\n', ...
    result.mu, expectedMu);
fprintf('  g : computed %.12g, expected %.12g\n', ...
    result.g, expectedG);
end

function Lk = linear_fourier(J, k)
if k == 0
    Lk = sparse(J);
else
    Lk = sparse(2, 2);
end
end

function value = cubic_local(g0, a, b, c)
value = (g0 / 3.0) * ( ...
    (a.' * b) * c + (a.' * c) * b + (b.' * c) * a);
end

function value = detuning_block(nOut, nIn, direction)
if nOut == nIn
    value = direction * speye(2);
else
    value = sparse(2, 2);
end
end
