function spec = wnl_spec(m, s, harmonics, ndof, label)
%WNL_SPEC Construct a Floquet-block specification.
%
% harmonics may be an explicit integer vector or a nonnegative cutoff N.
% For s=1/2 and scalar N, the convention is -N-1:N, matching the existing
% Vibrated-Interface subharmonic code. Otherwise the convention is -N:N.

if nargin < 5 || isempty(label)
    label = sprintf('m%d_s%+.6g', m, s);
end

s = wnl_wrap_quasifrequency(s);
if isscalar(harmonics)
    N = harmonics;
    validateattributes(N, {'numeric'}, ...
        {'real', 'finite', 'nonnegative', 'integer', 'scalar'});
    if abs(s - 0.5) < 64 * eps
        harmonics = -N-1:N;
    else
        harmonics = -N:N;
    end
end

harmonics = harmonics(:).';
if any(abs(harmonics - round(harmonics)) > 64 * eps)
    error('wnl_spec:NonintegerHarmonic', ...
        'Temporal harmonic indices must be integers.');
end
if numel(unique(harmonics)) ~= numel(harmonics)
    error('wnl_spec:RepeatedHarmonic', ...
        'Temporal harmonic indices must be unique.');
end

spec = struct();
spec.m = m;
spec.s = s;
% lambda is the continuous Floquet exponent after the wrapped carrier
% i*(n+s)*omega has been removed.  Neutral-point calculations use zero.
% Operating-point calculations overwrite it with the exact modal exponent.
spec.lambda = 0;
spec.n = round(harmonics);
spec.ndof = ndof;
spec.label = char(label);
end
