function fieldBar = wnl_conjugate_field(model, field, label)
%WNL_CONJUGATE_FIELD Construct the physical complex-conjugate Floquet field.
%
% If q = phi_n exp(i(n+s)tau+i*m*theta), conjugation requires
%   m_bar = -m,
%   s_bar = wrap(-s),
%   n_bar = -n-s-s_bar.
% For s=1/2 this gives n_bar=-n-1.

if nargin < 3 || isempty(label)
    label = [field.spec.label, '_bar'];
end
if ~isfield(model, 'makeSpec') || ~isa(model.makeSpec, 'function_handle')
    error('wnl_conjugate_field:MissingMakeSpec', ...
        'model.makeSpec(m,s,label) is required.');
end

specIn = field.spec;
specOut = model.makeSpec(-specIn.m, ...
    wnl_wrap_quasifrequency(-specIn.s), label);
specOut.lambda = conj(wnl_spec_lambda(specIn));
coeffOut = complex(zeros(specOut.ndof, numel(specOut.n)));

if specOut.ndof ~= specIn.ndof
    error('wnl_conjugate_field:NdofMismatch', ...
        'Conjugate blocks must have matching spatial dimensions.');
end

for j = 1:numel(specIn.n)
    nOutReal = -specIn.n(j) - specIn.s - specOut.s;
    nOut = round(nOutReal);
    if abs(nOutReal - nOut) > 1.0e-10
        error('wnl_conjugate_field:FrequencyMismatch', ...
            'Could not map conjugate temporal frequency to an integer index.');
    end
    idx = find(specOut.n == nOut, 1);
    if ~isempty(idx)
        coeffOut(:, idx) = conj(field.coeff(:, j));
    end
end

fieldBar = wnl_make_field(specOut, coeffOut);
end
