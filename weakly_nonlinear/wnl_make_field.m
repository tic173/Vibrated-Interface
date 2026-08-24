function field = wnl_make_field(spec, coeff)
%WNL_MAKE_FIELD Store Fourier coefficients with their Floquet metadata.
%
% coeff has size [spec.ndof, numel(spec.n)]. A stacked vector is accepted
% and reshaped using MATLAB column-major ordering.

nt = numel(spec.n);
if isvector(coeff) && numel(coeff) == spec.ndof * nt
    coeff = reshape(coeff, spec.ndof, nt);
end
if ~isequal(size(coeff), [spec.ndof, nt])
    error('wnl_make_field:BadSize', ...
        'Expected coefficient size [%d,%d], received [%d,%d].', ...
        spec.ndof, nt, size(coeff, 1), size(coeff, 2));
end

field = struct();
field.spec = spec;
field.coeff = coeff;
end
