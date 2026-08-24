function coeff = wnl_apply_cubic(model, x, y, z, specOut)
%WNL_APPLY_CUBIC Evaluate the projected symmetric trilinear action D(x,y,z).

if ~isfield(model, 'cubic') || ~isa(model.cubic, 'function_handle')
    coeff = complex(zeros(specOut.ndof, numel(specOut.n)));
    return;
end
coeff = model.cubic(x, y, z, specOut);
if ~isequal(size(coeff), [specOut.ndof, numel(specOut.n)])
    error('wnl_apply_cubic:BadSize', ...
        'model.cubic returned an array with the wrong size.');
end
end
