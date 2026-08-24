function coeff = wnl_apply_quadratic(model, x, y, specOut)
%WNL_APPLY_QUADRATIC Evaluate the projected symmetric bilinear action C(x,y).

if ~isfield(model, 'quadratic') || ~isa(model.quadratic, 'function_handle')
    coeff = complex(zeros(specOut.ndof, numel(specOut.n)));
    return;
end
coeff = model.quadratic(x, y, specOut);
if ~isequal(size(coeff), [specOut.ndof, numel(specOut.n)])
    error('wnl_apply_quadratic:BadSize', ...
        'model.quadratic returned an array with the wrong size.');
end
end
