function coefficient = wnl_quadratic_projection(model, targetMode, ...
    modeB, modeC, identicalInputs)
%WNL_QUADRATIC_PROJECTION Direct resonant coefficient Q_target,BC.

if nargin < 5
    identicalInputs = false;
end
specOut = wnl_combine_spec(model, ...
    {modeB.spec, modeC.spec}, [1, 1], 'quadratic_product');
if ~wnl_equivalent_spec(specOut, targetMode.spec)
    coefficient = 0.0;
    return;
end

factor = 2.0;
if identicalInputs
    factor = 1.0;
end
forcing = factor * wnl_apply_quadratic(model, ...
    modeB.field, modeC.field, targetMode.spec);
coefficient = (targetMode.left' * forcing(:)) / ...
    targetMode.normalization;
end
