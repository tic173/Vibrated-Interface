function modeBar = wnl_conjugate_mode(model, mode, userOpts)
%WNL_CONJUGATE_MODE Construct and normalize the conjugate direct/left mode.

if nargin < 3
    userOpts = struct();
end
opts = wnl_options(userOpts);

directField = wnl_conjugate_field(model, mode.field, ...
    [mode.spec.label, '_bar']);
leftField = wnl_conjugate_field(model, mode.leftField, ...
    [mode.spec.label, '_bar_left']);
specBar = directField.spec;
adjointRowTransform = 'identity equation-row map';
if isfield(model,'conjugateAdjointRows') && ...
        isa(model.conjugateAdjointRows,'function_handle')
    transformed = model.conjugateAdjointRows(leftField,mode.spec);
    if isstruct(transformed)
        leftField = transformed;
    else
        leftField = wnl_make_field(specBar,transformed);
    end
    if ~wnl_equivalent_spec(leftField.spec,specBar)
        error('wnl_conjugate_mode:AdjointRowMapSpec', ...
            ['model.conjugateAdjointRows returned a field in a ', ...
             'different Floquet block.']);
    end
    adjointRowTransform = 'model-specific equation-row map';
end
blockBar = model.block(specBar);
direct = wnl_field_vector(directField);
left = wnl_field_vector(leftField);

normalization = left' * (blockBar.Bslow * direct);
if abs(normalization) < opts.nullTolerance * norm(left) * ...
        max(norm(blockBar.Bslow * direct), eps)
    error('wnl_conjugate_mode:ZeroNormalization', ...
        'The conjugate direct/left descriptor pairing is nearly zero.');
end
left = left / conj(normalization);
leftField = wnl_make_field(specBar, left);

modeBar = struct();
modeBar.spec = specBar;
modeBar.block = blockBar;
modeBar.vector = direct;
modeBar.field = directField;
modeBar.left = left;
modeBar.leftField = leftField;
modeBar.normalization = left' * (blockBar.Bslow * direct);
[modeBar.directResidual,modeBar.directResidualDetails] = ...
    wnl_descriptor_residual( ...
    blockBar.A,blockBar.Bslow,direct,'direct');
[modeBar.leftResidual,modeBar.leftResidualDetails] = ...
    wnl_descriptor_residual( ...
    blockBar.A,blockBar.Bslow,left,'adjoint');
modeBar.adjointConjugation = struct( ...
    'rowTransform',adjointRowTransform, ...
    'usedModelSpecificRowMap', ...
    strcmp(adjointRowTransform,'model-specific equation-row map'));

if opts.verbose
    fprintf(['WNL conjugate %s: direct residual %.3e, left residual ', ...
        '%.3e, adjoint rows = %s\n'], ...
        specBar.label, modeBar.directResidual, modeBar.leftResidual, ...
        adjointRowTransform);
end

end
