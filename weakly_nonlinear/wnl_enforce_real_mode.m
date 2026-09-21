function [mode, diagnostics] = wnl_enforce_real_mode( ...
        model, mode, modeBar, userOpts)
%WNL_ENFORCE_REAL_MODE Project a self-conjugate Floquet pair onto real data.

if nargin < 4
    userOpts = struct();
end
opts = wnl_options(userOpts);
if ~strcmp(wnl_mode_coordinate_type( ...
        mode,opts.realModeEigenvalueTolerance, ...
        opts.realModeMinimumConjugateOverlap),'real')
    error('wnl_enforce_real_mode:NotRealCoordinate', ...
        'Mode %s is not a self-conjugate real Floquet coordinate.', ...
        mode.spec.label);
end

direct = mode.vector(:);
directBar = modeBar.vector(:);
directOverlap = abs(direct'*directBar)/max(norm(direct)*norm(directBar),eps);
if directOverlap < opts.realModeMinimumConjugateOverlap
    error('wnl_enforce_real_mode:ConjugateMismatch', ...
        ['Mode %s occupies a self-conjugate block but its direct field ', ...
         'overlap with the physical conjugate is only %.6g. Treat this ', ...
         'as a complex pair or improve the eigensolve.'], ...
        mode.spec.label,directOverlap);
end

phaseRelation = (direct'*directBar)/max(real(direct'*direct),eps);
alignment = exp(0.5i*angle(phaseRelation));
direct = 0.5*(alignment*direct+conj(alignment)*directBar);
left = 0.5*(alignment*mode.left(:)+conj(alignment)*modeBar.left(:));

spec = mode.spec;
originalLambda = wnl_spec_lambda(spec);
spec.lambda = real(originalLambda);
field = wnl_make_field(spec,direct);
if ~isempty(opts.normalizeDirect)
    normalized = opts.normalizeDirect(field);
    if isstruct(normalized)
        direct = wnl_field_vector(normalized);
    else
        direct = normalized(:);
    end
end
block = model.block(spec);
normalization = left'*(block.Bslow*direct);
if abs(normalization) < opts.nullTolerance*norm(left)* ...
        max(norm(block.Bslow*direct),eps)
    error('wnl_enforce_real_mode:ZeroNormalization', ...
        'The real direct/adjoint descriptor pairing is nearly zero.');
end
left = left/conj(normalization);

mode.spec = spec;
mode.block = block;
mode.vector = direct;
mode.field = wnl_make_field(spec,direct);
mode.left = left;
mode.leftField = wnl_make_field(spec,left);
mode.normalization = left'*(block.Bslow*direct);
[mode.directResidual,mode.directResidualDetails] = ...
    wnl_descriptor_residual(block.A,block.Bslow,direct,'direct');
[mode.leftResidual,mode.leftResidualDetails] = ...
    wnl_descriptor_residual(block.A,block.Bslow,left,'adjoint');
mode.coordinateType = 'real';

conjugated = wnl_conjugate_field(model,mode.field, ...
    [mode.spec.label,'_reality_check']);
realityResidual = norm(wnl_field_vector(conjugated)-direct)/ ...
    max(norm(direct),eps);
diagnostics = struct('applied',true, ...
    'originalLambda',originalLambda,'projectedLambda',spec.lambda, ...
    'directConjugateOverlap',directOverlap, ...
    'fieldRealityResidual',realityResidual, ...
    'phaseAlignment',alignment);
if ~isfield(mode,'tracking') || ~isstruct(mode.tracking)
    mode.tracking = struct();
end
mode.tracking.realCoordinateProjection = diagnostics;

if opts.verbose
    fprintf(['WNL real-coordinate projection %s: Im(lambda) %.3e -> 0, ', ...
        'conjugate overlap %.9f, field reality residual %.3e\n'], ...
        spec.label,imag(originalLambda),directOverlap,realityResidual);
end
end
