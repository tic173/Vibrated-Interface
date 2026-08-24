function value = wnl_spec_lambda(spec)
%WNL_SPEC_LAMBDA Continuous exponent carried by a Floquet specification.
%
% Older neutral-point specifications do not contain lambda.  They remain
% valid and are interpreted as lambda=0.

value = 0;
if isfield(spec, 'lambda') && ~isempty(spec.lambda)
    value = spec.lambda;
end
if ~isnumeric(value) || ~isscalar(value) || ...
        ~isfinite(real(value)) || ~isfinite(imag(value))
    error('wnl_spec_lambda:BadLambda', ...
        'spec.lambda must be a finite numeric scalar.');
end
end
