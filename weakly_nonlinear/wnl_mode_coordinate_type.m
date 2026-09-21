function type = wnl_mode_coordinate_type( ...
        value, eigenvalueTolerance, minimumConjugateOverlap)
%WNL_MODE_COORDINATE_TYPE Classify a retained Floquet amplitude coordinate.
%
% A mode is represented by one real signed amplitude only when its
% azimuthal/Floquet block is invariant under physical conjugation and its
% continuous exponent is real to numerical accuracy. All other modes use a
% complex amplitude together with a distinct conjugate field.

if nargin < 2 || isempty(eigenvalueTolerance)
    eigenvalueTolerance = 1.0e-4;
end
if nargin < 3 || isempty(minimumConjugateOverlap)
    minimumConjugateOverlap = 0.99;
end
validateattributes(eigenvalueTolerance, {'numeric'}, ...
    {'scalar','real','nonnegative','finite'});
validateattributes(minimumConjugateOverlap, {'numeric'}, ...
    {'scalar','real','>=',0,'<=',1,'finite'});

if isstruct(value) && isfield(value,'coordinateType') && ...
        ~isempty(value.coordinateType)
    type = validatestring(value.coordinateType, {'real','complex'});
    return;
end
if isstruct(value) && isfield(value,'spec')
    spec = value.spec;
else
    spec = value;
end
if ~isstruct(spec) || ~isfield(spec,'m') || ~isfield(spec,'s')
    type = 'complex';
    return;
end

blockConjugacyError = abs(wnl_wrap_quasifrequency(-spec.s-spec.s));
selfConjugateBlock = spec.m == 0 && ...
    blockConjugacyError <= 256*eps(max(1,abs(spec.s)));
lambda = wnl_spec_lambda(spec);
approximatelyRealExponent = abs(imag(lambda)) <= ...
    eigenvalueTolerance*max(1,abs(lambda));

if selfConjugateBlock && approximatelyRealExponent
    type = 'real';
else
    type = 'complex';
end
if strcmp(type,'real') && isstruct(value) && isfield(value,'field') && ...
        isstruct(value.field) && isfield(value.field,'coeff')
    conjugateCoeff = self_block_conjugate(value.field);
    original = value.field.coeff(:);
    overlap = abs(original'*conjugateCoeff(:))/ ...
        max(norm(original)*norm(conjugateCoeff(:)),eps);
    if overlap < minimumConjugateOverlap
        type = 'complex';
    end
end
end

function coeffBar = self_block_conjugate(field)
spec = field.spec;
coeffBar = complex(zeros(size(field.coeff)));
for sourceIndex = 1:numel(spec.n)
    targetReal = -spec.n(sourceIndex)-2*spec.s;
    targetHarmonic = round(targetReal);
    if abs(targetReal-targetHarmonic) > 1.0e-10
        continue;
    end
    targetIndex = find(spec.n == targetHarmonic,1);
    if ~isempty(targetIndex)
        coeffBar(:,targetIndex) = conj(field.coeff(:,sourceIndex));
    end
end
end
