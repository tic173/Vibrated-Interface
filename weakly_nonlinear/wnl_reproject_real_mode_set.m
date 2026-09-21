function result = wnl_reproject_real_mode_set(result,userOpts)
%WNL_REPROJECT_REAL_MODE_SET Reapply the reality gate to completed m=0 data.
% This function changes no direct, adjoint, nonlinear-action, or forced
% field. It only projects already-computed raw coefficients onto the real
% signed coordinates required for self-conjugate axisymmetric modes.

if nargin < 2
    userOpts = struct();
end
opts = wnl_options(userOpts);
required = {'modes','self','cross','g','forcedSolvesValid', ...
    'forcedSolvesExploratoryUsable'};
for fieldIndex = 1:numel(required)
    if ~isfield(result,required{fieldIndex})
        error('wnl_reproject_real_mode_set:MissingField', ...
            'result.%s is required.',required{fieldIndex});
    end
end
modes = result.modes(:);
nModes = numel(modes);
coordinateTypes = cellfun(@(mode) wnl_mode_coordinate_type( ...
    mode,opts.realModeEigenvalueTolerance, ...
    opts.realModeMinimumConjugateOverlap),modes,'UniformOutput',false);
if any(~strcmp(coordinateTypes,'real'))
    error('wnl_reproject_real_mode_set:ComplexMode', ...
        ['Saved-coefficient reprojection is restricted to mode sets whose ', ...
         'coordinates are all self-conjugate and real.']);
end
if ~all(result.forcedSolvesValid(:))
    error('wnl_reproject_real_mode_set:UnconvergedForcedField', ...
        ['A forced field missed its strict residual gate. Reprojection ', ...
         'cannot repair an unconverged slaved field.']);
end
if isfield(result,'quadraticResonances') && ...
        ~isempty(result.quadraticResonances)
    error('wnl_reproject_real_mode_set:QuadraticResonance', ...
        'A quadratic resonance prevents a cubic-only reprojection.');
end

g = complex(nan(nModes));
valid = false(nModes);
for targetIndex = 1:nModes
    entry = result.self{targetIndex};
    entry = reproject_entry(entry,modes{targetIndex},opts,'self');
    result.self{targetIndex} = entry;
    g(targetIndex,targetIndex) = entry.g;
    valid(targetIndex,targetIndex) = entry.validCubicScaling;
    for sourceIndex = 1:nModes
        if sourceIndex == targetIndex
            continue;
        end
        entry = result.cross{targetIndex,sourceIndex};
        entry = reproject_entry(entry,modes{targetIndex},opts,'cross');
        result.cross{targetIndex,sourceIndex} = entry;
        g(targetIndex,sourceIndex) = entry.g;
        valid(targetIndex,sourceIndex) = entry.validCubicScaling;
    end
end

result.g = g;
result.gPhysicalPeak = wnl_physical_cubic_coefficients(g,modes);
result.validCubicScaling = valid;
result.coordinateTypes = coordinateTypes;
result.numericalCoefficientValidity = all(valid(:)) && ...
    all(isfinite(real(g(:)))) && all(isfinite(imag(g(:))));
result.exploratoryCoefficientAvailability = ...
    all(result.forcedSolvesExploratoryUsable(:)) && ...
    result.numericalCoefficientValidity;
slowEnvelopeValid = true;
if isfield(result,'slowEnvelopeValid')
    slowEnvelopeValid = logical(result.slowEnvelopeValid);
end
result.quantitativelyValid = ...
    result.numericalCoefficientValidity && slowEnvelopeValid;
operatingPoint = isfield(result,'referenceType') && ...
    strcmp(result.referenceType,'operating-point Floquet reduction');
result.smallAmplitudeTransientAvailable = ...
    operatingPoint && result.numericalCoefficientValidity;
result.smallAmplitudeTransientExploratoryAvailable = ...
    result.smallAmplitudeTransientAvailable;
result.coefficientRealityReprojection = struct( ...
    'applied',true, ...
    'realCoefficientImaginaryTolerance', ...
    opts.realCoefficientImaginaryTolerance, ...
    'doesNotRecomputeFields',true);
end

function entry = reproject_entry(entry,targetMode,opts,entryType)
if ~isstruct(entry) || ~isfield(entry,'gUnprojected') || ...
        ~isfinite(entry.gUnprojected)
    error('wnl_reproject_real_mode_set:MissingRawCoefficient', ...
        'The saved %s entry has no finite raw coefficient.',entryType);
end
[entry.g,reality] = wnl_project_modal_coefficient( ...
    entry.gUnprojected,targetMode,opts);
entry.coefficientReality = reality;
quadraticResonance = isfield(entry,'quadraticResonance') && ...
    logical(entry.quadraticResonance);
forcedValid = isfield(entry,'forcedSolvesValid') && ...
    logical(entry.forcedSolvesValid);
entry.validCubicScaling = ...
    ~quadraticResonance && forcedValid && reality.accepted;
if entry.validCubicScaling
    entry.message = '';
else
    entry.message = sprintf([ ...
        'A real modal coefficient had imaginary part %.3e, above the ', ...
        'allowed %.3e.'],abs(imag(entry.gUnprojected)), ...
        reality.allowedImaginaryPart);
end
end
