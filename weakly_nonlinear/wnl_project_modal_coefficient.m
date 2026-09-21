function [value, diagnostics] = wnl_project_modal_coefficient( ...
        rawValue, targetMode, userOpts)
%WNL_PROJECT_MODAL_COEFFICIENT Enforce reality of a real modal equation.

if nargin < 3
    userOpts = struct();
end
opts = wnl_options(userOpts);
coordinateType = wnl_mode_coordinate_type( ...
    targetMode,opts.realModeEigenvalueTolerance);
diagnostics = struct('coordinateType',coordinateType, ...
    'rawValue',rawValue,'imaginaryLeakage',0, ...
    'allowedImaginaryPart',Inf,'accepted',true);
value = rawValue;
if ~strcmp(coordinateType,'real')
    return;
end

diagnostics.imaginaryLeakage = abs(imag(rawValue))/max(abs(rawValue),eps);
diagnostics.allowedImaginaryPart = ...
    opts.realCoefficientImaginaryTolerance*max(1,abs(real(rawValue)));
diagnostics.accepted = isfinite(rawValue) && ...
    abs(imag(rawValue)) <= diagnostics.allowedImaginaryPart;
if diagnostics.accepted
    value = real(rawValue);
else
    value = NaN;
end
end
