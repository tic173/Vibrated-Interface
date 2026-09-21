function gPhysical = wnl_physical_cubic_coefficients(gInternal, modes)
%WNL_PHYSICAL_CUBIC_COEFFICIENTS Express g in peak-displacement units.
%
% If A_k=S_k*a_k, then
%   dA_j/dt = lambda_j*A_j + A_j*sum_k gInternal(j,k)
%              *abs(A_k/S_k)^2.
% Hence only column k is divided by S_k^2.

if iscell(modes)
    numberOfModes = numel(modes);
else
    numberOfModes = numel(modes);
    modes = num2cell(modes);
end
if ~isequal(size(gInternal),[numberOfModes,numberOfModes])
    error('wnl_physical_cubic_coefficients:Size', ...
        'g must be square with one row and column per retained mode.');
end
scales = wnl_mode_amplitude_scale(modes);
gPhysical = gInternal./(abs(scales(:).').^2);
end
