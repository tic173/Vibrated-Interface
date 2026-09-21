function scale = wnl_mode_amplitude_scale(value)
%WNL_MODE_AMPLITUDE_SCALE Map internal modal coordinates to peak zeta/h.
%
% Direct cylinder modes are normalized so the complex carrier has unit
% peak interface displacement. A self-conjugate real mode therefore has
% physical peak A, while a distinct Fourier/conjugate pair
%   a*phi + conj(a*phi)
% has peak 2*abs(a). The returned scale S satisfies A_peak=S*a_internal.

if iscell(value)
    scale = cellfun(@wnl_mode_amplitude_scale,value(:));
    return;
end
if numel(value) > 1
    scale = arrayfun(@wnl_mode_amplitude_scale,value(:));
    return;
end

if strcmp(wnl_mode_coordinate_type(value),'real')
    scale = 1.0;
else
    scale = 2.0;
end
end
