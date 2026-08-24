function label = vi_wnl_mode_label(mode)
%VI_WNL_MODE_LABEL Generate a readable label from (m,l,s).

required = {'m','radialIndex','s'};
for fieldIndex = 1:numel(required)
    if ~isfield(mode,required{fieldIndex})
        error('vi_wnl_mode_label:MissingField', ...
            'mode.%s is required.',required{fieldIndex});
    end
end
validateattributes(mode.m,{'numeric'}, ...
    {'scalar','integer','nonnegative','finite'});
validateattributes(mode.radialIndex,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(mode.s,{'numeric'}, ...
    {'scalar','real','finite'});

s = wnl_wrap_quasifrequency(mode.s);
tolerance = 64*eps(max(1,abs(s)));
if abs(s) <= tolerance
    branch = 'harmonic';
elseif abs(s-0.5) <= tolerance
    branch = 'subharmonic';
else
    branch = ['s',number_token(s)];
end
label = sprintf('m%d_l%d_%s',mode.m,mode.radialIndex,branch);
end

function token = number_token(value)
token = lower(sprintf('%.8g',value));
token = strrep(token,'+','');
token = strrep(token,'-','m');
token = strrep(token,'.','p');
end
