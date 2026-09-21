function output = vi_wnl_run_with_input(inputOverride)
%VI_WNL_RUN_WITH_INPUT Run the editable cylinder driver with overrides.

if nargin < 1
    inputOverride = struct();
end
validateattributes(inputOverride,{'struct'},{'scalar'});
viWnlInputOverride = inputOverride; %#ok<NASGU>
driver = fullfile(fileparts(mfilename('fullpath')), ...
    'vi_wnl_user_run_full_cylinder.m');
run(driver);
if ~exist('output','var')
    error('vi_wnl_run_with_input:NoOutput', ...
        'The cylinder WNL driver returned without an output structure.');
end
end
