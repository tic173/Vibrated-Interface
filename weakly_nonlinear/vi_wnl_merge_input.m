function merged = vi_wnl_merge_input(defaults,override)
%VI_WNL_MERGE_INPUT Recursively apply a programmatic input override.

if ~isstruct(defaults) || ~isstruct(override) || ...
        ~isscalar(defaults) || ~isscalar(override)
    error('vi_wnl_merge_input:ScalarStructs', ...
        'Both inputs must be scalar structures.');
end
merged = defaults;
names = fieldnames(override);
for fieldIndex = 1:numel(names)
    name = names{fieldIndex};
    value = override.(name);
    if isfield(merged,name) && isstruct(merged.(name)) && ...
            isscalar(merged.(name)) && isstruct(value) && isscalar(value)
        merged.(name) = vi_wnl_merge_input(merged.(name),value);
    else
        merged.(name) = value;
    end
end
end
