function matches = wnl_matching_modes(modes, spec)
%WNL_MATCHING_MODES Return retained Floquet modes in a requested block.

if isempty(modes)
    matches = {};
    return;
end
if ~iscell(modes)
    modes = num2cell(modes);
end

matches = {};
for j = 1:numel(modes)
    if wnl_equivalent_spec(modes{j}.spec, spec)
        matches{end + 1} = modes{j}; %#ok<AGROW>
    end
end
matches = wnl_unique_modes(matches);
end
