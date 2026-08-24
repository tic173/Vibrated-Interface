function resonances = wnl_find_quadratic_resonances(model, ...
    targetModes, augmentedModes, tolerance)
%WNL_FIND_QUADRATIC_RESONANCES Enumerate allowed nonzero Q coefficients.

if nargin < 4
    tolerance = 1.0e-8;
end
if ~iscell(targetModes)
    targetModes = num2cell(targetModes);
end
if ~iscell(augmentedModes)
    augmentedModes = num2cell(augmentedModes);
end

resonances = struct('target', {}, 'input1', {}, 'input2', {}, ...
    'coefficient', {});
for a = 1:numel(targetModes)
    for b = 1:numel(augmentedModes)
        for c = b:numel(augmentedModes)
            out = wnl_combine_spec(model, ...
                {augmentedModes{b}.spec, augmentedModes{c}.spec}, ...
                [1, 1], 'quadratic_candidate');
            if ~wnl_equivalent_spec(out, targetModes{a}.spec)
                continue;
            end
            coefficient = wnl_quadratic_projection(model, ...
                targetModes{a}, augmentedModes{b}, ...
                augmentedModes{c}, b == c);
            if abs(coefficient) <= tolerance
                continue;
            end
            item = struct();
            item.target = targetModes{a}.spec.label;
            item.input1 = augmentedModes{b}.spec.label;
            item.input2 = augmentedModes{c}.spec.label;
            item.coefficient = coefficient;
            resonances(end + 1) = item; %#ok<AGROW>
        end
    end
end
end
