function uniqueModes = wnl_unique_modes(modes, tolerance)
%WNL_UNIQUE_MODES Remove duplicate null vectors within each Floquet block.
%
% A real axisymmetric harmonic mode may equal its own conjugate. Retaining
% both copies would make a bordered system singular.

if nargin < 2
    tolerance = 1.0e-10;
end
if isempty(modes)
    uniqueModes = {};
    return;
end
if ~iscell(modes)
    modes = num2cell(modes);
end

uniqueModes = {};
for j = 1:numel(modes)
    candidate = modes{j};
    sameBlock = {};
    for k = 1:numel(uniqueModes)
        if wnl_equivalent_spec(uniqueModes{k}.spec, candidate.spec)
            sameBlock{end + 1} = uniqueModes{k}; %#ok<AGROW>
        end
    end
    if isempty(sameBlock)
        uniqueModes{end + 1} = candidate; %#ok<AGROW>
        continue;
    end

    phi = complex(zeros(numel(candidate.vector), numel(sameBlock)));
    for k = 1:numel(sameBlock)
        phi(:, k) = sameBlock{k}.vector / norm(sameBlock{k}.vector);
    end
    candidateVector = candidate.vector / norm(candidate.vector);
    oldRank = rank(full(phi), tolerance);
    newRank = rank(full([phi, candidateVector]), tolerance);
    if newRank > oldRank
        uniqueModes{end + 1} = candidate; %#ok<AGROW>
    end
end
end
