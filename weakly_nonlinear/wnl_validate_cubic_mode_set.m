function information = wnl_validate_cubic_mode_set(modes)
%WNL_VALIDATE_CUBIC_MODE_SET Check completeness of the reduced monomials.
%
% The implemented system contains A_j*|A_k|^2 and, when allowed by
% azimuthal/Floquet symmetry, conj(A_j)*A_k^2. Enumerate every cubic
% monomial admitted by those symmetries and reject a mode set when another
% monomial can enter a retained equation. Examples include independent
% radial modes in the same (m,s) class and 1:3 azimuthal resonances.

if ~iscell(modes)
    modes = num2cell(modes);
end
numberOfModes = numel(modes);
specs = cellfun(@mode_spec,modes,'UniformOutput',false);
coordinateTypes = cellfun(@wnl_mode_coordinate_type,modes, ...
    'UniformOutput',false);
isReal = strcmp(coordinateTypes,'real');
components = coordinate_components(specs,isReal);
numberOfComponents = numel(components);
numberOfTriples = nchoosek(numberOfComponents+2,3);
emptyEntry = struct('targetMode',0,'positivePowers',[], ...
    'conjugatePowers',[],'azimuthalWavenumber',0, ...
    'quasifrequency',0);
unsupported = repmat(emptyEntry,numberOfModes*numberOfTriples,1);
unsupportedCount = 0;
for first = 1:numberOfComponents
    for second = first:numberOfComponents
        for third = second:numberOfComponents
            selected = components([first,second,third]);
            positivePowers = zeros(1,numberOfModes);
            conjugatePowers = zeros(1,numberOfModes);
            for componentIndex = 1:3
                modeIndex = selected(componentIndex).modeIndex;
                if selected(componentIndex).conjugated
                    conjugatePowers(modeIndex) = ...
                        conjugatePowers(modeIndex)+1;
                else
                    positivePowers(modeIndex) = ...
                        positivePowers(modeIndex)+1;
                end
            end
            mOut = sum([selected.m]);
            sOut = wnl_wrap_quasifrequency(sum([selected.s]));
            for targetMode = 1:numberOfModes
                if mOut ~= specs{targetMode}.m || ...
                        ~same_quasifrequency(sOut,specs{targetMode}.s)
                    continue;
                end
                if supported_monomial(targetMode,positivePowers, ...
                        conjugatePowers,isReal)
                    continue;
                end
                unsupportedCount = unsupportedCount+1;
                unsupported(unsupportedCount) = struct( ...
                    'targetMode',targetMode, ...
                    'positivePowers',positivePowers, ...
                    'conjugatePowers',conjugatePowers, ...
                    'azimuthalWavenumber',mOut, ...
                    'quasifrequency',sOut);
            end
        end
    end
end
unsupported = unsupported(1:unsupportedCount);

information = struct();
information.supported = isempty(unsupported);
information.unsupportedMonomials = unsupported;
information.coordinateTypes = coordinateTypes;
information.reason = '';
if information.supported
    return;
end

firstUnsupported = unsupported(1);
monomial = format_monomial(firstUnsupported.positivePowers, ...
    firstUnsupported.conjugatePowers,isReal);
information.reason = sprintf( ...
    'Symmetry permits the omitted cubic monomial %s in mode %d.', ...
    monomial,firstUnsupported.targetMode);
error('wnl_validate_cubic_mode_set:IncompleteCubicTensor', ...
    ['The compact g/h amplitude equation is incomplete: %s ', ...
     'Use a mode set without this internal resonance or extend the ', ...
     'reduction to the full cubic coefficient tensor.'], ...
    information.reason);
end

function components = coordinate_components(specs,isReal)
numberOfComponents = numel(specs)+nnz(~isReal);
components = repmat(struct('modeIndex',0,'conjugated',false, ...
    'm',0,'s',0),numberOfComponents,1);
componentIndex = 0;
for modeIndex = 1:numel(specs)
    componentIndex = componentIndex+1;
    components(componentIndex) = struct( ...
        'modeIndex',modeIndex,'conjugated',false, ...
        'm',specs{modeIndex}.m,'s',specs{modeIndex}.s);
    if ~isReal(modeIndex)
        componentIndex = componentIndex+1;
        components(componentIndex) = struct( ...
            'modeIndex',modeIndex,'conjugated',true, ...
            'm',-specs{modeIndex}.m,'s',-specs{modeIndex}.s);
    end
end
end

function tf = supported_monomial(target,positive,conjugate,isReal)
numberOfModes = numel(positive);
tf = false;
for source = 1:numberOfModes
    expectedPositive = zeros(1,numberOfModes);
    expectedConjugate = zeros(1,numberOfModes);
    expectedPositive(target) = expectedPositive(target)+1;
    if isReal(source)
        expectedPositive(source) = expectedPositive(source)+2;
    else
        expectedPositive(source) = expectedPositive(source)+1;
        expectedConjugate(source) = expectedConjugate(source)+1;
    end
    if isequal(positive,expectedPositive) && ...
            isequal(conjugate,expectedConjugate)
        tf = true;
        return;
    end

    if ~isReal(target) && ~isReal(source)
        expectedPositive = zeros(1,numberOfModes);
        expectedConjugate = zeros(1,numberOfModes);
        expectedConjugate(target) = 1;
        expectedPositive(source) = 2;
        if isequal(positive,expectedPositive) && ...
                isequal(conjugate,expectedConjugate)
            tf = true;
            return;
        end
    end
end
end

function tf = same_quasifrequency(first,second)
tf = abs(wnl_wrap_quasifrequency(first-second)) <= ...
    256*eps(max([1,abs(first),abs(second)]));
end

function value = format_monomial(positive,conjugate,isReal)
factors = {};
for modeIndex = 1:numel(positive)
    for power = 1:positive(modeIndex)
        if isReal(modeIndex)
            factors{end+1} = sprintf('x_%d',modeIndex); %#ok<AGROW>
        else
            factors{end+1} = sprintf('A_%d',modeIndex); %#ok<AGROW>
        end
    end
    for power = 1:conjugate(modeIndex)
        factors{end+1} = sprintf('conj(A_%d)',modeIndex); %#ok<AGROW>
    end
end
value = strjoin(factors,'*');
end

function spec = mode_spec(value)
if isstruct(value) && isfield(value,'spec')
    spec = value.spec;
else
    spec = value;
end
if ~isstruct(spec) || ~isfield(spec,'m') || ~isfield(spec,'s')
    error('wnl_validate_cubic_mode_set:MissingSpec', ...
        'Every retained mode must define m and s.');
end
end
