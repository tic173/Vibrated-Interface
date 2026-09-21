function report = vi_compare_shared_mode_invariance(files,options)
%VI_COMPARE_SHARED_MODE_INVARIANCE Compare duplicated modes across caches.

if nargin < 2
    options = struct();
end
options = defaults(options);
if ischar(files) || isstring(files)
    files = cellstr(files);
end
validateattributes(files,{'cell'},{'vector','nonempty'});
records = cell(numel(files),1);
for fileIndex = 1:numel(files)
    filename = char(files{fileIndex});
    if ~isfile(filename)
        error('vi_compare_shared_mode_invariance:MissingFile', ...
            'Missing coefficient cache: %s',filename);
    end
    saved = load(filename,'output');
    records{fileIndex} = saved.output;
end

operatorComparisons = compare_operators(records,files,options);
occurrences = collect_occurrences(records,files);
labels = unique({occurrences.label},'stable');
shared = struct([]);
for labelIndex = 1:numel(labels)
    selected = occurrences(strcmp({occurrences.label},labels{labelIndex}));
    if numel(selected) < 2
        continue;
    end
    for firstIndex = 1:numel(selected)-1
        for secondIndex = firstIndex+1:numel(selected)
            entry = compare_occurrences( ...
                selected(firstIndex),selected(secondIndex),options);
            if isempty(shared)
                shared = entry;
            else
                shared(end+1,1) = entry; %#ok<AGROW>
            end
        end
    end
end
if isempty(shared)
    error('vi_compare_shared_mode_invariance:NoSharedModes', ...
        'No mode label occurs in more than one supplied cache.');
end
report = struct();
report.schemaVersion = 'V1-common-operator-shared-mode-audit';
report.files = files(:);
report.operatorComparisons = operatorComparisons;
report.sharedModes = shared;
report.tolerances = options;
report.operatorInvariant = all([operatorComparisons.passed]);
report.sharedModesInvariant = all([shared.passed]);
report.passed = report.operatorInvariant && report.sharedModesInvariant;
print_report(report);
if options.assertPassed && ~report.passed
    error('vi_compare_shared_mode_invariance:Failed', ...
        'The common-operator shared-mode invariance audit failed.');
end
end

function options = defaults(options)
defaultValues.operatorTolerance = 1.0e-13;
defaultValues.linearSymmetricRelativeTolerance = 1.0e-7;
defaultValues.selfSymmetricRelativeTolerance = 5.0e-3;
defaultValues.minimumModeOverlap = 1-1.0e-7;
defaultValues.assertPassed = false;
names = fieldnames(defaultValues);
for nameIndex = 1:numel(names)
    name = names{nameIndex};
    if ~isfield(options,name)
        options.(name) = defaultValues.(name);
    end
end
end

function comparisons = compare_operators(records,files,options)
reference = records{1};
comparisons = repmat(struct(),numel(records)-1,1);
for recordIndex = 2:numel(records)
    candidate = records{recordIndex};
    radialA = reference.operatorMetadata.radialGrid;
    radialB = candidate.operatorMetadata.radialGrid;
    entry.referenceFile = files{1};
    entry.candidateFile = files{recordIndex};
    entry.radialNodeDifference = relative_difference( ...
        radialA.r,radialB.r);
    entry.firstDerivativeDifference = relative_difference( ...
        radialA.D,radialB.D);
    entry.secondDerivativeDifference = relative_difference( ...
        radialA.D2,radialB.D2);
    entry.temporalCutoffEqual = reference.parameters.numerics.N == ...
        candidate.parameters.numerics.N;
    entry.azimuthalGridEqual = ...
        reference.parameters.numerics.Ntheta == ...
        candidate.parameters.numerics.Ntheta;
    entry.lowerGridDifference = relative_difference( ...
        reference.operatorMetadata.verticalGrid.lower.z, ...
        candidate.operatorMetadata.verticalGrid.lower.z);
    entry.upperGridDifference = relative_difference( ...
        reference.operatorMetadata.verticalGrid.upper.z, ...
        candidate.operatorMetadata.verticalGrid.upper.z);
    entry.boundaryEqual = strcmp( ...
        reference.operatorMetadata.sidewallTangentialCondition, ...
        candidate.operatorMetadata.sidewallTangentialCondition) && ...
        strcmp(reference.operatorMetadata.contactLine, ...
        candidate.operatorMetadata.contactLine);
    numericDifferences = [entry.radialNodeDifference, ...
        entry.firstDerivativeDifference, ...
        entry.secondDerivativeDifference,entry.lowerGridDifference, ...
        entry.upperGridDifference];
    entry.passed = all(numericDifferences <= ...
        options.operatorTolerance) && entry.temporalCutoffEqual && ...
        entry.azimuthalGridEqual && entry.boundaryEqual;
    comparisons(recordIndex-1,1) = entry;
end
end

function occurrences = collect_occurrences(records,files)
occurrences = struct([]);
for fileIndex = 1:numel(records)
    wnl = records{fileIndex}.weaklyNonlinear;
    for modeIndex = 1:numel(wnl.modes)
        if ~isfield(wnl,'self') || numel(wnl.self) < modeIndex || ...
                ~isstruct(wnl.self{modeIndex}) || ...
                ~isfield(wnl.self{modeIndex},'g')
            continue;
        end
        entry.label = wnl.modes{modeIndex}.spec.label;
        entry.file = files{fileIndex};
        entry.modeIndex = modeIndex;
        entry.lambda = wnl.linearCoefficients(modeIndex);
        entry.g = wnl.self{modeIndex}.g;
        entry.gPhysicalPeak = wnl.gPhysicalPeak(modeIndex,modeIndex);
        entry.vector = wnl.modes{modeIndex}.vector(:);
        if isempty(occurrences)
            occurrences = entry;
        else
            occurrences(end+1,1) = entry; %#ok<AGROW>
        end
    end
end
end

function entry = compare_occurrences(a,b,options)
entry.label = a.label;
entry.firstFile = a.file;
entry.secondFile = b.file;
entry.lambda = [a.lambda,b.lambda];
entry.g = [a.g,b.g];
entry.gPhysicalPeak = [a.gPhysicalPeak,b.gPhysicalPeak];
entry.linearSymmetricRelativeDifference = symmetric_difference( ...
    a.lambda,b.lambda);
entry.selfSymmetricRelativeDifference = symmetric_difference(a.g,b.g);
entry.physicalSelfSymmetricRelativeDifference = symmetric_difference( ...
    a.gPhysicalPeak,b.gPhysicalPeak);
entry.modeOverlap = abs(a.vector'*b.vector) / ...
    max(norm(a.vector)*norm(b.vector),eps);
entry.passed = entry.linearSymmetricRelativeDifference <= ...
    options.linearSymmetricRelativeTolerance && ...
    entry.selfSymmetricRelativeDifference <= ...
    options.selfSymmetricRelativeTolerance && ...
    entry.physicalSelfSymmetricRelativeDifference <= ...
    options.selfSymmetricRelativeTolerance && ...
    entry.modeOverlap >= options.minimumModeOverlap;
end

function value = relative_difference(a,b)
if ~isequal(size(a),size(b))
    value = Inf;
    return;
end
value = norm(a(:)-b(:))/max([norm(a(:)),norm(b(:)),eps]);
end

function value = symmetric_difference(a,b)
value = 2*abs(a-b)/max(abs(a)+abs(b),eps);
end

function print_report(report)
fprintf('\nCommon-operator shared-mode invariance audit\n');
fprintf('  complete operator invariant = %d\n',report.operatorInvariant);
for index = 1:numel(report.operatorComparisons)
    item = report.operatorComparisons(index);
    fprintf(['  operator comparison %d: D %.3e, D2 %.3e, ', ...
        'lower/upper z %.3e/%.3e, passed=%d\n'],index, ...
        item.firstDerivativeDifference,item.secondDerivativeDifference, ...
        item.lowerGridDifference,item.upperGridDifference,item.passed);
end
for index = 1:numel(report.sharedModes)
    item = report.sharedModes(index);
    fprintf(['  %-22s lambda %.12g / %.12g (spread %.3e), ', ...
        'g_self %.12g / %.12g (spread %.3e), overlap %.12g, ', ...
        'passed=%d\n'],item.label,real(item.lambda(1)), ...
        real(item.lambda(2)),item.linearSymmetricRelativeDifference, ...
        real(item.g(1)),real(item.g(2)), ...
        item.selfSymmetricRelativeDifference,item.modeOverlap,item.passed);
end
fprintf('  overall passed = %d\n\n',report.passed);
end
