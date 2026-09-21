function report = vi_compact_saved_output( ...
        sourceFile,destinationFile,profile)
%VI_COMPACT_SAVED_OUTPUT Create a smaller restart/postprocessing cache.
%
% REPORT = VI_COMPACT_SAVED_OUTPUT(SOURCE) writes SOURCE_compact.mat.
% The source file is never deleted or overwritten. Supply a distinct
% destination path to choose another name.

if nargin < 1 || isempty(sourceFile)
    error('vi_compact_saved_output:SourceRequired', ...
        'A source MAT-file is required.');
end
sourceFile = char(sourceFile);
if ~isfile(sourceFile)
    error('vi_compact_saved_output:SourceMissing', ...
        'The source MAT-file does not exist: %s',sourceFile);
end
if nargin < 3 || isempty(profile)
    profile = 'compact';
end
profile = validatestring(profile,{'compact','full'});

[sourceDirectory,sourceName,sourceExtension] = fileparts(sourceFile);
if isempty(sourceExtension)
    sourceExtension = '.mat';
end
if nargin < 2 || isempty(destinationFile)
    destinationFile = fullfile(sourceDirectory, ...
        [sourceName,'_compact',sourceExtension]);
else
    destinationFile = char(destinationFile);
end
if strcmp(canonical_path(sourceFile),canonical_path(destinationFile))
    error('vi_compact_saved_output:RefuseOverwrite', ...
        ['The compact destination must differ from the source. Verify the ', ...
         'new cache before replacing or deleting the original.']);
end
if isfile(destinationFile)
    error('vi_compact_saved_output:DestinationExists', ...
        'The destination already exists: %s',destinationFile);
end

loaded = load(sourceFile,'output');
if ~isfield(loaded,'output') || ~isstruct(loaded.output)
    error('vi_compact_saved_output:MissingOutput', ...
        'The source MAT-file does not contain an output structure.');
end
output = loaded.output;
if ~isfield(output,'input') || ~isstruct(output.input)
    output.input = struct();
end
if ~isfield(output.input,'run') || ~isstruct(output.input.run)
    output.input.run = struct();
end
output.input.run.saveProfile = profile;

repositoryRoot = sourceDirectory;
[~,savedFile,saveInformation] = vi_save_output_record( ...
    output,destinationFile,repositoryRoot);
sourceDetails = dir(sourceFile);
destinationDetails = dir(savedFile);

verification = load(savedFile,'output');
if ~isfield(verification.output,'storage') || ...
        ~strcmp(verification.output.storage.profile,profile)
    error('vi_compact_saved_output:VerificationFailed', ...
        'The saved cache does not contain the requested storage profile.');
end
equivalence = verify_retained_data(output,verification.output);

report = struct();
report.sourceFile = sourceFile;
report.destinationFile = savedFile;
report.profile = profile;
report.sourceBytes = sourceDetails.bytes;
report.destinationBytes = destinationDetails.bytes;
report.reductionFraction = ...
    1-destinationDetails.bytes/max(sourceDetails.bytes,1);
report.saveInformation = saveInformation;
report.storage = verification.output.storage;
report.verification = equivalence;
fprintf('Compact cache saved: %s\n',savedFile);
fprintf('  %.1f MiB -> %.1f MiB (%.1f%% smaller)\n', ...
    report.sourceBytes/2^20,report.destinationBytes/2^20, ...
    100*report.reductionFraction);
end

function value = canonical_path(value)
value = char(java.io.File(value).getCanonicalPath());
end

function report = verify_retained_data(original,compact)
compare_field(original,compact,'input','output');
compare_field(original,compact,'parameters','output');
compare_field(original,compact,'linear','output');
if ~isfield(original,'weaklyNonlinear') || ...
        ~isfield(compact,'weaklyNonlinear')
    error('vi_compact_saved_output:VerificationFailed', ...
        'The source or compact output lacks weakly nonlinear data.');
end
first = original.weaklyNonlinear;
second = compact.weaklyNonlinear;
coefficientFields = {'linearCoefficients','linearCoefficient','mu', ...
    'g','gPhysicalPeak','phaseSensitiveG', ...
    'phaseSensitiveGPhysicalPeak','forcedSolvesValid', ...
    'forcedSolvesExploratoryUsable','validCubicScaling'};
for fieldIndex = 1:numel(coefficientFields)
    compare_field(first,second,coefficientFields{fieldIndex}, ...
        'weaklyNonlinear');
end

firstModes = retained_modes(first);
secondModes = retained_modes(second);
if numel(firstModes) ~= numel(secondModes)
    verification_error('retained mode count changed');
end
for modeIndex = 1:numel(firstModes)
    compare_field(firstModes{modeIndex},secondModes{modeIndex}, ...
        'vector',sprintf('mode %d',modeIndex));
    compare_field(firstModes{modeIndex},secondModes{modeIndex}, ...
        'left',sprintf('mode %d',modeIndex));
    compare_field(firstModes{modeIndex}.field, ...
        secondModes{modeIndex}.field,'coeff', ...
        sprintf('mode %d field',modeIndex));
    specFields = {'m','s','lambda','n','ndof','label', ...
        'radialIndex','betaStar'};
    for fieldIndex = 1:numel(specFields)
        compare_field(firstModes{modeIndex}.spec, ...
            secondModes{modeIndex}.spec,specFields{fieldIndex}, ...
            sprintf('mode %d spec',modeIndex));
    end
end

verify_coefficient_cells(first,second,'self',true);
verify_coefficient_cells(first,second,'cross',false);
report = struct('passed',true,'numberOfModes',numel(firstModes), ...
    'coefficientFieldsChecked',{coefficientFields}, ...
    'exactNumericEquality',true);
end

function verify_coefficient_cells(first,second,name,isSelf)
if ~isfield(first,name)
    return;
end
if ~isfield(second,name) || ~isequal(size(first.(name)),size(second.(name)))
    verification_error(sprintf('%s coefficient-cell dimensions changed',name));
end
for entryIndex = 1:numel(first.(name))
    originalEntry = first.(name){entryIndex};
    compactEntry = second.(name){entryIndex};
    if isempty(originalEntry)
        if ~isempty(compactEntry)
            verification_error(sprintf('%s{%d} changed from empty', ...
                name,entryIndex));
        end
        continue;
    end
    scalarFields = {'g','gUnprojected','validCubicScaling', ...
        'quadraticResonance','forcedSolvesValid', ...
        'forcedSolvesExploratoryUsable'};
    for fieldIndex = 1:numel(scalarFields)
        compare_field(originalEntry,compactEntry,scalarFields{fieldIndex}, ...
            sprintf('%s{%d}',name,entryIndex));
    end
    if isSelf
        retainedForcedFields = {'qAA','qAbarA'};
    else
        retainedForcedFields = {'qAB','qAbarB'};
    end
    for fieldIndex = 1:numel(retainedForcedFields)
        forcedName = retainedForcedFields{fieldIndex};
        if ~isfield(originalEntry,forcedName) || ...
                isempty(originalEntry.(forcedName))
            continue;
        end
        if ~isfield(compactEntry,forcedName) || ...
                ~isfield(compactEntry.(forcedName),'field')
            verification_error(sprintf('%s{%d}.%s was not retained', ...
                name,entryIndex,forcedName));
        end
        compare_field(originalEntry.(forcedName).field, ...
            compactEntry.(forcedName).field,'coeff', ...
            sprintf('%s{%d}.%s.field',name,entryIndex,forcedName));
    end
end
end

function modes = retained_modes(wnl)
if isfield(wnl,'modes') && iscell(wnl.modes)
    modes = wnl.modes(:);
elseif isfield(wnl,'mode') && isstruct(wnl.mode)
    modes = {wnl.mode};
else
    modes = {};
end
end

function compare_field(first,second,name,context)
firstHas = isfield(first,name);
secondHas = isfield(second,name);
if firstHas ~= secondHas
    verification_error(sprintf('%s.%s presence changed',context,name));
end
if firstHas && ~isequaln(first.(name),second.(name))
    verification_error(sprintf('%s.%s changed',context,name));
end
end

function verification_error(message)
error('vi_compact_saved_output:VerificationFailed','%s',message);
end
