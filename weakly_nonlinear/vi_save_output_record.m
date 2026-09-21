function [output, savedFile, information] = vi_save_output_record( ...
        output, requestedFile, repositoryRoot)
%VI_SAVE_OUTPUT_RECORD Save a result without relying on MATLAB's pwd.
%
% Bare filenames go to REPOSITORYROOT/weakly_nonlinear/data.
% Explicit relative subdirectories are resolved from REPOSITORYROOT. This matters when the
% driver was launched from a temporary folder that has subsequently been
% renamed or removed: SAVE with a bare filename then fails even though the
% script itself is still running. If the requested destination cannot be
% created or written, the function tries a repository data directory and
% finally MATLAB's temporary directory. The chosen absolute path is stored
% in output.save before the MAT file is written.

if nargin < 3
    error('vi_save_output_record:Inputs', ...
        'output, requestedFile, and repositoryRoot are required.');
end
if ~(ischar(requestedFile) || ...
        (isstring(requestedFile) && isscalar(requestedFile)))
    error('vi_save_output_record:RequestedFileType', ...
        'requestedFile must be a character vector or string scalar.');
end
if ~(ischar(repositoryRoot) || ...
        (isstring(repositoryRoot) && isscalar(repositoryRoot)))
    error('vi_save_output_record:RepositoryRootType', ...
        'repositoryRoot must be a character vector or string scalar.');
end
requestedFile = char(requestedFile);
repositoryRoot = char(repositoryRoot);
if isempty(strtrim(requestedFile))
    error('vi_save_output_record:EmptyRequestedFile', ...
        'The requested output filename is empty.');
end
if ~isfolder(repositoryRoot)
    error('vi_save_output_record:RepositoryRootMissing', ...
        'The repository root does not exist: %s', repositoryRoot);
end

if is_absolute_path(requestedFile)
    primaryFile = requestedFile;
    % Redirect legacy root-level output paths retained in old MAT records.
    if strcmp(fileparts(requestedFile),repositoryRoot)
        [~,name,extension]=fileparts(requestedFile);
        primaryFile=fullfile(repositoryRoot,'weakly_nonlinear','data',[name extension]);
    end
elseif isempty(fileparts(requestedFile))
    primaryFile = fullfile(repositoryRoot,'weakly_nonlinear','data',requestedFile);
else
    primaryFile = fullfile(repositoryRoot, requestedFile);
end
[~,baseName,extension] = fileparts(primaryFile);
if isempty(baseName)
    error('vi_save_output_record:MissingFilename', ...
        'The requested output path must include a filename.');
end
if isempty(extension)
    extension = '.mat';
    primaryFile = [primaryFile, extension];
end
outputName = [baseName, extension];

candidateFiles = {primaryFile, ...
    fullfile(repositoryRoot,'weakly_nonlinear','data',outputName), ...
    fullfile(tempdir,outputName)};
candidateFiles = unique_paths(candidateFiles);
errors = repmat(struct('file','','identifier','','message',''),0,1);

information = struct();
information.requestedFile = requestedFile;
information.primaryFile = primaryFile;
information.savedFile = '';
information.usedFallback = false;
information.failedCandidates = errors;

saveProfile = 'full';
if isfield(output,'input') && isstruct(output.input) && ...
        isfield(output.input,'run') && isstruct(output.input.run) && ...
        isfield(output.input.run,'saveProfile') && ...
        ~isempty(output.input.run.saveProfile)
    saveProfile = output.input.run.saveProfile;
end
[output,storageInformation] = ...
    vi_compact_output_record(output,saveProfile);
if strcmp(storageInformation.profile,'compact')
    fprintf(['Compact save profile: estimated in-memory record %.1f -> ', ...
        '%.1f MiB (%.1f%% removed before MAT compression).\n'], ...
        storageInformation.estimatedBytesBefore/2^20, ...
        storageInformation.estimatedBytesAfter/2^20, ...
        100*storageInformation.estimatedReductionFraction);
end

for candidateIndex = 1:numel(candidateFiles)
    candidateFile = candidateFiles{candidateIndex};
    candidateDirectory = fileparts(candidateFile);
    try
        if ~isfolder(candidateDirectory)
            [madeDirectory,message] = mkdir(candidateDirectory);
            if ~madeDirectory
                error('vi_save_output_record:CreateDirectory', ...
                    'Could not create %s (%s).',candidateDirectory,message);
            end
        end
        information.savedFile = candidateFile;
        information.usedFallback = candidateIndex > 1;
        information.failedCandidates = errors;
        output.save = information;
        save(candidateFile,'output','-v7.3');
        savedFile = candidateFile;
        if information.usedFallback
            fprintf(['Requested MAT-file destination was unavailable; ', ...
                'saved the reproducible record to %s instead.\n'], ...
                candidateFile);
        end
        return;
    catch saveError
        failed.file = candidateFile;
        failed.identifier = saveError.identifier;
        failed.message = saveError.message;
        errors(end+1,1) = failed; %#ok<AGROW>
    end
end

message = 'No candidate output destination could be written.';
for errorIndex = 1:numel(errors)
    message = sprintf('%s\n  %s: %s',message, ...
        errors(errorIndex).file,errors(errorIndex).message); %#ok<AGROW>
end
error('vi_save_output_record:AllDestinationsFailed','%s',message);
end

function tf = is_absolute_path(pathValue)
tf = startsWith(pathValue,'/') || startsWith(pathValue,'\') || ...
    ~isempty(regexp(pathValue,'^[A-Za-z]:[\\/]', 'once'));
end

function paths = unique_paths(paths)
keep = true(size(paths));
for pathIndex = 2:numel(paths)
    for earlierIndex = 1:pathIndex-1
        if keep(earlierIndex) && ...
                strcmpi(paths{pathIndex},paths{earlierIndex})
            keep(pathIndex) = false;
            break;
        end
    end
end
paths = paths(keep);
end
