function resolved = vi_wnl_resolve_data_file(requestedFile,repositoryRoot)
%VI_WNL_RESOLVE_DATA_FILE Find current or relocated WNL MAT-file records.
% Explicit existing paths are respected. Bare names prefer data/. Old
% absolute paths embedded in saved records can be relocated by basename.
if nargin<2 || isempty(repositoryRoot)
    repositoryRoot=fileparts(fileparts(mfilename('fullpath')));
end
requestedFile=char(requestedFile); repositoryRoot=char(repositoryRoot);
[parent,name,extension]=fileparts(requestedFile);
if isempty(extension), extension='.mat'; requestedFile=[requestedFile extension]; end
base=[name extension];
dataFile=fullfile(repositoryRoot,'weakly_nonlinear','data',base);
isAbsolute=startsWith(requestedFile,'/') || startsWith(requestedFile,'\') || ...
    ~isempty(regexp(requestedFile,'^[A-Za-z]:[\\/]','once'));
if isAbsolute
    candidates={requestedFile,dataFile};
elseif isempty(parent)
    candidates={dataFile,fullfile(repositoryRoot,requestedFile)};
else
    candidates={fullfile(repositoryRoot,requestedFile),dataFile};
end
candidates=[candidates,{fullfile(repositoryRoot,base), ...
    fullfile(repositoryRoot,'weakly_nonlinear','results',base), ...
    fullfile(tempdir,base)}];
for i=1:numel(candidates)
    if isfile(candidates{i}), resolved=candidates{i}; return; end
end
resolved=candidates{1};
end
