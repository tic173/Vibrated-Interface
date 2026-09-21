function [output, information] = vi_compact_output_record(output,profile)
%VI_COMPACT_OUTPUT_RECORD Remove reconstructible data before saving.
%
% 'compact' retains the direct/adjoint restart vectors, physical
% coefficients, primary interface fields, slaved interface fields, and
% scalar coefficient diagnostics. It removes assembled operator matrices,
% duplicate conjugate/neutral modes, and full nonlinear forcing vectors.
% 'full' leaves the result unchanged for low-level forensic debugging.

if nargin < 2 || isempty(profile)
    profile = 'compact';
end
profile = validatestring(profile,{'compact','full'});
beforeBytes = value_bytes(output);

information = base_information(profile,beforeBytes);
if strcmp(profile,'full') || ~isfield(output,'weaklyNonlinear') || ...
        ~isstruct(output.weaklyNonlinear) || ...
        isempty(output.weaklyNonlinear)
    information.estimatedBytesAfter = beforeBytes;
    information.estimatedReductionFraction = 0;
    output.storage = information;
    return;
end

wnl = output.weaklyNonlinear;
if isfield(wnl,'modes') && iscell(wnl.modes)
    for modeIndex = 1:numel(wnl.modes)
        wnl.modes{modeIndex} = compact_mode(wnl.modes{modeIndex});
    end
elseif isfield(wnl,'mode') && isstruct(wnl.mode)
    wnl.mode = compact_mode(wnl.mode);
end

duplicateFields = {'conjugateModes','neutralModes'};
for fieldIndex = 1:numel(duplicateFields)
    if isfield(wnl,duplicateFields{fieldIndex})
        wnl = rmfield(wnl,duplicateFields{fieldIndex});
    end
end

targets = retained_modes(wnl);
if isfield(wnl,'self') && iscell(wnl.self)
    for targetIndex = 1:min(numel(wnl.self),numel(targets))
        if isstruct(wnl.self{targetIndex})
            wnl.self{targetIndex} = compact_coefficient_entry( ...
                wnl.self{targetIndex},targets{targetIndex},true);
        end
    end
end
if isfield(wnl,'cross') && iscell(wnl.cross)
    for targetIndex = 1:min(size(wnl.cross,1),numel(targets))
        for sourceIndex = 1:size(wnl.cross,2)
            if isstruct(wnl.cross{targetIndex,sourceIndex}) && ...
                    ~isempty(wnl.cross{targetIndex,sourceIndex})
                wnl.cross{targetIndex,sourceIndex} = ...
                    compact_coefficient_entry( ...
                    wnl.cross{targetIndex,sourceIndex}, ...
                    targets{targetIndex},false);
            end
        end
    end
end
if isfield(wnl,'phaseSensitive') && iscell(wnl.phaseSensitive)
    for targetIndex = 1:min(size(wnl.phaseSensitive,1),numel(targets))
        for sourceIndex = 1:size(wnl.phaseSensitive,2)
            if isstruct(wnl.phaseSensitive{targetIndex,sourceIndex}) && ...
                    ~isempty(wnl.phaseSensitive{targetIndex,sourceIndex})
                wnl.phaseSensitive{targetIndex,sourceIndex} = ...
                    compact_phase_sensitive_entry( ...
                    wnl.phaseSensitive{targetIndex,sourceIndex}, ...
                    targets{targetIndex});
            end
        end
    end
end

output.weaklyNonlinear = wnl;
information.estimatedBytesAfter = value_bytes(output);
information.estimatedReductionFraction = max(0, ...
    1-information.estimatedBytesAfter/max(beforeBytes,1));
output.storage = information;
end

function information = base_information(profile,beforeBytes)
information = struct();
information.schemaVersion = 'V2-phase-sensitive-dedup';
information.profile = profile;
information.estimatedBytesBefore = beforeBytes;
information.estimatedBytesAfter = beforeBytes;
information.estimatedReductionFraction = 0;
information.capabilities = struct( ...
    'savedCoefficientPostprocessing',true, ...
    'interfaceReconstruction',true, ...
    'modeRecoveryRestart',true, ...
    'coefficientReprojection',true, ...
    'fullOperatorForensics',strcmp(profile,'full'));
end

function mode = compact_mode(mode)
if ~isstruct(mode)
    return;
end
if isfield(mode,'block') && isstruct(mode.block)
    block = mode.block;
    removable = intersect(fieldnames(block), ...
        {'A','spec','solve','solveAdjoint'});
    if ~isempty(removable)
        block = rmfield(block,removable);
    end
    mode.block = block;
end
if isfield(mode,'spec') && isstruct(mode.spec)
    mode.spec = compact_spec(mode.spec);
end
if isfield(mode,'field') && isstruct(mode.field) && ...
        isfield(mode.field,'spec')
    mode.field.spec = compact_spec(mode.field.spec);
end
if isfield(mode,'leftField')
    mode = rmfield(mode,'leftField');
end
if isfield(mode,'left') && isfield(mode,'vector') && ...
        isfield(mode,'block') && isstruct(mode.block) && ...
        isfield(mode.block,'Bslow') && isfield(mode,'normalization')
    descriptorVector = mode.block.Bslow*mode.vector(:);
    mode.storageDiagnostics = struct( ...
        'directNorm',norm(mode.vector), ...
        'leftNorm',norm(mode.left), ...
        'descriptorDirectNorm',norm(descriptorVector), ...
        'eigenpairConditionIndicator', ...
        norm(mode.left)*norm(descriptorVector)/ ...
        max(abs(mode.normalization),eps));
end
end

function spec = compact_spec(spec)
largeRestartCopies = {'directSeed','direct','left', ...
    'directLiftInitialVector'};
removable = intersect(fieldnames(spec),largeRestartCopies);
if ~isempty(removable)
    spec = rmfield(spec,removable);
end
end

function entry = compact_coefficient_entry(entry,target,isSelf)
termFields = {'termMean','termSecondHarmonic','termDifference', ...
    'termSum','termMixed','termDirectCubic','cubicForcing'};
diagnostics = repmat(struct('name','','norm',NaN,'projection',NaN),0,1);
for fieldIndex = 1:numel(termFields)
    name = termFields{fieldIndex};
    if ~isfield(entry,name) || ~isnumeric(entry.(name)) || ...
            isempty(entry.(name))
        continue;
    end
    value = entry.(name)(:);
    item = struct('name',name,'norm',norm(value),'projection',NaN);
    if isfield(target,'left') && numel(target.left) == numel(value) && ...
            isfield(target,'normalization') && ...
            abs(target.normalization) > 0
        item.projection = ...
            (target.left(:)'*value)/target.normalization;
    end
    diagnostics(end+1,1) = item; %#ok<AGROW>
end
if ~isempty(diagnostics) || ~isfield(entry,'storageProjectionDiagnostics')
    entry.storageProjectionDiagnostics = diagnostics;
end
removableTerms = intersect(fieldnames(entry),termFields);
if ~isempty(removableTerms)
    entry = rmfield(entry,removableTerms);
end

forcedFields = {'qAA','qAbarA','qBB','qBbarB','qAB','qAbarB', ...
    'qSourceSource','qTargetBarSource'};
forcedDiagnostics = repmat(struct('name','','stateNorm',NaN, ...
    'forcingNorm',NaN,'relativeResidual',NaN, ...
    'originalRelativeResidual',NaN,'valid',false),0,1);
for fieldIndex = 1:numel(forcedFields)
    name = forcedFields{fieldIndex};
    if ~isfield(entry,name) || ~isstruct(entry.(name)) || ...
            isempty(entry.(name)) || isempty(fieldnames(entry.(name)))
        continue;
    end
    solution = entry.(name);
    item = struct('name',name,'stateNorm',forced_state_norm(solution), ...
        'forcingNorm',forced_forcing_norm(solution), ...
        'relativeResidual',forced_relative_residual(solution), ...
        'originalRelativeResidual', ...
            original_forced_relative_residual(solution), ...
        'valid',logical_field(solution,'valid'));
    forcedDiagnostics(end+1,1) = item; %#ok<AGROW>
    entry.(name) = compact_forced_solution(solution);
end
entry.storageForcedFieldDiagnostics = forcedDiagnostics;

% Cross-source self fields duplicate fields already retained in self{j}.
if ~isSelf
    duplicateSourceFields = intersect(fieldnames(entry),{'qBB','qBbarB'});
    if ~isempty(duplicateSourceFields)
        entry = rmfield(entry,duplicateSourceFields);
    end
end
end

function entry = compact_phase_sensitive_entry(entry,target)
entry = compact_coefficient_entry(entry,target,false);

% Both fields are aliases of solutions already retained in self/cross.
% Their scalar residual summaries were copied above, so retaining a second
% interface-field copy adds no restart or postprocessing capability.
duplicateFields = intersect(fieldnames(entry), ...
    {'qSourceSource','qTargetBarSource'});
if ~isempty(duplicateFields)
    entry = rmfield(entry,duplicateFields);
end
end

function solution = compact_forced_solution(solution)
removable = intersect(fieldnames(solution),{'vector','forcing'});
if ~isempty(removable)
    solution = rmfield(solution,removable);
end
if isfield(solution,'field') && isstruct(solution.field) && ...
        isfield(solution.field,'spec')
    solution.field.spec = compact_spec(solution.field.spec);
end
solution.storageVectorReconstructibleFromField = ...
    isfield(solution,'field') && isstruct(solution.field) && ...
    isfield(solution.field,'coeff');
end

function value = forced_state_norm(solution)
if isfield(solution,'vector') && isnumeric(solution.vector)
    value = norm(solution.vector);
elseif isfield(solution,'field') && isstruct(solution.field) && ...
        isfield(solution.field,'coeff')
    value = norm(solution.field.coeff(:));
else
    value = NaN;
end
end

function value = forced_forcing_norm(solution)
if isfield(solution,'forcingNorm') && isnumeric(solution.forcingNorm) && ...
        isscalar(solution.forcingNorm)
    value = solution.forcingNorm;
elseif isfield(solution,'forcing') && isnumeric(solution.forcing)
    value = norm(solution.forcing);
else
    value = NaN;
end
end

function value = forced_relative_residual(solution)
names = {'acceptanceRelativeEquationResidual', ...
    'relativeEquationResidual','forcingRelativeResidual', ...
    'relativeResidual'};
value = NaN;
for fieldIndex = 1:numel(names)
    if isfield(solution,names{fieldIndex}) && ...
            isnumeric(solution.(names{fieldIndex})) && ...
            isscalar(solution.(names{fieldIndex}))
        value = solution.(names{fieldIndex});
        return;
    end
end
end

function value = original_forced_relative_residual(solution)
value = NaN;
for name = {'relativeEquationResidual','forcingRelativeResidual', ...
        'relativeResidual'}
    if isfield(solution,name{1}) && isnumeric(solution.(name{1})) && ...
            isscalar(solution.(name{1}))
        value = solution.(name{1});
        return;
    end
end
end

function value = logical_field(source,name)
value = false;
if isfield(source,name) && isscalar(source.(name))
    value = logical(source.(name));
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

function bytes = value_bytes(value)
details = whos('value');
bytes = details.bytes;
end
