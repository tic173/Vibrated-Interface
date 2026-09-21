function [specs, information] = vi_wnl_apply_recovery_cache( ...
        cacheFile, repositoryRoot, input, parameters, specs, ndof)
%VI_WNL_APPLY_RECOVERY_CACHE Reuse validated direct/adjoint Floquet modes.
%
% [SPECS,INFO] = VI_WNL_APPLY_RECOVERY_CACHE(FILE,ROOT,INPUT,P,SPECS,NDOF)
% loads a previous driver output. Strictly converged direct/adjoint vectors
% bypass primitive-variable recovery. A finite matching incomplete recovery
% instead seeds a recovery-only continuation automatically.
% The current block always recomputes the applicable residual and validity
% gates before coefficients are evaluated.
%
% The cache is rejected unless the physical parameters, analysis amplitude,
% boundary condition, complete spatial/Floquet grid, mode branch, and vector
% dimensions agree.  This is a restart mechanism, not permission to move a
% mode between different operating points or discretizations.

information = struct('requested',true,'used',false,'valid',false, ...
    'warmStarted',false, ...
    'file','','reason','','sourceCodeRelease','', ...
    'sourceLinearOperatorVersion','', ...
    'sourceAdjointSolverVersion','', ...
    'numberOfModes',0,'labels',{{}});
if nargin < 1 || isempty(cacheFile)
    information.reason = 'no recovery-cache file was specified';
    return;
end
resolvedFile = resolve_cache_file(cacheFile,repositoryRoot);
information.file = resolvedFile;
if ~isfile(resolvedFile)
    information.reason = sprintf('cache file does not exist: %s', ...
        resolvedFile);
    return;
end

try
    loaded = load(resolvedFile,'output');
catch cacheError
    information.reason = sprintf('cache could not be loaded: %s', ...
        cacheError.message);
    return;
end
if ~isfield(loaded,'output') || ~isstruct(loaded.output)
    information.reason = 'cache does not contain an output structure';
    return;
end
cachedOutput = loaded.output;
if isfield(cachedOutput,'codeRelease')
    information.sourceCodeRelease = cachedOutput.codeRelease;
end
currentLinearOperatorVersion = expected_linear_operator_version(input);
if isfield(cachedOutput,'operatorMetadata') && ...
        isfield(cachedOutput.operatorMetadata,'linearOperatorVersion')
    information.sourceLinearOperatorVersion = ...
        cachedOutput.operatorMetadata.linearOperatorVersion;
end
operatorVersionMatches = strcmp( ...
    information.sourceLinearOperatorVersion,currentLinearOperatorVersion);
currentAdjointSolverVersion = expected_adjoint_solver_version();
if isfield(cachedOutput,'operatorMetadata') && ...
        isfield(cachedOutput.operatorMetadata,'adjointSolverVersion')
    information.sourceAdjointSolverVersion = ...
        cachedOutput.operatorMetadata.adjointSolverVersion;
end
adjointVersionMatches = strcmp( ...
    information.sourceAdjointSolverVersion,currentAdjointSolverVersion);

[compatible,reason] = compatible_run( ...
    cachedOutput,input,parameters,ndof);
if ~compatible
    information.reason = reason;
    return;
end

[cachedModes,reason] = extract_modes(cachedOutput);
if isempty(cachedModes)
    information.reason = reason;
    return;
end
if numel(cachedModes) < numel(specs)
    information.reason = sprintf( ...
        'cache contains %d mode(s), but this run requires %d', ...
        numel(cachedModes),numel(specs));
    return;
end
labels = cell(numel(specs),1);

strictTolerance = input.options.coefficientModeResidualTolerance;
cachedModesPassStrictGate = true(numel(specs),1);
for modeIndex = 1:numel(specs)
    cachedMode = cachedModes{modeIndex};
    currentSpec = specs{modeIndex};
    [compatible,reason] = compatible_mode( ...
        cachedMode,currentSpec,ndof);
    if ~compatible
        information.reason = sprintf('mode %d cache mismatch: %s', ...
            modeIndex,reason);
        return;
    end
    cachedModesPassStrictGate(modeIndex) = ...
        operatorVersionMatches && adjointVersionMatches && ...
        cached_mode_passes_gate(cachedMode,strictTolerance);
end

if ~all(cachedModesPassStrictGate)
    for modeIndex = 1:numel(specs)
        cachedMode = cachedModes{modeIndex};
        currentSpec = specs{modeIndex};
        hasReducedReference = isfield(currentSpec, ...
            'reducedReferenceLambda') && ...
            ~isempty(currentSpec.reducedReferenceLambda);
        if hasReducedReference
            reducedReferenceLambda = ...
                currentSpec.reducedReferenceLambda;
        end
        useExactReducedExponent = hasReducedReference && ...
            isfield(currentSpec,'exactSeparableBesselInterface') && ...
            logical(currentSpec.exactSeparableBesselInterface);
        if useExactReducedExponent
            currentSpec.lambda = reducedReferenceLambda;
        else
            currentSpec.lambda = wnl_spec_lambda(cachedMode.spec);
        end
        if hasReducedReference
            currentSpec.reducedReferenceLambda = reducedReferenceLambda;
        end
        currentSpec.directLiftInitialVector = cachedMode.vector(:);
        currentSpec.recoveryWarmStartFile = resolvedFile;
        currentSpec.recoveryCacheCodeRelease = ...
            information.sourceCodeRelease;
        specs{modeIndex} = currentSpec;
        labels{modeIndex} = currentSpec.label;
    end
    information.valid = true;
    information.warmStarted = true;
    if operatorVersionMatches && adjointVersionMatches
        information.reason = sprintf([ ...
            'matching incomplete recovery supplied as a warm start; ', ...
            'strict-ready modes %d/%d'], ...
            nnz(cachedModesPassStrictGate),numel(specs));
    elseif ~operatorVersionMatches
        sourceVersion = information.sourceLinearOperatorVersion;
        if isempty(sourceVersion)
            sourceVersion = 'legacy/unspecified';
        end
        information.reason = sprintf([ ...
            'linear operator changed from %s to %s; matching saved ', ...
            'vectors are supplied only as recovery warm starts'], ...
            sourceVersion,currentLinearOperatorVersion);
    else
        sourceVersion = information.sourceAdjointSolverVersion;
        if isempty(sourceVersion)
            sourceVersion = 'legacy/unspecified';
        end
        information.reason = sprintf([ ...
            'adjoint solver changed from %s to %s; saved direct vectors ', ...
            'are supplied only as recovery warm starts'], ...
            sourceVersion,currentAdjointSolverVersion);
    end
    information.numberOfModes = numel(specs);
    information.labels = labels;
    return;
end

for modeIndex = 1:numel(specs)
    cachedMode = cachedModes{modeIndex};
    currentSpec = specs{modeIndex};
    hasReducedReference = isfield(currentSpec, ...
        'reducedReferenceLambda') && ...
        ~isempty(currentSpec.reducedReferenceLambda);
    if hasReducedReference
        reducedReferenceLambda = currentSpec.reducedReferenceLambda;
    end
    currentSpec.direct = cachedMode.vector(:);
    currentSpec.left = cachedMode.left(:);
    currentSpec.lambda = wnl_spec_lambda(cachedMode.spec);
    if hasReducedReference
        currentSpec.reducedReferenceLambda = reducedReferenceLambda;
    end
    currentSpec.recoveryCacheFile = resolvedFile;
    currentSpec.recoveryCacheCodeRelease = ...
        information.sourceCodeRelease;
    specs{modeIndex} = currentSpec;
    labels{modeIndex} = currentSpec.label;
end

information.used = true;
information.valid = true;
information.reason = 'validated cached modes supplied to the current blocks';
information.numberOfModes = numel(specs);
information.labels = labels;
end

function version = expected_linear_operator_version(input)
sidewallCondition = 'legacyFreeSlide';
if isfield(input,'boundary') && ...
        isfield(input.boundary,'sidewallTangentialCondition')
    sidewallCondition = validatestring( ...
        input.boundary.sidewallTangentialCondition, ...
        {'legacyFreeSlide','stressFree'});
end
if strcmp(sidewallCondition,'stressFree')
    version = 'V3-pressure-compatible-stress-free-sidewall';
else
    version = 'V2-pressure-compatible-sidewall';
end
end

function version = expected_adjoint_solver_version()
version = 'V2-pressure-nullspace-bordered-temporal-schur';
end

function tf = cached_mode_passes_gate(mode,tolerance)
tf = isfinite(mode.directResidual) && isfinite(mode.leftResidual) && ...
    max(mode.directResidual,mode.leftResidual) <= tolerance;
end

function [compatible,reason] = compatible_run(output,input,parameters,ndof)
compatible = false;
reason = '';
required = {'input','parameters','operatorMetadata','weaklyNonlinear'};
for fieldIndex = 1:numel(required)
    if ~isfield(output,required{fieldIndex})
        reason = sprintf('saved output is missing field %s', ...
            required{fieldIndex});
        return;
    end
end
if ~isfield(output.operatorMetadata,'ndof') || ...
        output.operatorMetadata.ndof ~= ndof
    reason = 'unknowns per harmonic do not match';
    return;
end

parameterNames = {'omegaStar','R0','C','Bd','At','eta','g_sgn', ...
    'aAnalysis','phase'};
for fieldIndex = 1:numel(parameterNames)
    fieldName = parameterNames{fieldIndex};
    if ~isfield(output.parameters,fieldName) || ...
            ~isfield(parameters,fieldName) || ...
            ~same_numeric(output.parameters.(fieldName), ...
                parameters.(fieldName))
        reason = sprintf('physical/forcing parameter %s changed',fieldName);
        return;
    end
end

savedInput = output.input;
if ~isfield(savedInput,'numerics') || ...
        ~same_numerics(savedInput.numerics,input.numerics)
    reason = 'spatial, azimuthal, or Floquet discretization changed';
    return;
end
if ~isfield(savedInput,'boundary') || ...
        (~isequaln(savedInput.boundary,input.boundary) && ...
        ~axisymmetric_sidewall_warm_start_compatible(savedInput,input))
    reason = 'boundary/contact-line model changed';
    return;
end
if ~isfield(savedInput,'modes') || ...
        numel(savedInput.modes) ~= numel(input.modes)
    reason = 'the retained mode set changed';
    return;
end
for modeIndex = 1:numel(input.modes)
    savedMode = savedInput.modes(modeIndex);
    currentMode = input.modes(modeIndex);
    requiredModeFields = {'m','radialIndex','s'};
    for fieldIndex = 1:numel(requiredModeFields)
        fieldName = requiredModeFields{fieldIndex};
        if ~isfield(savedMode,fieldName) || ...
                ~isfield(currentMode,fieldName) || ...
                ~same_numeric(savedMode.(fieldName), ...
                    currentMode.(fieldName))
            reason = sprintf('retained mode %d field %s changed', ...
                modeIndex,fieldName);
            return;
        end
    end
end
if ~isfield(savedInput,'run') || ...
        ~isfield(savedInput.run,'operatorFactory') || ...
        ~strcmp(savedInput.run.operatorFactory,input.run.operatorFactory)
    reason = 'operator factory changed';
    return;
end
compatible = true;
end

function tf = axisymmetric_sidewall_warm_start_compatible(savedInput,input)
% The two tangential-wall formulas differ only in the u_theta row. For a
% free-line m=0 branch u_theta is identically zero, so changing only this wall
% selector does not change the physical axisymmetric mode. The caller still
% treats an operator-version mismatch as a warm start and rechecks residuals.
tf = false;
if ~isfield(savedInput,'modes') || ~isfield(input,'modes') || ...
        isempty(input.modes) || numel(savedInput.modes) ~= numel(input.modes)
    return;
end
if any([savedInput.modes.m] ~= 0) || any([input.modes.m] ~= 0)
    return;
end
savedBoundary = savedInput.boundary;
currentBoundary = input.boundary;
if ~isfield(savedBoundary,'contactLine') || ...
        ~isfield(currentBoundary,'contactLine') || ...
        ~strcmpi(savedBoundary.contactLine,'free') || ...
        ~strcmpi(currentBoundary.contactLine,'free')
    return;
end
if isfield(savedBoundary,'sidewallTangentialCondition')
    savedBoundary = rmfield(savedBoundary,'sidewallTangentialCondition');
end
if isfield(currentBoundary,'sidewallTangentialCondition')
    currentBoundary = rmfield(currentBoundary, ...
        'sidewallTangentialCondition');
end
tf = isequaln(savedBoundary,currentBoundary);
end

function [modes,reason] = extract_modes(output)
modes = {};
reason = '';
result = output.weaklyNonlinear;
if ~isstruct(result) || isempty(result)
    reason = 'saved output has no weakly nonlinear mode result';
    return;
end
if isfield(result,'modes') && iscell(result.modes)
    modes = result.modes(:);
elseif isfield(result,'mode') && isstruct(result.mode)
    modes = {result.mode};
else
    reason = 'saved WNL result contains neither mode nor modes';
end
end

function [compatible,reason] = compatible_mode(mode,spec,ndof)
compatible = false;
reason = '';
required = {'spec','vector','left','directResidual','leftResidual'};
for fieldIndex = 1:numel(required)
    if ~isfield(mode,required{fieldIndex})
        reason = sprintf('saved mode is missing field %s', ...
            required{fieldIndex});
        return;
    end
end
expectedLength = ndof*numel(spec.n);
if numel(mode.vector) ~= expectedLength || ...
        numel(mode.left) ~= expectedLength
    reason = 'direct/adjoint vector dimensions changed';
    return;
end
if mode.spec.m ~= spec.m || ...
        abs(wnl_wrap_quasifrequency(mode.spec.s)- ...
            wnl_wrap_quasifrequency(spec.s)) > 1.0e-12 || ...
        ~isequal(mode.spec.n(:),spec.n(:))
    reason = 'azimuthal or Floquet block changed';
    return;
end
if isfield(spec,'radialIndex')
    if ~isfield(mode.spec,'radialIndex') || ...
            spec.radialIndex ~= mode.spec.radialIndex
        reason = 'radial branch index changed or is absent';
        return;
    end
end
if isfield(spec,'betaStar')
    if ~isfield(mode.spec,'betaStar') || ...
            ~same_numeric(spec.betaStar,mode.spec.betaStar)
        reason = 'radial eigenvalue changed or is absent';
        return;
    end
end
if any(~isfinite(real(mode.vector))) || ...
        any(~isfinite(imag(mode.vector))) || ...
        any(~isfinite(real(mode.left))) || ...
        any(~isfinite(imag(mode.left)))
    reason = 'saved direct/adjoint vectors contain nonfinite values';
    return;
end
compatible = true;
end

function tf = same_numeric(a,b)
tf = isnumeric(a) && isnumeric(b) && isequal(size(a),size(b));
if ~tf
    return;
end
scale = max([1;abs(a(:));abs(b(:))]);
tf = all(abs(a(:)-b(:)) <= 1.0e-12*scale);
end

function tf = same_numerics(a,b)
% An empty referenceModes field is the explicit spelling of the legacy
% retained-mode default; it does not change the radial operator.
a = remove_empty_reference_modes(a);
b = remove_empty_reference_modes(b);
tf = isequaln(a,b);
end

function value = remove_empty_reference_modes(value)
if isstruct(value) && isfield(value,'radialGrid') && ...
        isstruct(value.radialGrid) && ...
        isfield(value.radialGrid,'referenceModes') && ...
        isempty(value.radialGrid.referenceModes)
    value.radialGrid = rmfield(value.radialGrid,'referenceModes');
end
if isstruct(value) && isfield(value,'radialGrid') && ...
        isstruct(value.radialGrid) && ...
        isfield(value.radialGrid,'preferredBasisLabels') && ...
        isempty(value.radialGrid.preferredBasisLabels)
    value.radialGrid = rmfield(value.radialGrid,'preferredBasisLabels');
end
end

function resolved = resolve_cache_file(fileName,repositoryRoot)
fileName = char(fileName);
isAbsolute = startsWith(fileName,filesep) || ...
    ~isempty(regexp(fileName,'^[A-Za-z]:[\\/]','once'));
[~,name,extension] = fileparts(fileName);
if isempty(extension)
    extension = '.mat';
    fileName = [fileName,extension];
end
outputName = [name,extension];
if isAbsolute
    % Relocate stale absolute paths when the repository has been moved.
    candidates = {fileName,fullfile(repositoryRoot,outputName), ...
        fullfile(repositoryRoot,'weakly_nonlinear','results',outputName), ...
        fullfile(tempdir,outputName)};
else
    candidates = {fullfile(repositoryRoot,fileName), ...
        fullfile(repositoryRoot,'weakly_nonlinear','results',outputName), ...
        fullfile(tempdir,outputName)};
end
candidates = unique(candidates,'stable');
for candidateIndex = 1:numel(candidates)
    if isfile(candidates{candidateIndex})
        resolved = candidates{candidateIndex};
        return;
    end
end
resolved = candidates{1};
end
