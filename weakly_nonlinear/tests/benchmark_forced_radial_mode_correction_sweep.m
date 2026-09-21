function result = benchmark_forced_radial_mode_correction_sweep( ...
        outputDirectory,userOptions)
%BENCHMARK_FORCED_RADIAL_MODE_CORRECTION_SWEEP Forced c_n(r) mode sweep.
%
% This is a branch-tracking smoke/convergence test rather than a production
% WNL coefficient run. Each requested (m,l,s) case uses the reduced Floquet
% mode only as a shift and tracking seed, solves the full stress-free-wall
% descriptor eigenproblem, and projects the recovered interface harmonics as
%
%   zeta_n(r) = b_n*J_m(beta*r) + c_n(r).
%
% The default sweep covers m=1:8 and l=1:10 at a/g=3 and 30 Hz. Harmonic
% and subharmonic reduced searches are compared first, and the dominant
% physical Floquet mode is recovered once for each spatial branch. The
% deliberately modest N/Nr/Nz values make this a practical smoke test;
% quantitative results require refinement.

if nargin < 1 || isempty(outputDirectory)
    outputDirectory = pwd;
end
if nargin < 2
    userOptions = struct();
end
options = default_options();
names = fieldnames(userOptions);
for nameIndex = 1:numel(names)
    options.(names{nameIndex}) = userOptions.(names{nameIndex});
end
validate_options(options);
if ~isfolder(outputDirectory)
    mkdir(outputDirectory);
end

parameters = physical_parameters(options);
classes = struct('name','dominantFloquet','shortName','dominant','s',NaN);
numberOfM = numel(options.mValues);
numberOfL = numel(options.lValues);
numberOfClasses = numel(classes);

result = initialize_result(parameters,options,classes);
if strlength(string(options.resumeFile)) > 0
    result = resume_result(result,options.resumeFile,parameters,options);
end
baseName = sprintf('vi_radial_correction_sweep_ag%s_%gHz', ...
    number_token(options.accelerationOverG),options.frequencyHz);
checkpointFile = fullfile(outputDirectory,[baseName,'_partial.mat']);
totalClock = tic;
numberOfCases = numberOfM*numberOfL*numberOfClasses;
caseCounter = 0;
for classIndex = 1:numberOfClasses
    classValue = classes(classIndex);
    for mPosition = 1:numberOfM
        m = options.mValues(mPosition);
        for lPosition = 1:numberOfL
            l = options.lValues(lPosition);
            caseCounter = caseCounter+1;
            caseClock = tic;
            fprintf('[%3d/%3d] m=%d, l=%d, %s: ', ...
                caseCounter,numberOfCases,m,l,classValue.shortName);
            hasRequestedProfile = ...
                ~options.retainRadialCorrectionAnalyses || ...
                ~isempty(result.radialCorrectionAnalyses{ ...
                mPosition,lPosition,classIndex});
            if result.success(mPosition,lPosition,classIndex) && ...
                    hasRequestedProfile
                if strlength(result.trackingMethod( ...
                        mPosition,lPosition,classIndex)) == 0
                    result.trackingMethod( ...
                        mPosition,lPosition,classIndex) = "direct";
                end
                fprintf('reused validated result\n');
                continue;
            end
            try
                result.failureMessage( ...
                    mPosition,lPosition,classIndex) = "";
                value = solve_case(parameters,options,m,l,classValue);
                result.success(mPosition,lPosition,classIndex) = true;
                result.betaStar(mPosition,lPosition,classIndex) = ...
                    value.betaStar;
                result.reducedLambda(mPosition,lPosition,classIndex) = ...
                    value.reducedLambda;
                result.fullLambda(mPosition,lPosition,classIndex) = ...
                    value.fullLambda;
                result.reducedResidual(mPosition,lPosition,classIndex) = ...
                    value.reducedResidual;
                result.fullResidual(mPosition,lPosition,classIndex) = ...
                    value.fullResidual;
                result.trackingOverlap(mPosition,lPosition,classIndex) = ...
                    value.trackingOverlap;
                result.continuationMinimumStepOverlap( ...
                    mPosition,lPosition,classIndex) = ...
                    value.continuationMinimumStepOverlap;
                result.trackingMethod(mPosition,lPosition,classIndex) = ...
                    string(value.trackingMethod);
                if options.retainRadialCorrectionAnalyses
                    result.radialCorrectionAnalyses{ ...
                        mPosition,lPosition,classIndex} = ...
                        value.radialCorrectionAnalysis;
                end
                result.relativeCorrectionL2( ...
                    mPosition,lPosition,classIndex) = ...
                    value.relativeCorrectionL2;
                result.peakCorrectionFraction( ...
                    mPosition,lPosition,classIndex) = ...
                    value.peakCorrectionFraction;
                result.includeInSurfacePattern( ...
                    mPosition,lPosition,classIndex) = ...
                    value.includeInSurfacePattern;
                result.radialBasisCondition( ...
                    mPosition,lPosition,classIndex) = ...
                    value.radialBasisCondition;
                result.eigsFlag(mPosition,lPosition,classIndex) = ...
                    value.eigsFlag;
                result.quasifrequency(mPosition,lPosition,classIndex) = ...
                    value.quasifrequency;
                result.sourceModeType(mPosition,lPosition,classIndex) = ...
                    string(value.sourceModeType);
                result.caseSeconds(mPosition,lPosition,classIndex) = ...
                    toc(caseClock);
                fprintf(['lambda=%+.4g%+.4gi, ||c||/||zeta||=%.3e, ', ...
                    'overlap=%.4f, include=%d, track=%s, %.2f s\n'], ...
                    real(value.fullLambda),imag(value.fullLambda), ...
                    value.relativeCorrectionL2,value.trackingOverlap, ...
                    value.includeInSurfacePattern,value.trackingMethod, ...
                    result.caseSeconds(mPosition,lPosition,classIndex));
            catch caseError
                result.failureMessage(mPosition,lPosition,classIndex) = ...
                    string(caseError.message);
                result.caseSeconds(mPosition,lPosition,classIndex) = ...
                    toc(caseClock);
                fprintf('FAILED after %.2f s: %s\n', ...
                    result.caseSeconds(mPosition,lPosition,classIndex), ...
                    caseError.message);
            end
        end
        result.elapsedSeconds = toc(totalClock);
        save(checkpointFile,'result','-v7.3');
    end
end
result.elapsedSeconds = toc(totalClock);
result.numberOfSuccessfulCases = nnz(result.success);
result.numberOfFailedCases = numberOfCases-result.numberOfSuccessfulCases;
result.table = result_table(result);

matFile = fullfile(outputDirectory,[baseName,'.mat']);
csvFile = fullfile(outputDirectory,[baseName,'.csv']);
pngFile = fullfile(outputDirectory,[baseName,'.png']);
figureFile = fullfile(outputDirectory,[baseName,'.fig']);
result.outputFiles = struct('mat',matFile,'csv',csvFile, ...
    'png',pngFile,'figure',figureFile,'checkpoint',checkpointFile);
save(matFile,'result','-v7.3');
writetable(result.table,csvFile);
figureHandle = plot_result(result,options);
exportgraphics(figureHandle,pngFile,'Resolution',220);
savefig(figureHandle,figureFile);
if ~options.figureVisible
    close(figureHandle);
end
if isfile(checkpointFile)
    delete(checkpointFile);
end
fprintf(['Forced radial-correction sweep complete: %d/%d cases, ', ...
    '%.1f s\n  %s\n  %s\n'], ...
    result.numberOfSuccessfulCases,numberOfCases,result.elapsedSeconds, ...
    matFile,pngFile);
end

function options = default_options()
options.mValues = 1:8;
options.lValues = 1:10;
options.frequencyHz = 30;
options.accelerationOverG = 3;
options.temporalCutoff = 2;
options.Nr = 16;
options.Nz = 9;
options.numberOfEigenpairs = 16;
options.continuationNumberOfEigenpairs = 24;
options.eigsTolerance = 1.0e-8;
options.eigsMaximumIterations = 500;
options.directContinuationOverlapThreshold = 0.80;
options.minimumContinuationStepOverlap = 0.10;
options.enableAccelerationContinuationFallback = true;
options.continuationAccelerationStep = 0.5;
options.radialCorrectionRelativeL2Threshold = 1.0e-2;
options.radialCorrectionTemporalSamples = 128;
options.retainRadialCorrectionAnalyses = false;
options.figureVisible = false;
options.resumeFile = "";
end

function validate_options(options)
validateattributes(options.mValues,{'numeric'}, ...
    {'vector','integer','positive','finite','nonempty'});
validateattributes(options.lValues,{'numeric'}, ...
    {'vector','integer','positive','finite','nonempty'});
validateattributes(options.frequencyHz,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(options.accelerationOverG,{'numeric'}, ...
    {'scalar','real','nonnegative','finite'});
validateattributes(options.temporalCutoff,{'numeric'}, ...
    {'scalar','integer','>=',1,'finite'});
validateattributes(options.Nr,{'numeric'}, ...
    {'scalar','integer','>=',7,'finite'});
validateattributes(options.Nz,{'numeric'}, ...
    {'scalar','integer','>=',7,'finite'});
validateattributes(options.numberOfEigenpairs,{'numeric'}, ...
    {'scalar','integer','>=',4,'finite'});
validateattributes(options.continuationNumberOfEigenpairs,{'numeric'}, ...
    {'scalar','integer','>=',4,'finite'});
validateattributes(options.eigsTolerance,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(options.eigsMaximumIterations,{'numeric'}, ...
    {'scalar','integer','positive','finite'});
validateattributes(options.directContinuationOverlapThreshold,{'numeric'}, ...
    {'scalar','real','>',0,'<=',1,'finite'});
validateattributes(options.minimumContinuationStepOverlap,{'numeric'}, ...
    {'scalar','real','>',0,'<=',1,'finite'});
validateattributes(options.enableAccelerationContinuationFallback, ...
    {'logical','numeric'},{'scalar'});
validateattributes(options.continuationAccelerationStep,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(options.radialCorrectionRelativeL2Threshold, ...
    {'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(options.radialCorrectionTemporalSamples,{'numeric'}, ...
    {'scalar','integer','>=',16,'finite'});
validateattributes(options.retainRadialCorrectionAnalyses, ...
    {'logical','numeric'},{'scalar'});
validateattributes(options.figureVisible,{'logical','numeric'},{'scalar'});
if ~(ischar(options.resumeFile) || ...
        (isstring(options.resumeFile) && isscalar(options.resumeFile)))
    error('benchmark_forced_sweep:ResumeFileType', ...
        'resumeFile must be a character vector or scalar string.');
end
end

function parameters = physical_parameters(options)
g0 = 9.81;
h = 22.0e-3;
R = 35.0e-3;
rhoUpper = 1.20;
rhoLower = 997.0;
muUpper = 1.81e-5;
muLower = 1.00e-3;
sigma = 72.0e-3;
timeScale = sqrt(h/g0);

parameters.omegaStar = 2*pi*options.frequencyHz*timeScale;
parameters.R0 = R/h;
parameters.C = (muLower/rhoLower)/sqrt(g0*h^3);
parameters.Bd = rhoLower*g0*h^2/sigma;
parameters.At = (rhoLower-rhoUpper)/(rhoLower+rhoUpper);
parameters.eta = muUpper/muLower;
parameters.g_sgn = -1;
parameters.aCritical = options.accelerationOverG;
parameters.aAnalysis = options.accelerationOverG;
parameters.phase = 0;
parameters.numerics.Nr = options.Nr;
parameters.numerics.NzUpper = options.Nz;
parameters.numerics.NzLower = options.Nz;
parameters.numerics.verticalGrid.type = 'single';
parameters.numerics.radialGrid.type = 'besselEnriched';
parameters.numerics.radialGrid.maximumProductOrder = 2;
parameters.numerics.radialGrid.maximumConditionNumber = 1.0e14;
parameters.numerics.radialGrid.independenceTolerance = 1.0e-12;
parameters.numerics.radialGrid.fallbackToChebyshev = false;
parameters.numerics.Ntheta = 12;
parameters.numerics.quadraticStep = 2.0e-4;
parameters.numerics.cubicStep = 2.0e-3;
parameters.boundary.contactLine = 'free';
parameters.boundary.sidewallTangentialCondition = 'stressFree';
parameters.dimensional = struct('g0',g0,'h',h,'R',R, ...
    'rhoUpper',rhoUpper,'rhoLower',rhoLower, ...
    'muUpper',muUpper,'muLower',muLower,'sigma',sigma, ...
    'frequencyHz',options.frequencyHz, ...
    'accelerationOverG',options.accelerationOverG);
end

function result = initialize_result(parameters,options,classes)
arraySize = [numel(options.mValues),numel(options.lValues),numel(classes)];
result = struct();
result.description = ['Coarse full-cylinder forced-Floquet radial ', ...
    'nonseparability sweep; refine N, Nr, and Nz before quantitative use.'];
result.parameters = parameters;
result.options = options;
result.mValues = options.mValues(:).';
result.lValues = options.lValues(:).';
result.classNames = {classes.name};
result.classShortNames = {classes.shortName};
result.sValues = [classes.s];
result.success = false(arraySize);
result.betaStar = nan(arraySize);
result.reducedLambda = complex(nan(arraySize));
result.fullLambda = complex(nan(arraySize));
result.reducedResidual = nan(arraySize);
result.fullResidual = nan(arraySize);
result.trackingOverlap = nan(arraySize);
result.continuationMinimumStepOverlap = nan(arraySize);
result.trackingMethod = strings(arraySize);
result.radialCorrectionAnalyses = cell(arraySize);
result.relativeCorrectionL2 = nan(arraySize);
result.peakCorrectionFraction = nan(arraySize);
result.includeInSurfacePattern = false(arraySize);
result.radialBasisCondition = nan(arraySize);
result.eigsFlag = nan(arraySize);
result.quasifrequency = nan(arraySize);
result.sourceModeType = strings(arraySize);
result.caseSeconds = nan(arraySize);
result.failureMessage = strings(arraySize);
result.elapsedSeconds = 0;
end

function result = resume_result(result,resumeFile,parameters,options)
resumeFile = char(string(resumeFile));
if ~isfile(resumeFile)
    error('benchmark_forced_sweep:MissingResumeFile', ...
        'Could not find resume file: %s',resumeFile);
end
saved = load(resumeFile,'result');
if ~isfield(saved,'result')
    error('benchmark_forced_sweep:InvalidResumeFile', ...
        'Resume file does not contain a result structure: %s',resumeFile);
end
previous = saved.result;
compatible = isequal(previous.mValues,result.mValues) && ...
    isequal(previous.lValues,result.lValues) && ...
    previous.options.temporalCutoff == options.temporalCutoff && ...
    previous.options.Nr == options.Nr && ...
    previous.options.Nz == options.Nz && ...
    abs(previous.parameters.omegaStar-parameters.omegaStar) <= ...
        1.0e-12*max(1,abs(parameters.omegaStar)) && ...
    abs(previous.parameters.aAnalysis-parameters.aAnalysis) <= ...
        1.0e-12*max(1,abs(parameters.aAnalysis));
if ~compatible
    error('benchmark_forced_sweep:IncompatibleResumeFile', ...
        'The saved sweep grid, forcing, or discretization is incompatible.');
end
copyFields = {'success','betaStar','reducedLambda','fullLambda', ...
    'reducedResidual','fullResidual','trackingOverlap', ...
    'relativeCorrectionL2','peakCorrectionFraction', ...
    'includeInSurfacePattern','radialBasisCondition','eigsFlag', ...
    'quasifrequency','sourceModeType','caseSeconds','failureMessage'};
for fieldIndex = 1:numel(copyFields)
    fieldName = copyFields{fieldIndex};
    if isfield(previous,fieldName) && ...
            isequal(size(previous.(fieldName)),size(result.(fieldName)))
        result.(fieldName) = previous.(fieldName);
    end
end
if isfield(previous,'continuationMinimumStepOverlap')
    result.continuationMinimumStepOverlap = ...
        previous.continuationMinimumStepOverlap;
end
if isfield(previous,'trackingMethod')
    result.trackingMethod = previous.trackingMethod;
else
    result.trackingMethod(result.success) = "direct";
end
if options.retainRadialCorrectionAnalyses && ...
        isfield(previous,'radialCorrectionAnalyses') && ...
        isequal(size(previous.radialCorrectionAnalyses), ...
        size(result.radialCorrectionAnalyses))
    result.radialCorrectionAnalyses = ...
        previous.radialCorrectionAnalyses;
end
weakDirect = result.success & result.trackingMethod == "direct" & ...
    result.trackingOverlap < options.directContinuationOverlapThreshold;
if any(weakDirect,'all')
    result.success(weakDirect) = false;
    result.includeInSurfacePattern(weakDirect) = false;
    result.failureMessage(weakDirect) = ...
        "Recomputing weak direct match by acceleration continuation.";
    fprintf(['Invalidated %d direct match(es) below the %.2f ', ...
        'continuation threshold.\n'],nnz(weakDirect), ...
        options.directContinuationOverlapThreshold);
end
fprintf('Resuming %d validated case(s) from %s\n', ...
    nnz(result.success),resumeFile);
end

function value = solve_case(parameters,options,m,l,classValue)
caseParameters = parameters;
roots = bessel_derivative_root(m,l);
betaStar = roots(l)/parameters.R0;
caseParameters.modes = struct('m',m,'radialIndex',l, ...
    'betaStar',betaStar);
[reduced,sourceModeType] = dominant_reduced_mode( ...
    caseParameters,options,m,l);
targetLambda = reduced.lambda;
try
    [selectedVector,selectedLambda,trackingOverlap,fullResidual, ...
        block,metadata,spec,eigsFlag] = direct_full_mode( ...
        caseParameters,options,m,l,classValue,reduced,betaStar);
    trackingMethod = 'direct';
    continuationMinimumStepOverlap = NaN;
catch directError
    canContinue = options.enableAccelerationContinuationFallback && ...
        strcmp(directError.identifier, ...
        'benchmark_forced_sweep:TrackingOverlap');
    if ~canContinue
        rethrow(directError);
    end
    [selectedVector,selectedLambda,trackingOverlap,fullResidual, ...
        block,metadata,spec,eigsFlag,continuationMinimumStepOverlap] = ...
        continued_full_mode(caseParameters,options,m,l,classValue, ...
        reduced,betaStar,sourceModeType);
    trackingMethod = 'accelerationContinuation';
end

spec.lambda = selectedLambda;
field = metadata.normalizeDirect(wnl_make_field(spec,selectedVector));
mode = struct('spec',spec,'field',field);
correctionSettings.radialCorrectionPolicy = 'auto';
correctionSettings.radialCorrectionRelativeL2Threshold = ...
    options.radialCorrectionRelativeL2Threshold;
correctionSettings.radialCorrectionTemporalSamples = ...
    options.radialCorrectionTemporalSamples;
correction = vi_radial_mode_correction( ...
    mode,metadata,caseParameters,correctionSettings);

value = struct();
value.betaStar = betaStar;
value.reducedLambda = targetLambda;
value.fullLambda = selectedLambda;
value.reducedResidual = reduced.diagnostics.vectorResidual;
value.fullResidual = fullResidual;
value.trackingOverlap = trackingOverlap;
value.continuationMinimumStepOverlap = continuationMinimumStepOverlap;
value.trackingMethod = trackingMethod;
value.radialCorrectionAnalysis = correction;
value.relativeCorrectionL2 = correction.relativeCorrectionL2;
value.peakCorrectionFraction = correction.peakCorrectionFraction;
value.includeInSurfacePattern = correction.includeInSurfacePattern;
value.radialBasisCondition = metadata.radialGrid.conditionNumber;
value.eigsFlag = eigsFlag;
value.quasifrequency = wnl_wrap_quasifrequency( ...
    reduced.s+imag(selectedLambda)/parameters.omegaStar);
value.sourceModeType = sourceModeType;
end

function [selectedVector,selectedLambda,trackingOverlap,fullResidual, ...
        block,metadata,spec,eigsFlag] = direct_full_mode( ...
        parameters,options,m,l,classValue,reduced,betaStar)
[block,metadata,spec,radialSeed] = full_block( ...
    parameters,options,m,l,classValue,reduced,betaStar);
interfaceSeed = radialSeed*reduced.zeta(:).';
initialVector = interface_seed_vector(interfaceSeed,metadata,spec);
[selectedVector,selectedLambda,trackingOverlap,fullResidual,eigsFlag] = ...
    solve_and_select(block,metadata,interfaceSeed,initialVector, ...
    reduced.lambda,options.numberOfEigenpairs, ...
    options.directContinuationOverlapThreshold,options);
end

function [selectedVector,selectedLambda,besselOverlap,fullResidual, ...
        block,metadata,spec,eigsFlag,minimumStepOverlap] = ...
        continued_full_mode(parameters,options,m,l,classValue, ...
        targetReduced,betaStar,sourceModeType)
targetAcceleration = parameters.aAnalysis;
accelerationValues = 0:options.continuationAccelerationStep: ...
    targetAcceleration;
if isempty(accelerationValues) || accelerationValues(end) < targetAcceleration
    accelerationValues(end+1) = targetAcceleration;
end
accelerationValues = unique(accelerationValues,'stable');

previousVector = [];
previousLambda = [];
previousInterface = [];
stepOverlaps = nan(size(accelerationValues));
for stepIndex = 1:numel(accelerationValues)
    stepParameters = parameters;
    stepParameters.aCritical = accelerationValues(stepIndex);
    stepParameters.aAnalysis = accelerationValues(stepIndex);
    if stepIndex == 1
        if strcmp(sourceModeType,'H')
            sectorS = 0;
        else
            sectorS = 0.5;
        end
        stepReduced = reduced_mode_in_sector( ...
            stepParameters,options,m,l,sectorS);
    else
        stepReduced = targetReduced;
    end
    [block,metadata,spec,radialSeed] = full_block( ...
        stepParameters,options,m,l,classValue,targetReduced,betaStar);
    if stepIndex == 1
        temporalSeed = align_temporal_coefficients( ...
            stepReduced.zeta,stepReduced.harmonicIndices,spec.n);
        interfaceSeed = radialSeed*temporalSeed(:).';
        initialVector = interface_seed_vector( ...
            interfaceSeed,metadata,spec);
        eigenvalueShift = stepReduced.lambda;
    else
        interfaceSeed = previousInterface;
        initialVector = previousVector/max(norm(previousVector),eps);
        eigenvalueShift = previousLambda;
    end
    [selectedVector,selectedLambda,stepOverlaps(stepIndex), ...
        fullResidual,eigsFlag] = solve_and_select( ...
        block,metadata,interfaceSeed,initialVector,eigenvalueShift, ...
        options.continuationNumberOfEigenpairs, ...
        options.minimumContinuationStepOverlap,options);
    selectedField = wnl_make_field(spec,selectedVector);
    previousInterface = selectedField.coeff(metadata.layout.zeta,:);
    previousVector = selectedVector;
    previousLambda = selectedLambda;
end
minimumStepOverlap = min(stepOverlaps,[],'omitnan');
targetSeed = radialSeed*targetReduced.zeta(:).';
besselOverlap = normalized_interface_overlap( ...
    previousInterface,targetSeed,metadata);
end

function [block,metadata,spec,radialSeed] = full_block( ...
        parameters,options,m,l,classValue,reduced,betaStar)
[operators,metadata] = cylinder_wnl_operators(parameters);
config.omega = parameters.omegaStar;
config.N = options.temporalCutoff;
config.ndof = metadata.ndof;
model = vi_wnl_model(config,operators);
spec = wnl_spec(m,reduced.s,reduced.harmonicIndices,metadata.ndof, ...
    sprintf('m%d_l%d_%s',m,l,classValue.name));
spec.betaStar = betaStar;
spec.radialIndex = l;
block = model.block(spec);
radialSeed = besselj(abs(m),betaStar*metadata.discretization.r(:));
radialSeed = radialSeed/max(abs(radialSeed));
end

function initialVector = interface_seed_vector(interfaceSeed,metadata,spec)
seedCoefficients = complex(zeros(metadata.ndof,numel(spec.n)));
if size(interfaceSeed,2) ~= numel(spec.n)
    error('benchmark_forced_sweep:ReducedHarmonics', ...
        'The reduced and full temporal harmonic counts differ.');
end
seedCoefficients(metadata.layout.zeta,:) = interfaceSeed;
initialVector = seedCoefficients(:);
initialVector = initialVector/max(norm(initialVector),eps);
end

function [selectedVector,selectedLambda,trackingOverlap,fullResidual, ...
        eigsFlag] = solve_and_select(block,metadata,interfaceSeed, ...
        initialVector,eigenvalueShift,numberOfEigenpairs, ...
        minimumOverlap,options)
eigenOptions = struct('tol',options.eigsTolerance, ...
    'maxit',options.eigsMaximumIterations,'disp',0, ...
    'v0',initialVector);
numberRequested = min(numberOfEigenpairs,size(block.A,1)-2);
[vectors,diagonal,eigsFlag] = eigs( ...
    -block.A,block.Bslow,numberRequested,eigenvalueShift,eigenOptions);
[selectedVector,selectedLambda,trackingOverlap,fullResidual] = ...
    select_target_mode(vectors,diag(diagonal),block,metadata, ...
    interfaceSeed,eigenvalueShift,minimumOverlap);
end

function aligned = align_temporal_coefficients(values,sourceN,targetN)
aligned = complex(zeros(size(targetN)));
[present,sourcePosition] = ismember(targetN,sourceN);
aligned(present) = values(sourcePosition(present));
end

function [selected,sourceModeType] = dominant_reduced_mode( ...
        parameters,options,m,l)
trialS = [0,0.5];
trialNames = {'H','SH'};
trials = cell(2,1);
failures = strings(2,1);
for trialIndex = 1:2
    try
        trials{trialIndex} = reduced_mode_in_sector( ...
            parameters,options,m,l,trialS(trialIndex));
    catch trialError
        failures(trialIndex) = string(trialError.message);
    end
end
available = find(~cellfun(@isempty,trials));
if isempty(available)
    error('benchmark_forced_sweep:ReducedRootFailure', ...
        'Both reduced root searches failed: H: %s; SH: %s', ...
        failures(1),failures(2));
end
growth = cellfun(@(value) real(value.gamma),trials(available));
maximumGrowth = max(growth);
growthTolerance = 1.0e-7*max(1,parameters.omegaStar);
nearlyDominant = available(growth >= maximumGrowth-growthTolerance);
% Prefer the harmonic truncation when both searches converge to the same
% physical multiplier. This avoids reporting a duplicate RT root shifted
% by half a forcing frequency as a separate subharmonic mode.
if any(nearlyDominant == 1)
    selectedIndex = 1;
else
    selectedIndex = nearlyDominant(1);
end
selected = trials{selectedIndex};
sourceModeType = trialNames{selectedIndex};
end

function mode = reduced_mode_in_sector(parameters,options,m,l,s)
if s == 0
    modeType = 'H';
else
    modeType = 'SH';
end
modeInput = struct('m',m,'radialIndex',l,'s',s, ...
    'label',sprintf('m%d_l%d_%s',m,l,modeType));
rootOptions = reduced_root_options(parameters.omegaStar,false);
try
    mode = vi_operating_point_floquet_mode( ...
        parameters,modeInput,options.temporalCutoff,rootOptions);
catch firstError
    rootOptions = reduced_root_options(parameters.omegaStar,true);
    try
        mode = vi_operating_point_floquet_mode( ...
            parameters,modeInput,options.temporalCutoff,rootOptions);
    catch secondError
        if strlength(string(secondError.message)) > 0
            rethrow(secondError);
        end
        rethrow(firstError);
    end
end
end

function options = reduced_root_options(omegaStar,useWideSearch)
options.verbose = false;
options.suppressSolverOutput = true;
options.residualTolerance = 1.0e-7;
options.duplicateTolerance = 1.0e-6;
options.realGuesses = omegaStar*[-1,-0.5,-0.2,0,0.2,0.5,1];
if useWideSearch
    options.imaginaryOffsets = [-0.05,0,0.05];
else
    options.imaginaryOffsets = 0;
end
end

function [vector,lambda,overlap,residual] = select_target_mode( ...
        vectors,eigenvalues,block,metadata,interfaceSeed, ...
        targetLambda,minimumOverlap)
weights = metadata.discretization.radial.quadratureWeights(:).* ...
    metadata.discretization.r(:);
seed = interfaceSeed;
seedNorm = weighted_norm(seed,weights);
numberOfCandidates = numel(eigenvalues);
overlaps = nan(numberOfCandidates,1);
distances = inf(numberOfCandidates,1);
residuals = inf(numberOfCandidates,1);
for candidateIndex = 1:numberOfCandidates
    candidateLambda = eigenvalues(candidateIndex);
    candidateVector = vectors(:,candidateIndex);
    if ~isfinite(real(candidateLambda)) || ...
            ~isfinite(imag(candidateLambda)) || ...
            any(~isfinite(real(candidateVector))) || ...
            any(~isfinite(imag(candidateVector)))
        continue;
    end
    candidateField = wnl_make_field(block.spec,candidateVector);
    zeta = candidateField.coeff(metadata.layout.zeta,:);
    zetaNorm = weighted_norm(zeta,weights);
    if zetaNorm <= 1.0e-12*max(norm(candidateVector),1)
        continue;
    end
    overlaps(candidateIndex) = abs(sum(conj(seed).*(weights.*zeta),'all')) / ...
        max(seedNorm*zetaNorm,eps);
    distances(candidateIndex) = abs(candidateLambda-targetLambda) / ...
        max([1,abs(targetLambda),abs(imag(targetLambda))]);
    residuals(candidateIndex) = norm( ...
        -block.A*candidateVector- ...
        candidateLambda*(block.Bslow*candidateVector)) / max( ...
        norm(block.A*candidateVector)+ ...
        abs(candidateLambda)*norm(block.Bslow*candidateVector),eps);
end
valid = isfinite(overlaps) & isfinite(distances) & isfinite(residuals);
if ~any(valid)
    error('benchmark_forced_sweep:NoFiniteEigenmode', ...
        'No finite full-cylinder interface eigenmode was returned.');
end
scores = -Inf(numberOfCandidates,1);
scores(valid) = overlaps(valid)./(1+distances(valid));
[~,position] = max(scores);
if overlaps(position) < minimumOverlap
    error('benchmark_forced_sweep:TrackingOverlap', ...
        ['Best full/reduced interface overlap %.3f is below the ', ...
         'required %.3f.'],overlaps(position),minimumOverlap);
end
vector = vectors(:,position);
lambda = eigenvalues(position);
overlap = overlaps(position);
residual = residuals(position);
end

function overlap = normalized_interface_overlap(candidate,seed,metadata)
weights = metadata.discretization.radial.quadratureWeights(:).* ...
    metadata.discretization.r(:);
overlap = abs(sum(conj(seed).*(weights.*candidate),'all')) / max( ...
    weighted_norm(seed,weights)*weighted_norm(candidate,weights),eps);
end

function value = weighted_norm(matrix,weights)
value = sqrt(max(real(sum(conj(matrix).*(weights.*matrix),'all')),0));
end

function tableValue = result_table(result)
numberOfRows = numel(result.success);
className = strings(numberOfRows,1);
m = zeros(numberOfRows,1);
l = zeros(numberOfRows,1);
s = zeros(numberOfRows,1);
position = 0;
for classIndex = 1:numel(result.classNames)
    for mPosition = 1:numel(result.mValues)
        for lPosition = 1:numel(result.lValues)
            position = position+1;
            className(position) = string(result.classNames{classIndex});
            m(position) = result.mValues(mPosition);
            l(position) = result.lValues(lPosition);
            s(position) = result.quasifrequency( ...
                mPosition,lPosition,classIndex);
        end
    end
end
linear = @(array) reshape(permute(array,[2,1,3]),[],1);
tableValue = table(className,m,l,s, ...
    linear(result.betaStar), ...
    real(linear(result.reducedLambda)), ...
    imag(linear(result.reducedLambda)), ...
    real(linear(result.fullLambda)), ...
    imag(linear(result.fullLambda)), ...
    linear(result.relativeCorrectionL2), ...
    linear(result.peakCorrectionFraction), ...
    linear(result.includeInSurfacePattern), ...
    linear(result.trackingOverlap), ...
    linear(result.continuationMinimumStepOverlap), ...
    linear(result.trackingMethod),linear(result.fullResidual), ...
    linear(result.radialBasisCondition),linear(result.caseSeconds), ...
    linear(result.success),linear(result.sourceModeType), ...
    linear(result.failureMessage), ...
    'VariableNames',{'class','m','l','s','betaStar', ...
    'reducedLambdaReal','reducedLambdaImag','fullLambdaReal', ...
    'fullLambdaImag','relativeCorrectionL2','peakCorrectionFraction', ...
    'included','trackingOverlap','continuationMinimumStepOverlap', ...
    'trackingMethod','fullResidual','radialBasisCondition','caseSeconds', ...
    'success','sourceModeType','failureMessage'});
end

function figureHandle = plot_result(result,options)
visibility = 'off';
if logical(options.figureVisible)
    visibility = 'on';
end
figureHandle = figure('Name','Forced radial mode-correction sweep', ...
    'Color','w','Visible',visibility,'Position',[100,100,1600,480]);
layout = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
correction = result.relativeCorrectionL2(:,:,1);
correctionPercent = 100*max(correction,1.0e-8);
included = result.includeInSurfacePattern(:,:,1);
growth = real(result.fullLambda(:,:,1))/result.parameters.omegaStar;
quasifrequency = result.quasifrequency(:,:,1);
overlap = result.trackingOverlap(:,:,1);

    correctionAxis = nexttile(layout,1);
    imagesc(correctionAxis,result.lValues,result.mValues, ...
        correctionPercent);
    axis(correctionAxis,'xy');
    set(correctionAxis,'ColorScale','log');
    colormap(correctionAxis,turbo(256));
    colorbar(correctionAxis);
    hold(correctionAxis,'on');
    [includedM,includedL] = find(included);
    plot(correctionAxis,result.lValues(includedL), ...
        result.mValues(includedM),'ws','MarkerSize',7,'LineWidth',1.2);
    title(correctionAxis,'100||c||/||\zeta|| [%]');
    xlabel(correctionAxis,'radial branch l');
    ylabel(correctionAxis,'azimuthal mode m');
    grid(correctionAxis,'on');

    growthAxis = nexttile(layout,2);
    imagesc(growthAxis,result.lValues,result.mValues,growth);
    axis(growthAxis,'xy');
    growthLimit = max(abs(growth(isfinite(growth))),[],'all');
    if isempty(growthLimit) || growthLimit <= eps
        growthLimit = 1;
    end
    clim(growthAxis,[-growthLimit,growthLimit]);
    colormap(growthAxis,blue_red_map(256));
    colorbar(growthAxis);
    title(growthAxis,'Re(\lambda)/\omega^*');
    xlabel(growthAxis,'radial branch l');
    ylabel(growthAxis,'azimuthal mode m');
    grid(growthAxis,'on');

    frequencyAxis = nexttile(layout,3);
    imagesc(frequencyAxis,result.lValues,result.mValues, ...
        quasifrequency,[0,0.5]);
    axis(frequencyAxis,'xy');
    colormap(frequencyAxis,parula(256));
    colorbar(frequencyAxis);
    title(frequencyAxis,'Floquet quasifrequency s');
    xlabel(frequencyAxis,'radial branch l');
    ylabel(frequencyAxis,'azimuthal mode m');
    grid(frequencyAxis,'on');

    overlapAxis = nexttile(layout,4);
    imagesc(overlapAxis,result.lValues,result.mValues,overlap,[0,1]);
    axis(overlapAxis,'xy');
    colormap(overlapAxis,parula(256));
    colorbar(overlapAxis);
    title(overlapAxis,'reduced/full overlap');
    xlabel(overlapAxis,'radial branch l');
    ylabel(overlapAxis,'azimuthal mode m');
    grid(overlapAxis,'on');
sgtitle(layout,sprintf([ ...
    'Stress-free wall, a/g=%.3g, f=%.3g Hz, N=%d, Nr=%d, Nz=%d; ', ...
    'white squares: c_n retained'], ...
    options.accelerationOverG,options.frequencyHz, ...
    options.temporalCutoff,options.Nr,options.Nz));
end

function map = blue_red_map(numberOfColors)
half = ceil(numberOfColors/2);
blue = [linspace(0.1,1,half).',linspace(0.25,1,half).',ones(half,1)];
red = [ones(half,1),linspace(1,0.2,half).',linspace(1,0.1,half).'];
map = [blue;red(2:end,:)];
map = interp1(linspace(0,1,size(map,1)),map, ...
    linspace(0,1,numberOfColors));
end

function token = number_token(value)
token = lower(sprintf('%.8g',value));
token = strrep(token,'+','');
token = strrep(token,'-','m');
token = strrep(token,'.','p');
end
