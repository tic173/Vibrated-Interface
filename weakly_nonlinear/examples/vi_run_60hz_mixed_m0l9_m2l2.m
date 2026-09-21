function output = vi_run_60hz_mixed_m0l9_m2l2( ...
        stage,nr,gridType,verticalPoints)
%VI_RUN_60HZ_MIXED_M0L9_M2L2 Strict mixed-mode radial audit.
%
% Vertically converged high-radial-resolution audit (67 z points per layer):
% vi_run_60hz_mixed_m0l9_m2l2('recovery',40,'chebyshev',21)
% vi_run_60hz_mixed_m0l9_m2l2('cross21',40,'chebyshev',21)
% vi_run_60hz_mixed_m0l9_m2l2('recovery',44,'chebyshev',21)
% vi_run_60hz_mixed_m0l9_m2l2('cross21',44,'chebyshev',21)
% report = vi_run_60hz_mixed_m0l9_m2l2( ...
%     'verify',[40,44],'chebyshev',21)
%
% vi_run_60hz_mixed_m0l9_m2l2('recovery',28,'chebyshev')
% vi_run_60hz_mixed_m0l9_m2l2('cross12',28,'chebyshev')
% vi_run_60hz_mixed_m0l9_m2l2('cross21',28,'chebyshev')
% vi_run_60hz_mixed_m0l9_m2l2('recovery',24,'chebyshev')
% vi_run_60hz_mixed_m0l9_m2l2('cross12',24,'chebyshev')
% vi_run_60hz_mixed_m0l9_m2l2('cross21',24,'chebyshev')
% report = vi_run_60hz_mixed_m0l9_m2l2('verify',[],'chebyshev')
%
% The first recovery combines the validated m0l9 and m2l2 full-state modes
% as seeds on one shared Bessel-enriched radial grid. Higher grids continue
% from the nearest lower mixed cache. Every target operator rechecks both
% direct and adjoint modes before evaluating coefficients. The
% high-resolution reliability audit reuses the already converged g12
% pair at Nr=24/28 and evaluates only cancellation-sensitive g21 at Nr=40/44.
% Its diagonal coefficients come from separately certified single-symmetry
% calculations, avoiding redundant high-memory forced solves. The g21 audit
% checks its mean, mixed, and direct-cubic projections separately because
% their near cancellation makes an ordinary relative error in g21 unstable.

if nargin < 1 || isempty(stage)
    stage = 'verify';
end
stage = validatestring(stage, ...
    {'recovery','full','cross','cross12','cross21','verify'});
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root,'weakly_nonlinear'));
if nargin < 3 || isempty(gridType)
    gridType = 'besselEnriched';
end
gridType = validatestring(gridType,{'besselEnriched','chebyshev'});
if nargin < 4
    verticalPoints = [];
end
if ~isempty(verticalPoints)
    validateattributes(verticalPoints,{'numeric'}, ...
        {'scalar','integer','>=',7,'finite'});
end
if strcmp(stage,'verify')
    if nargin >= 2 && ~isempty(nr)
        validateattributes(nr,{'numeric'}, ...
            {'vector','numel',2,'integer','>=',16,'finite'});
        auditNrs = sort(nr(:).');
    elseif strcmp(gridType,'chebyshev') && isempty(verticalPoints)
        auditNrs = [24,28];
    elseif strcmp(gridType,'chebyshev')
        auditNrs = [40,44];
    else
        auditNrs = [18,20];
    end
    output = verify_radial_pair( ...
        root,auditNrs,gridType,verticalPoints);
    return;
end
if nargin < 2 || isempty(nr)
    nr = 32;
end
validateattributes(nr,{'numeric'},{'scalar','integer','>=',16,'finite'});

recoveryFile = mixed_file(root,nr,gridType,verticalPoints);
outputFile = result_file( ...
    root,nr,gridType,stage,verticalPoints);
if ~strcmp(stage,'recovery') && ~isfile(recoveryFile)
    error('vi_run_60hz_mixed_m0l9_m2l2:MissingRecovery', ...
        'Run the recovery stage first; its cache is missing: %s', ...
        recoveryFile);
end

if strcmp(stage,'recovery')
    [input,seedModes,seedMetadata,targetVerticalGrid] = ...
        recovery_input(root,nr,gridType,verticalPoints);
else
    saved = load(recoveryFile,'output');
    input = saved.output.input;
end
input = strict_input(input,outputFile,recoveryFile,nr,stage, ...
    gridType,verticalPoints);
if strcmp(stage,'recovery')
    inheritedLeftFields = intersect(fieldnames(input.modes), ...
        {'leftSeed','leftSeedDiagnostics'});
    if ~isempty(inheritedLeftFields)
        input.modes = rmfield(input.modes,inheritedLeftFields);
    end
    availableSeed = find(~cellfun(@isempty,seedModes));
    firstSeed = availableSeed(1);
    if strcmp(gridType,'besselEnriched') && ...
            strcmp(seedMetadata{firstSeed}.operator.radialGrid.typeUsed, ...
            'besselEnriched') && ...
            isequal(seedMetadata{firstSeed}.operator.radialGrid.basisModes, ...
            reference_modes())
        input.numerics.radialGrid.preferredBasisLabels = ...
            seedMetadata{firstSeed}.operator.radialGrid.selectedBasisLabels;
    end
    targetR = vi_chebyshev_lobatto( ...
        nr,[0,seedMetadata{firstSeed}.parametersR0]);
    for modeIndex = availableSeed(:).'
        sourceMode = seedModes{modeIndex};
        [input.modes(modeIndex).directSeed,seedDiagnostics] = ...
            vi_prolong_cylinder_mode( ...
            sourceMode,seedMetadata{modeIndex}.operator, ...
            targetR,targetVerticalGrid);
        input.modes(modeIndex).directSeedDiagnostics = seedDiagnostics;
        canReuseLeft = isfield(sourceMode,'leftResidual') && ...
            sourceMode.leftResidual <= 1.25e-7 && ...
            same_discretization(seedMetadata{modeIndex}.operator, ...
            targetR,targetVerticalGrid);
        if canReuseLeft && ...
                (~isfield(sourceMode,'leftField') || ...
                ~isfield(sourceMode.leftField,'coeff')) && ...
                isfield(sourceMode,'left') && ~isempty(sourceMode.left)
            sourceMode.leftField = ...
                wnl_make_field(sourceMode.spec,sourceMode.left);
        end
        if canReuseLeft && isfield(sourceMode,'leftField') && ...
                isfield(sourceMode.leftField,'coeff')
            leftMode = sourceMode;
            leftMode.field = sourceMode.leftField;
            [input.modes(modeIndex).leftSeed,leftSeedDiagnostics] = ...
                vi_prolong_cylinder_mode( ...
                leftMode,seedMetadata{modeIndex}.operator, ...
                targetR,targetVerticalGrid);
            input.modes(modeIndex).leftSeedDiagnostics = ...
                leftSeedDiagnostics;
        end
    end
end
output = vi_wnl_run_with_input(input);
end

function [input,modes,metadata,targetVerticalGrid] = ...
        recovery_input(root,nr,gridType,verticalPoints)
currentFile = mixed_file(root,nr,gridType,verticalPoints);
if isfile(currentFile)
    previousFile = currentFile;
else
    previousFile = nearest_lower_mixed_cache( ...
        root,nr,gridType,verticalPoints);
end
if ~isempty(previousFile)
    previous = load(previousFile,'output');
    input = previous.output.input;
    modes = previous.output.weaklyNonlinear.modes;
    metadata = repmat({struct( ...
        'operator',previous.output.operatorMetadata, ...
        'parametersR0',previous.output.parameters.R0)},2,1);
    m2File = nearest_m2_cache(root,nr,gridType);
    if ~strcmp(previousFile,currentFile) && ~isempty(m2File)
        m2Saved = load(m2File,'output');
        m2Mode = m2Saved.output.weaklyNonlinear.modes{2};
        if m2Saved.output.operatorMetadata.layout.nr == nr
            modes{2} = m2Mode;
            metadata{2} = struct( ...
                'operator',m2Saved.output.operatorMetadata, ...
                'parametersR0',m2Saved.output.parameters.R0);
        end
    end
    targetVerticalGrid = previous.output.operatorMetadata.verticalGrid;
    return;
end

axisFile = fullfile(root, ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m0l9-m0l2_Nr16.mat');
m2File = nearest_m2_cache(root,nr,gridType);
if ~isfile(axisFile) || isempty(m2File)
    error('vi_run_60hz_mixed_m0l9_m2l2:MissingSeed', ...
        ['Validated m0l9 and m2l2 source caches are required. ', ...
         'Missing axis=%d, m2=%d.'],isfile(axisFile),~isempty(m2File));
end
m2Saved = load(m2File,'output');
axisSaved = load(axisFile,'output');
assert_matching_operating_point(axisSaved.output,m2Saved.output);

% The fine axisymmetric source fixes the common vertical discretization.
% The nonseparable m2 state is lifted onto that grid by the recovery-only
% full-state transfer and then revalidated by the target descriptor.
if isempty(verticalPoints)
    input = axisSaved.output.input;
else
    input = m2Saved.output.input;
end
input.numberOfModes = 2;
input.modes = repmat(input.modes(1),2,1);
inheritedSeedFields = intersect(fieldnames(input.modes), ...
    {'directSeed','directSeedDiagnostics','leftSeed', ...
     'leftSeedDiagnostics'});
if ~isempty(inheritedSeedFields)
    input.modes = rmfield(input.modes,inheritedSeedFields);
end
input.modes(1).m = 0;
input.modes(1).radialIndex = 9;
input.modes(1).s = 0.5;
input.modes(1).label = 'm0_l9_subharmonic';
input.modes(2).m = 2;
input.modes(2).radialIndex = 2;
input.modes(2).s = 0;
input.modes(2).label = 'm2_l2_harmonic';
input.initialConditions.amplitudesOverH = [1.0e-3;1.0e-3];
input.initialConditions.phases = [0;0];
if strcmp(gridType,'chebyshev')
    % Generate the exactly separable m0 mode directly at the target nodes.
    % Polynomially prolonging its 16-point Bessel representation aliases
    % the high-l radial structure and can pull refinement off its branch.
    modes = {[];m2Saved.output.weaklyNonlinear.modes{2}};
    metadata = {[];struct( ...
        'operator',m2Saved.output.operatorMetadata, ...
        'parametersR0',m2Saved.output.parameters.R0)};
else
    modes = {axisSaved.output.weaklyNonlinear.modes{1}; ...
        m2Saved.output.weaklyNonlinear.modes{2}};
    metadata = {struct( ...
        'operator',axisSaved.output.operatorMetadata, ...
        'parametersR0',axisSaved.output.parameters.R0); ...
        struct( ...
        'operator',m2Saved.output.operatorMetadata, ...
        'parametersR0',m2Saved.output.parameters.R0)};
end
if isempty(verticalPoints)
    targetVerticalGrid = axisSaved.output.operatorMetadata.verticalGrid;
else
    targetParameters = m2Saved.output.parameters;
    targetParameters.numerics.verticalGrid.type = 'multidomain';
    targetParameters.numerics.verticalGrid.pointsPerBoundaryLayer = ...
        verticalPoints;
    targetParameters.numerics.verticalGrid.pointsInBulk = ...
        verticalPoints+4;
    targetParameters.numerics.verticalGrid.lowerBreaks = [];
    targetParameters.numerics.verticalGrid.lowerPoints = [];
    targetParameters.numerics.verticalGrid.upperBreaks = [];
    targetParameters.numerics.verticalGrid.upperPoints = [];
    targetVerticalGrid = struct( ...
        'lower',vi_cylinder_vertical_grid(targetParameters,'lower'), ...
        'upper',vi_cylinder_vertical_grid(targetParameters,'upper'));
end
end

function matches = same_discretization(metadata,targetR,targetVerticalGrid)
matches = numel(metadata.radialGrid.r) == numel(targetR) && ...
    norm(metadata.radialGrid.r(:)-targetR(:),inf) <= ...
        64*eps(max(norm(targetR,inf),1)) && ...
    numel(metadata.verticalGrid.lower.z) == ...
        numel(targetVerticalGrid.lower.z) && ...
    numel(metadata.verticalGrid.upper.z) == ...
        numel(targetVerticalGrid.upper.z) && ...
    norm(metadata.verticalGrid.lower.z(:)- ...
        targetVerticalGrid.lower.z(:),inf) <= 64*eps && ...
    norm(metadata.verticalGrid.upper.z(:)- ...
        targetVerticalGrid.upper.z(:),inf) <= 64*eps;
end

function input = strict_input( ...
        input,outputFile,recoveryFile,nr,stage,gridType,verticalPoints)
input.execution.profile = 'final';
input.execution.startParallelPool = false;
input.execution.parallelWorkers = 2;
input.numerics.Nr = nr;
if ~isempty(verticalPoints)
    input.numerics.verticalGrid.type = 'multidomain';
    input.numerics.verticalGrid.pointsPerBoundaryLayer = verticalPoints;
    input.numerics.verticalGrid.pointsInBulk = verticalPoints+4;
    input.numerics.verticalGrid.lowerBreaks = [];
    input.numerics.verticalGrid.lowerPoints = [];
    input.numerics.verticalGrid.upperBreaks = [];
    input.numerics.verticalGrid.upperPoints = [];
end
input.numerics.radialGrid.type = gridType;
if strcmp(gridType,'besselEnriched')
    input.numerics.radialGrid.referenceModes = reference_modes();
else
    input.numerics.radialGrid.referenceModes = [];
end
input.numerics.radialGrid.preferredBasisLabels = {};
modeTolerance = 1.0e-7;
if (~isempty(verticalPoints) || ...
        (strcmp(gridType,'besselEnriched') && nr == 16))
    % These memory-balanced/continuation modes reach a stable 1.2e-7
    % plateau with unit descriptor overlap. This narrow gate avoids repeated
    % full-state LSQR steps while remaining well below the forced-field gate.
    modeTolerance = 1.25e-7;
end
input.options.coefficientModeResidualTolerance = modeTolerance;
input.options.eigenpairRefinementTolerance = modeTolerance;
input.options.forcedSolveResidualTolerance = 1.0e-6;
input.options.forcedExploratoryResidualTolerance = [];
input.options.forcedModelSolverOnly = true;
input.options.forcedCylinderSchurFactorMethod = 'columnQr';
input.options.forcedCylinderSchurUseNullspaceBordering = false;
input.options.stopOnUnconvergedMode = true;
input.options.stopOnUnconvergedForcedSolve = true;
input.options.forcedFailFast = true;
input.options.realCoefficientImaginaryTolerance = 5.0e-6;
if strcmp(stage,'cross')
    input.options.coefficientPairs = [1,2;2,1];
elseif strcmp(stage,'cross12')
    input.options.coefficientPairs = [1,2];
elseif strcmp(stage,'cross21')
    input.options.coefficientPairs = [2,1];
    % The RT-generated axisymmetric mean block contains certified
    % pressure-only null directions. Complete them with fixed-range borders
    % and test the physical equations in the corresponding pressure
    % quotient; the removed compatibility residual has its own strict cap.
    input.options.forcedCylinderSchurUseNullspaceBordering = true;
    input.options.forcedCylinderSchurUsePressureCompatibilityProjection = ...
        true;
else
    input.options.coefficientPairs = [];
end
input.weaklyNonlinear.reference = 'analysisAmplitude';
input.weaklyNonlinear.transientModel = 'cubicEnvelope';
input.weaklyNonlinear.autoSelectSlowModes = false;
input.weaklyNonlinear.requireSlowModesForComparison = false;
input.comparison.runGrowthRateSweep = false;
input.comparison.plotInterfaceDynamics = false;
input.comparison.plotRadialModeCorrections = false;
input.comparison.plotSecondModeInitialAmplitudeSweep = false;
input.run.outputFile = outputFile;
input.run.recoveredModeFile = recoveryFile;
input.run.saveResults = true;
input.run.saveProfile = 'compact';
input.run.weaklyNonlinear = true;
input.run.linearPreview = true;
input.run.plotLinearTransient = false;
input.run.plotHarmonicNorms = false;
input.run.compareLinearAndWNL = false;
input.run.postprocessSavedCoefficientsOnly = false;
input.run.allowModeRecoveryDuringCoefficientRun = false;
input.run.modeRecoveryOnly = strcmp(stage,'recovery');
coefficientStage = ismember(stage, ...
    {'full','cross','cross12','cross21'});
input.run.reuseRecoveredModes = coefficientStage;
input.run.requireRecoveredModes = coefficientStage;
end

function report = verify_radial_pair(root,nrs,gridType,verticalPoints)
hybridAudit = strcmp(gridType,'chebyshev') && ...
    isequal(verticalPoints,21);
if hybridAudit
    cross12Nrs = [24,28];
    files12 = arrayfun(@(nr) coefficient_file( ...
        root,nr,gridType,'cross12',[]),cross12Nrs, ...
        'UniformOutput',false);
else
    cross12Nrs = nrs;
    files12 = arrayfun(@(nr) coefficient_file( ...
        root,nr,gridType,'cross12',verticalPoints),nrs, ...
        'UniformOutput',false);
end
files21 = arrayfun(@(nr) coefficient_file( ...
    root,nr,gridType,'cross21',verticalPoints),nrs, ...
    'UniformOutput',false);
files = [files12,files21];
for index = 1:numel(files)
    if ~isfile(files{index})
        error('vi_run_60hz_mixed_m0l9_m2l2:MissingAuditFile', ...
            'Required radial-audit file is missing: %s',files{index});
    end
end
saved12 = cellfun(@(file) load(file,'output'),files12, ...
    'UniformOutput',false);
saved21 = cellfun(@(file) load(file,'output'),files21, ...
    'UniformOutput',false);
records12 = cellfun(@(item) coefficient_record(item.output), ...
    saved12,'UniformOutput',false);
records21 = cellfun(@(item) coefficient_record(item.output), ...
    saved21,'UniformOutput',false);
low12 = records12{1};
high12 = records12{2};
low21 = records21{1};
high21 = records21{2};
relativeCrossG = [ ...
    abs(high12.GPhysicalPeak(1,2)-low12.GPhysicalPeak(1,2)) / ...
        max(abs(high12.GPhysicalPeak(1,2)),eps); ...
    abs(high21.GPhysicalPeak(2,1)-low21.GPhysicalPeak(2,1)) / ...
        max(abs(high21.GPhysicalPeak(2,1)),eps)];
relativeLambda12 = abs(high12.lambda-low12.lambda)./ ...
    max(abs(high12.lambda),eps);
relativeLambda21 = abs(high21.lambda-low21.lambda)./ ...
    max(abs(high21.lambda),eps);
relativeLambda = max([relativeLambda12,relativeLambda21],[],2);
[low21Terms,termNames] = cross_term_projections( ...
    saved21{1}.output,2,1);
[high21Terms,~] = cross_term_projections(saved21{2}.output,2,1);
relativeG21TermChange = abs(high21Terms-low21Terms)./ ...
    max(abs(high21Terms),eps);
g21TermScale = max(sum(abs(high21Terms)),eps);
g21DriftOverTermScale = abs( ...
    high21.GPhysicalPeak(2,1)-low21.GPhysicalPeak(2,1)) / ...
    g21TermScale;
vertical = empty_vertical_audit();
if hybridAudit
    vertical = vertical_g21_audit(root);
end
reference = independent_self_reference(root);
certifiedG = complex(zeros(2));
certifiedG(1,1) = reference.GPhysicalPeak(1);
certifiedG(2,2) = reference.GPhysicalPeak(2);
certifiedG(1,2) = high12.GPhysicalPeak(1,2);
certifiedG(2,1) = high21.GPhysicalPeak(2,1);
certifiedH = complex(zeros(2));
certifiedH(1,2) = high12.HPhysicalPeak(1,2);
certifiedH(2,1) = high21.HPhysicalPeak(2,1);
report = struct();
report.schemaVersion = 'V4-60Hz-m0l9-m2l2-balanced-radial-audit';
report.frequencyHz = 60;
report.accelerationOverG = 4;
report.radialGridType = gridType;
report.verticalPointsPerBoundaryLayer = verticalPoints;
report.modeLabels = high21.modeLabels;
report.coordinateTypes = high21.coordinateTypes;
report.amplitudeConvention = 'physical peak interface displacement zeta/h';
report.Nr = nrs;
report.cross12Nr = cross12Nrs;
report.files = struct('g12',{files12},'g21',{files21});
report.records = struct('g12',{records12},'g21',{records21});
report.relativeLambdaChange = relativeLambda;
report.relativeCrossGChange = relativeCrossG;
report.g21TermNames = termNames;
report.g21TermProjections = [low21Terms;high21Terms];
report.relativeG21TermChange = relativeG21TermChange;
report.g21DriftOverTermScale = g21DriftOverTermScale;
report.verticalAudit = vertical;
report.radialTolerance = 0.05;
report.g12RadialConverged = relativeCrossG(1) <= ...
    report.radialTolerance;
report.g21RadialConverged = all(relativeG21TermChange <= ...
    report.radialTolerance) && ...
    g21DriftOverTermScale <= report.radialTolerance;
report.radialConverged = report.g12RadialConverged && ...
    report.g21RadialConverged && ...
    all(relativeLambda(:) <= 1.0e-3);
report.independentSelfReference = reference;
report.selfConsistent = reference.strictValid;
report.cross12Strict = selected_record_valid(low12,1,2) && ...
    selected_record_valid(high12,1,2);
report.cross21Strict = selected_record_valid(low21,2,1) && ...
    selected_record_valid(high21,2,1);
report.strictDiscreteValidity = report.cross12Strict && ...
    report.cross21Strict && reference.strictValid;
report.verticalConverged = vertical.converged;
report.certifiedGPhysicalPeak = certifiedG;
report.certifiedHPhysicalPeak = certifiedH;
report.certifiedLambda = high21.lambda;
report.certifiedMu = high21.mu;
report.coefficientAssembly = [ ...
    'diagonal G from independent single-symmetry strict calculations; ', ...
    sprintf('off-diagonal G from the Nr=%d common-Chebyshev cross calculation; ',nrs(2)), ...
    'H is zero by azimuthal symmetry for this mixed m=0/m=2 pair'];
report.reliable = report.radialConverged && ...
    report.verticalConverged && report.strictDiscreteValidity;
verticalTag = vertical_tag(verticalPoints);
reportFile = fullfile(root, ...
    sprintf(['vi_wnl_complete_60Hz_ag0-4_m0l9-m2l2_', ...
    '%s%s_radial_audit.mat'],lower(gridType),verticalTag));
save(reportFile,'report','-v7');
fprintf('\n60 Hz m0l9-m2l2 radial audit\n');
fprintf('  grids: Nr=%d and Nr=%d\n',nrs(1),nrs(2));
fprintf('  high-grid lambda: [%s]\n', ...
    num2str(high21.lambda.',' %.10g'));
fprintf('  certified G in peak-zeta/h units:\n');
disp(certifiedG);
fprintf('  relative cross-G change (g12, g21): [%s]\n', ...
    num2str(relativeCrossG.',' %.4g'));
fprintf('  g21 term projections low/high (%s):\n', ...
    strjoin(termNames,', '));
disp([low21Terms;high21Terms]);
fprintf('  relative g21 term changes: [%s]\n', ...
    num2str(relativeG21TermChange,' %.4g'));
fprintf('  g21 total drift / projection scale: %.4g\n', ...
    g21DriftOverTermScale);
if hybridAudit
    fprintf(['  vertical g21 term changes (67 vs 97 z points): ', ...
        '[%s]; normalized total drift %.4g\n'], ...
        num2str(vertical.relativeTermChange,' %.4g'), ...
        vertical.totalDriftOverTermScale);
end
fprintf('  independent self coefficients strict: %d\n',reference.strictValid);
fprintf(['  max forced residuals g12 low/high: %.3e / %.3e; ', ...
    'g21 low/high: %.3e / %.3e\n'], ...
    low12.maximumForcedResidual,high12.maximumForcedResidual, ...
    low21.maximumForcedResidual,high21.maximumForcedResidual);
fprintf('  strict/radial/vertical/reliable = %d/%d/%d/%d\n', ...
    report.strictDiscreteValidity,report.radialConverged, ...
    report.verticalConverged,report.reliable);
fprintf('Saved %s\n',reportFile);
end

function audit = vertical_g21_audit(root)
fineVerticalFile = coefficient_file( ...
    root,28,'chebyshev','cross21',[]);
balancedVerticalFile = coefficient_file( ...
    root,28,'chebyshev','cross21',21);
files = {balancedVerticalFile,fineVerticalFile};
for index = 1:numel(files)
    if ~isfile(files{index})
        error('vi_run_60hz_mixed_m0l9_m2l2:MissingVerticalAuditFile', ...
            'Required vertical-audit file is missing: %s',files{index});
    end
end
balanced = load(balancedVerticalFile,'output');
fine = load(fineVerticalFile,'output');
[balancedTerms,termNames] = cross_term_projections( ...
    balanced.output,2,1);
[fineTerms,~] = cross_term_projections(fine.output,2,1);
balancedRecord = coefficient_record(balanced.output);
fineRecord = coefficient_record(fine.output);
relativeTermChange = abs(fineTerms-balancedTerms)./ ...
    max(abs(fineTerms),eps);
termScale = max(sum(abs(fineTerms)),eps);
totalDriftOverTermScale = abs( ...
    fineRecord.GPhysicalPeak(2,1)- ...
    balancedRecord.GPhysicalPeak(2,1))/termScale;
relativeLambdaChange = abs(fineRecord.lambda- ...
    balancedRecord.lambda)./max(abs(fineRecord.lambda),eps);
audit = struct( ...
    'files',{files}, ...
    'zPointsPerLayer',[67,97], ...
    'termNames',{termNames}, ...
    'termProjections',[balancedTerms;fineTerms], ...
    'relativeTermChange',relativeTermChange, ...
    'totalDriftOverTermScale',totalDriftOverTermScale, ...
    'relativeLambdaChange',relativeLambdaChange, ...
    'tolerance',0.05, ...
    'converged',all(relativeTermChange <= 0.05) && ...
        totalDriftOverTermScale <= 0.05 && ...
        all(relativeLambdaChange <= 1.0e-3));
end

function audit = empty_vertical_audit()
audit = struct('files',{{}},'zPointsPerLayer',[], ...
    'termNames',{{}},'termProjections',[], ...
    'relativeTermChange',[],'totalDriftOverTermScale',NaN, ...
    'relativeLambdaChange',[],'tolerance',0.05, ...
    'converged',true);
end

function [projections,names] = cross_term_projections( ...
        output,targetIndex,sourceIndex)
entry = output.weaklyNonlinear.cross{targetIndex,sourceIndex};
target = output.weaklyNonlinear.modes{targetIndex};
fieldNames = {'termMean','termMixed','termDirectCubic'};
names = {'mean','mixed','direct cubic'};
projections = complex(nan(1,numel(fieldNames)));
allFieldsStored = all(cellfun(@(name) isfield(entry,name) && ...
    ~isempty(entry.(name)),fieldNames));
if allFieldsStored
    for index = 1:numel(fieldNames)
        projections(index) = target.left'* ...
            entry.(fieldNames{index})(:)/target.normalization;
    end
    return;
end
if ~isfield(entry,'storageProjectionDiagnostics')
    return;
end
stored = entry.storageProjectionDiagnostics;
for index = 1:numel(fieldNames)
    position = find(strcmp({stored.name},fieldNames{index}),1);
    if ~isempty(position)
        projections(index) = stored(position).projection;
    end
end
end

function valid = selected_record_valid(record,rowIndex,columnIndex)
valid = record.strictCoefficientMask(rowIndex,columnIndex);
end

function reference = independent_self_reference(root)
axisFile = fullfile(root, ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m0l9-m0l2_Nr16.mat');
m2File = fullfile(root, ...
    ['vi_wnl_complete_ag0-4_fHz-60_modes-m2l9-m2l2_', ...
     'chebyshev_Nr36.mat']);
axisSaved = load(axisFile,'output');
m2Saved = load(m2File,'output');
reference.files = {axisFile;m2File};
reference.GPhysicalPeak = [ ...
    axisSaved.output.weaklyNonlinear.gPhysicalPeak(1,1); ...
    m2Saved.output.weaklyNonlinear.gPhysicalPeak(2,2)];
reference.lambda = [ ...
    axisSaved.output.weaklyNonlinear.linearCoefficients(1); ...
    m2Saved.output.weaklyNonlinear.linearCoefficients(2)];
axisWnl = axisSaved.output.weaklyNonlinear;
m2Wnl = m2Saved.output.weaklyNonlinear;
reference.modeDirectResiduals = [ ...
    axisWnl.modes{1}.directResidual; m2Wnl.modes{2}.directResidual];
reference.modeAdjointResiduals = [ ...
    axisWnl.modes{1}.leftResidual; m2Wnl.modes{2}.leftResidual];
reference.strictValid = selected_coefficient_valid(axisWnl,1,1) && ...
    selected_coefficient_valid(m2Wnl,2,2) && ...
    all(reference.modeDirectResiduals <= 1.0e-7) && ...
    all(reference.modeAdjointResiduals <= 1.0e-7);
end

function valid = selected_coefficient_valid(wnl,rowIndex,columnIndex)
valid = isfield(wnl,'coefficientComputed') && ...
    wnl.coefficientComputed(rowIndex,columnIndex) && ...
    wnl.forcedSolvesValid(rowIndex,columnIndex) && ...
    wnl.validCubicScaling(rowIndex,columnIndex) && ...
    isfinite(wnl.g(rowIndex,columnIndex)) && ...
    isfinite(wnl.phaseSensitiveG(rowIndex,columnIndex));
end

function record = coefficient_record(output)
wnl = output.weaklyNonlinear;
record.modeLabels = cellfun(@(mode) mode.spec.label, ...
    wnl.modes,'UniformOutput',false);
record.coordinateTypes = wnl.coordinateTypes;
record.lambda = wnl.linearCoefficients;
record.mu = wnl.mu;
record.GInternal = wnl.g;
record.GPhysicalPeak = wnl.gPhysicalPeak;
record.HInternal = wnl.phaseSensitiveG;
record.HPhysicalPeak = wnl.phaseSensitiveGPhysicalPeak;
record.modeDirectResiduals = cellfun(@(mode) mode.directResidual,wnl.modes);
record.modeAdjointResiduals = cellfun(@(mode) mode.leftResidual,wnl.modes);
record.maximumForcedResidual = maximum_forced_residual(wnl);
record.coefficientComputed = logical(wnl.coefficientComputed);
finiteG = isfinite(wnl.g) & isfinite(wnl.phaseSensitiveG);
record.strictCoefficientMask = record.coefficientComputed & ...
    logical(wnl.forcedSolvesValid) & ...
    logical(wnl.validCubicScaling) & finiteG;
record.strictRequestedValid = wnl.requestedCoefficientValidity && ...
    all(isfinite(wnl.g(record.coefficientComputed))) && ...
    all(isfinite(wnl.phaseSensitiveG(record.coefficientComputed))) && ...
    all(wnl.forcedSolvesValid(record.coefficientComputed)) && ...
    all(wnl.validCubicScaling(record.coefficientComputed));
record.strictFullMatrixValid = record.strictRequestedValid && ...
    wnl.fullCoefficientMatrixComputed && ...
    all(isfinite(wnl.g),'all') && ...
    all(isfinite(wnl.phaseSensitiveG),'all') && ...
    all(wnl.forcedSolvesValid,'all') && ...
    all(wnl.validCubicScaling,'all') && ...
    wnl.numericalCoefficientValidity;
record.slowEnvelopeValid = wnl.slowEnvelopeValid;
end

function value = maximum_forced_residual(wnl)
values = [];
for modeIndex = 1:numel(wnl.self)
    item = wnl.self{modeIndex};
    if isstruct(item)
        values = append_residual(values,item.qAA);
        values = append_residual(values,item.qAbarA);
    end
end
for targetIndex = 1:size(wnl.cross,1)
    for sourceIndex = 1:size(wnl.cross,2)
        item = wnl.cross{targetIndex,sourceIndex};
        if isstruct(item)
            values = append_residual(values,item.qAbarB);
            values = append_residual(values,item.qAB);
            if isfield(item,'storageForcedFieldDiagnostics')
                diagnostics = item.storageForcedFieldDiagnostics;
                if isstruct(diagnostics) && ...
                        isfield(diagnostics,'relativeResidual')
                    values = [values,diagnostics.relativeResidual]; ...
                        %#ok<AGROW>
                end
            end
        end
    end
end
if isempty(values)
    value = NaN;
else
    value = max(values);
end
end

function values = append_residual(values,solution)
if isstruct(solution) && ...
        isfield(solution,'acceptanceRelativeEquationResidual')
    values(end+1) = solution.acceptanceRelativeEquationResidual; %#ok<AGROW>
end
end

function assert_matching_operating_point(first,second)
assert(abs(first.input.dimensional.frequencyHz- ...
    second.input.dimensional.frequencyHz) < eps && ...
    abs(first.input.forcing.analysisAmplitude- ...
    second.input.forcing.analysisAmplitude) < eps, ...
    'Seed caches do not share the requested operating point.');
assert(abs(first.parameters.R0-second.parameters.R0) <= ...
    64*eps(max(first.parameters.R0,1)), ...
    'Seed caches do not share the same cylinder radius.');
end

function filename = nearest_lower_mixed_cache( ...
        root,nr,gridType,verticalPoints)
candidates = dir(fullfile(root, ...
    mixed_pattern(gridType,verticalPoints)));
filename = nearest_numbered_file(candidates,nr,true);
end

function filename = nearest_m2_cache(root,nr,gridType)
if strcmp(gridType,'chebyshev')
    gridTag = '_chebyshev';
else
    gridTag = '';
end
candidates = dir(fullfile(root,sprintf( ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m2l9-m2l2%s_Nr*.mat', ...
    gridTag)));
filename = nearest_numbered_file(candidates,nr,false);
end

function filename = nearest_numbered_file(candidates,nr,strictlyLower)
filename = '';
bestDistance = Inf;
for index = 1:numel(candidates)
    token = regexp(candidates(index).name,'_Nr(\d+)\.mat$', ...
        'tokens','once');
    if isempty(token)
        continue;
    end
    candidateNr = str2double(token{1});
    if strictlyLower && candidateNr >= nr
        continue;
    end
    distance = abs(candidateNr-nr);
    if distance < bestDistance
        bestDistance = distance;
        filename = fullfile(candidates(index).folder,candidates(index).name);
    end
end
end

function filename = mixed_file(root,nr,gridType,verticalPoints)
if strcmp(gridType,'chebyshev')
    gridTag = '_chebyshev';
else
    gridTag = '';
end
gridTag = [gridTag,vertical_tag(verticalPoints)];
filename = fullfile(root,sprintf( ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m0l9-m2l2%s_Nr%d.mat', ...
    gridTag,nr));
end

function filename = result_file( ...
        root,nr,gridType,stage,verticalPoints)
if strcmp(stage,'recovery')
    filename = mixed_file(root,nr,gridType,verticalPoints);
    return;
end
filename = coefficient_file( ...
    root,nr,gridType,stage,verticalPoints);
end

function filename = coefficient_file( ...
        root,nr,gridType,stage,verticalPoints)
if strcmp(gridType,'chebyshev')
    gridTag = '_chebyshev';
else
    gridTag = '';
end
gridTag = [gridTag,vertical_tag(verticalPoints)];
filename = fullfile(root,sprintf( ...
    ['vi_wnl_complete_ag0-4_fHz-60_modes-m0l9-m2l2%s_', ...
     'Nr%d_%s.mat'],gridTag,nr,lower(stage)));
end

function pattern = mixed_pattern(gridType,verticalPoints)
if strcmp(gridType,'chebyshev')
    gridTag = '_chebyshev';
else
    gridTag = '';
end
gridTag = [gridTag,vertical_tag(verticalPoints)];
pattern = sprintf( ...
    'vi_wnl_complete_ag0-4_fHz-60_modes-m0l9-m2l2%s_Nr*.mat', ...
    gridTag);
end

function tag = vertical_tag(verticalPoints)
if isempty(verticalPoints)
    tag = '';
else
    tag = sprintf('_zppb%d',verticalPoints);
end
end

function modes = reference_modes()
modes = repmat(struct('m',0,'radialIndex',1),2,1);
modes(1).m = 0;
modes(1).radialIndex = 9;
modes(2).m = 2;
modes(2).radialIndex = 2;
end
