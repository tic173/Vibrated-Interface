function result = vi_cylinder_wnl_forced_schur_solve( ...
        blocks,omega,spec,forcing,userOpts)
%VI_CYLINDER_WNL_FORCED_SCHUR_SOLVE Exact low-rank temporal block solve.
%
% The vertically vibrated cylinder couples neighboring Floquet harmonics
% only through the interface-displacement columns of Lplus and Lminus. If
% D_j denotes one uncoupled primitive-variable temporal block, then
%
%   D_j q_j + A_+ q_{j-1} + A_- q_{j+1} = f_j,
%
% where A_+ and A_- have only a few nonzero columns. Eliminating each full
% q_j therefore produces a small block-tridiagonal system for those active
% interface coordinates. This routine solves that Schur system and then
% reconstructs velocity, pressure, and interface displacement. It is
% algebraically equivalent to solving the complete Floquet matrix; no row
% is projected, reweighted, or removed.

if nargin < 5
    userOpts = struct();
end
opts = wnl_options(userOpts);
validateattributes(opts.forcedUseCylinderSchur,{'logical','numeric'}, ...
    {'scalar'});
validateattributes(opts.forcedCylinderSchurUseEquilibratedBlocks, ...
    {'logical','numeric'},{'scalar'});
validatestring(opts.forcedCylinderSchurFactorMethod, ...
    {'equilibratedLu','columnQr'});
validateattributes(opts.forcedCylinderSchurBlockRefinementSteps, ...
    {'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(opts.forcedCylinderSchurReducedRefinementSteps, ...
    {'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(opts.forcedCylinderSchurBlockRegularization, ...
    {'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(opts.forcedCylinderSchurBlockRegularizationGrowth, ...
    {'numeric'},{'scalar','real','>',1,'finite'});
validateattributes(opts.forcedCylinderSchurBlockRegularizationAttempts, ...
    {'numeric'},{'scalar','integer','nonnegative','finite'});
validateattributes(opts.forcedCylinderSchurBlockMaximumInverseGain, ...
    {'numeric'},{'scalar','real','positive','finite'});
validateattributes(opts.forcedCylinderSchurUseNullspaceBordering, ...
    {'logical','numeric'},{'scalar'});
validateattributes(opts.forcedCylinderSchurNullTolerance, ...
    {'numeric'},{'scalar','real','positive','finite'});
validateattributes(opts.forcedCylinderSchurMaximumNullity, ...
    {'numeric'},{'scalar','integer','positive','finite'});
validateattributes( ...
    opts.forcedCylinderSchurUsePressureCompatibilityProjection, ...
    {'logical','numeric'},{'scalar'});
validateattributes( ...
    opts.forcedCylinderSchurPressureCompatibilityTolerance, ...
    {'numeric'},{'scalar','real','positive','finite'});
forcing = forcing(:);
numberOfHarmonics = numel(spec.n);
ndof = size(blocks.B0,1);
expectedLength = ndof*numberOfHarmonics;
if numel(forcing) ~= expectedLength
    error('vi_cylinder_wnl_forced_schur_solve:ForcingSize', ...
        'Forcing has %d entries; the requested block needs %d.', ...
        numel(forcing),expectedLength);
end

information = initial_information();
information.attempted = true;
information.numberOfHarmonics = numberOfHarmonics;
result = struct('vector',complex(zeros(expectedLength,1)), ...
    'diagnostics',information);
if ~opts.forcedUseCylinderSchur
    result.diagnostics.stopReason = ...
        'the cylinder temporal-Schur solver was disabled';
    return;
end

solverClock = tic;
try
    validate_blocks(blocks,ndof);
    previousCoupling = -blocks.Lplus;
    nextCoupling = -blocks.Lminus;
    activeColumns = active_coupling_columns( ...
        previousCoupling,nextCoupling);
    numberOfActiveColumns = numel(activeColumns);
    information.numberOfActiveColumns = numberOfActiveColumns;

    lambda = wnl_spec_lambda(spec);
    [~,referenceIndex] = min(abs(spec.n(:)+spec.s));
    referenceDiagonal = temporal_diagonal( ...
        blocks,omega,spec,lambda,referenceIndex);
    pressureBasis = pressure_nullspace_basis( ...
        referenceDiagonal,blocks.B0,spec,opts);
    information.pressureNullspaceDetected = pressureBasis.nullity > 0;
    information.pressureNullspaceCertified = pressureBasis.certified;
    information.pressureNullspaceNullity = pressureBasis.nullity;
    information.pressureNullspaceReferenceHarmonic = ...
        spec.n(referenceIndex)+spec.s;
    information.pressureNullspaceRightResidual = ...
        pressureBasis.rightResidual;
    information.pressureNullspaceMassResidual = pressureBasis.massResidual;
    information.pressureNullspaceLeftResidual = ...
        pressureBasis.leftResidual;
    information.pressureNullspaceLeftMassResidual = ...
        pressureBasis.leftMassResidual;
    diagonals = cell(numberOfHarmonics,1);
    previousResponse = cell(numberOfHarmonics,1);
    nextResponse = cell(numberOfHarmonics,1);
    maximumResponseResidual = 0;
    blockRegularization = zeros(numberOfHarmonics,1);
    blockInverseGain = nan(numberOfHarmonics,1);
    blockNullity = zeros(numberOfHarmonics,1);
    blockFactorMethod = cell(numberOfHarmonics,1);
    leftCompatibilityBasis = cell(numberOfHarmonics,1);
    nullspaceBasis = [];
    if pressureBasis.certified
        nullspaceBasis = pressureBasis;
    end
    reducedRegularization = NaN;
    reducedInverseGain = NaN;
    forcingByHarmonic = reshape(forcing,ndof,numberOfHarmonics);
    particular = complex(zeros(ndof,numberOfHarmonics));

    factorClock = tic;
    for harmonicIndex = 1:numberOfHarmonics
        diagonal = temporal_diagonal( ...
            blocks,omega,spec,lambda,harmonicIndex);
        diagonals{harmonicIndex} = diagonal;
        [diagonalFactor, ...
            blockRegularization(harmonicIndex), ...
            blockInverseGain(harmonicIndex)] = ...
            factor_temporal_system(diagonal,opts,sprintf( ...
            'temporal block %d',harmonicIndex),nullspaceBasis);
        blockNullity(harmonicIndex) = ...
            diagonalFactor.nullity;
        blockFactorMethod{harmonicIndex} = ...
            diagonalFactor.method;
        leftCompatibilityBasis{harmonicIndex} = ...
            harmonic_left_compatibility_basis( ...
            diagonalFactor,diagonal,opts);
        if numberOfActiveColumns > 0
            responseRightHandSide = [ ...
                previousCoupling(:,activeColumns), ...
                nextCoupling(:,activeColumns), ...
                forcingByHarmonic(:,harmonicIndex)];
            response = refined_factor_solve( ...
                diagonalFactor,diagonal, ...
                responseRightHandSide, ...
                opts.forcedCylinderSchurBlockRefinementSteps);
            previousResponse{harmonicIndex} = ...
                response(:,1:numberOfActiveColumns);
            nextResponse{harmonicIndex} = ...
                response(:,numberOfActiveColumns+1: ...
                2*numberOfActiveColumns);
            particular(:,harmonicIndex) = response(:,end);
            responseResidual = norm( ...
                diagonal*response(:,1:2*numberOfActiveColumns)- ...
                responseRightHandSide(:,1:2*numberOfActiveColumns), ...
                'fro') / max(norm(responseRightHandSide( ...
                :,1:2*numberOfActiveColumns),'fro'),eps);
            maximumResponseResidual = max( ...
                maximumResponseResidual,responseResidual);
        else
            previousResponse{harmonicIndex} = ...
                complex(zeros(ndof,0));
            nextResponse{harmonicIndex} = complex(zeros(ndof,0));
            particular(:,harmonicIndex) = refined_factor_solve( ...
                diagonalFactor,diagonal, ...
                forcingByHarmonic(:,harmonicIndex), ...
                opts.forcedCylinderSchurBlockRefinementSteps);
        end
        clear diagonalFactor;
    end
    information.factorSeconds = toc(factorClock);
    information.maximumResponseResidual = maximumResponseResidual;
    information.blockRegularizationRange = ...
        finite_range(blockRegularization);
    information.maximumBlockInverseGain = maximum_finite(blockInverseGain);
    information.blockNullityRange = finite_range(blockNullity);
    information.blockFactorMethods = unique(blockFactorMethod,'stable');

    information.particularSeconds = 0;

    reducedClock = tic;
    if numberOfActiveColumns == 0
        stateByHarmonic = particular;
        information.reducedDimension = 0;
        information.reducedResidual = 0;
    else
        [reduced,reducedRightHandSide,previousPositions, ...
            nextPositions] = assemble_reduced_system( ...
            spec,activeColumns,particular,previousResponse,nextResponse);
        information.reducedDimension = size(reduced,1);
        [reducedFactor,reducedRegularization,reducedInverseGain] = ...
            factor_temporal_system(reduced,opts, ...
            'reduced temporal Schur system',[]);
        activeState = refined_factor_solve( ...
            reducedFactor,reduced,reducedRightHandSide, ...
            opts.forcedCylinderSchurReducedRefinementSteps);
        information.reducedResidual = norm( ...
            reduced*activeState-reducedRightHandSide) / ...
            max(norm(reducedRightHandSide),eps);
        activeState = reshape(activeState,numberOfActiveColumns, ...
            numberOfHarmonics);
        stateByHarmonic = particular;
        for harmonicIndex = 1:numberOfHarmonics
            previousPosition = previousPositions(harmonicIndex);
            nextPosition = nextPositions(harmonicIndex);
            if previousPosition > 0
                stateByHarmonic(:,harmonicIndex) = ...
                    stateByHarmonic(:,harmonicIndex) - ...
                    previousResponse{harmonicIndex} * ...
                    activeState(:,previousPosition);
            end
            if nextPosition > 0
                stateByHarmonic(:,harmonicIndex) = ...
                    stateByHarmonic(:,harmonicIndex) - ...
                    nextResponse{harmonicIndex} * ...
                    activeState(:,nextPosition);
            end
        end
    end
    information.reducedRegularization = reducedRegularization;
    information.reducedInverseGain = reducedInverseGain;
    information.reducedSeconds = toc(reducedClock);

    % Opt-in physical conjugacy for a real axisymmetric forced field.
    % Validate the projected state with the same raw/quotient residual gates
    % below, rather than modifying an already-accepted state afterwards.
    if isfield(opts,'forcedCylinderSchurEnforceRealField') && ...
            opts.forcedCylinderSchurEnforceRealField
        assert(spec.m==0 && abs(imag(lambda))<1e-12 && ...
            abs(2*spec.s-round(2*spec.s))<1e-12, ...
            'Real-field enforcement requires a self-conjugate block.');
        [found,partners]=ismember(-spec.n-2*spec.s,spec.n);
        assert(all(found),'The harmonic window must be conjugacy closed.');
        forcingDefect=norm(forcingByHarmonic- ...
            conj(forcingByHarmonic(:,partners)),'fro')/max(norm(forcing),eps);
        assert(forcingDefect<1e-10,'The supplied forcing is not self-conjugate.');
        reflected=conj(stateByHarmonic(:,partners));
        information.realFieldRawConjugacyDefect= ...
            norm(stateByHarmonic-reflected,'fro')/ ...
            max(norm(stateByHarmonic,'fro'),eps);
        stateByHarmonic=0.5*(stateByHarmonic+reflected);
        information.realFieldProjectionApplied=true;
    end

    vector = stateByHarmonic(:);
    residual = complete_temporal_residual( ...
        blocks,omega,spec,lambda,stateByHarmonic,forcingByHarmonic);
    information.relativeResidual = norm(residual(:)) / ...
        max(norm(forcing),eps);
    quotient = pressure_compatible_residual( ...
        residual,diagonals,pressureBasis,leftCompatibilityBasis, ...
        norm(forcing),opts);
    information.pressureCompatibilityProjectionApplied = ...
        quotient.applied;
    information.pressureCompatibleRelativeResidual = ...
        quotient.relativeResidual;
    information.pressureCompatibilityRelativeResidual = ...
        quotient.compatibilityRelativeResidual;
    information.pressureCompatibilityAnalyzedHarmonics = ...
        quotient.analyzedHarmonics;
    information.pressureCompatibilityMaximumNullResidual = ...
        quotient.maximumNullResidual;
    information.pressureCompatibilityGatePassed = ...
        ~quotient.applied || ...
        quotient.compatibilityRelativeResidual <= ...
        opts.forcedCylinderSchurPressureCompatibilityTolerance;
    information.acceptanceRelativeResidual = quotient.relativeResidual;
    information.available = all(isfinite(real(vector))) && ...
        all(isfinite(imag(vector))) && ...
        isfinite(information.relativeResidual) && ...
        isfinite(information.acceptanceRelativeResidual);
    information.passedPhysicalGate = information.available && ...
        information.acceptanceRelativeResidual <= ...
        opts.forcedSolveResidualTolerance && ...
        information.pressureCompatibilityGatePassed;
    information.totalSeconds = toc(solverClock);
    if information.passedPhysicalGate && quotient.applied
        information.stopReason = ...
            'the certified pressure-quotient Floquet equations passed';
    elseif information.passedPhysicalGate
        information.stopReason = ...
            'the original unequilibrated Floquet equations passed';
    elseif information.available
        information.stopReason = ...
            'the Schur reconstruction missed the physical residual gate';
    else
        information.stopReason = ...
            'the Schur reconstruction returned a nonfinite state';
    end
    result.vector = vector;
    result.diagnostics = information;
catch solverError
    information.errorIdentifier = solverError.identifier;
    information.errorMessage = solverError.message;
    information.stopReason = sprintf( ...
        'cylinder temporal-Schur solve unavailable: %s', ...
        solverError.message);
    information.totalSeconds = toc(solverClock);
    result.diagnostics = information;
end
end

function information = initial_information()
information = struct('attempted',false,'available',false, ...
    'method','cylinder low-rank temporal Schur', ...
    'passedPhysicalGate',false,'relativeResidual',Inf, ...
    'acceptanceRelativeResidual',Inf, ...
    'pressureCompatibleRelativeResidual',Inf, ...
    'pressureCompatibilityRelativeResidual',0, ...
    'pressureCompatibilityProjectionApplied',false, ...
    'pressureCompatibilityGatePassed',true, ...
    'pressureCompatibilityAnalyzedHarmonics',0, ...
    'pressureCompatibilityMaximumNullResidual',NaN, ...
    'pressureNullspaceDetected',false, ...
    'pressureNullspaceCertified',false, ...
    'pressureNullspaceNullity',0, ...
    'pressureNullspaceReferenceHarmonic',NaN, ...
    'pressureNullspaceRightResidual',NaN, ...
    'pressureNullspaceMassResidual',NaN, ...
    'pressureNullspaceLeftResidual',NaN, ...
    'pressureNullspaceLeftMassResidual',NaN, ...
    'reducedResidual',NaN,'maximumResponseResidual',NaN, ...
    'blockRegularizationRange',[NaN,NaN], ...
    'blockNullityRange',[NaN,NaN],'blockFactorMethods',{{}}, ...
    'maximumBlockInverseGain',NaN, ...
    'reducedRegularization',NaN,'reducedInverseGain',NaN, ...
    'numberOfHarmonics',0,'numberOfActiveColumns',0, ...
    'reducedDimension',0,'factorSeconds',0,'particularSeconds',0, ...
    'reducedSeconds',0,'totalSeconds',0,'stopReason','not attempted', ...
    'errorIdentifier','','errorMessage','');
end

function validate_blocks(blocks,ndof)
required = {'B0','L0','Lplus','Lminus'};
for fieldIndex = 1:numel(required)
    fieldName = required{fieldIndex};
    if ~isfield(blocks,fieldName) || ...
            ~isequal(size(blocks.(fieldName)),[ndof,ndof])
        error('vi_cylinder_wnl_forced_schur_solve:BlockSize', ...
            'blocks.%s is missing or has the wrong size.',fieldName);
    end
end
end

function columns = active_coupling_columns(previousCoupling,nextCoupling)
[~,previousColumns] = find(previousCoupling);
[~,nextColumns] = find(nextCoupling);
columns = unique([previousColumns;nextColumns]);
columns = columns(:).';
end

function diagonal = temporal_diagonal( ...
        blocks,omega,spec,lambda,harmonicIndex)
frequency = spec.n(harmonicIndex)+spec.s;
diagonal = (lambda+1i*omega*frequency)*blocks.B0-blocks.L0;
end

function [factor,regularization,inverseGain] = ...
        factor_temporal_system(matrix,opts,description,nullspaceBasis)
if nargin < 4
    nullspaceBasis = [];
end
nullityHint = [];
if isstruct(nullspaceBasis) && isfield(nullspaceBasis,'nullity')
    nullityHint = nullspaceBasis.nullity;
elseif isnumeric(nullspaceBasis)
    nullityHint = nullspaceBasis;
end
matrixSize = size(matrix,1);
if size(matrix,2) ~= matrixSize
    error('vi_cylinder_wnl_forced_schur_solve:FactorSize', ...
        'The %s matrix is not square.',description);
end
identity = speye(matrixSize);
oneNorm = norm(matrix,1);
infinityNorm = norm(matrix,inf);
blockScale = max(sqrt(max(oneNorm,0)*max(infinityNorm,0)),1.0);
baseShift = opts.forcedCylinderSchurBlockRegularization*blockScale;
if baseShift > 0 && opts.forcedCylinderSchurBlockRegularizationAttempts > 0
    shiftGrowth = opts.forcedCylinderSchurBlockRegularizationGrowth;
    shiftAttempts = opts.forcedCylinderSchurBlockRegularizationAttempts;
    shifts = [0, baseShift*shiftGrowth.^(0:shiftAttempts-1)];
else
    shifts = 0;
end
probeIndex = (1:matrixSize).';
probe = [complex(ones(matrixSize,1)), ...
    exp(1i*sqrt(2)*probeIndex)];
probeNorm = sqrt(sum(abs(probe).^2,1));
factor = [];
regularization = NaN;
inverseGain = Inf;
lastMessage = 'no factor attempt was made';

% A certified pressure-nullspace basis takes precedence over structural
% rank. Spectral collocation can turn an exact pressure gauge freedom into
% a tiny nonzero pivot, so sprank(matrix)==matrixSize does not imply that
% an unconstrained LU is physically meaningful. Complete that known
% nullspace before solving any forcing or interface-response columns.
hasCertifiedNullityHint = ...
    opts.forcedCylinderSchurUseNullspaceBordering && ...
    isstruct(nullspaceBasis) && ...
    isfield(nullspaceBasis,'certified') && nullspaceBasis.certified && ...
    nullityHint > 0 && ...
    ~contains(description,'reduced temporal Schur system');
if hasCertifiedNullityHint
    try
        [trialFactor,trialGain,relativeTrialResidual] = ...
            construct_fixed_range_border_factor( ...
            matrix,probe,probeNorm,opts,nullspaceBasis,description);
        if trialGain <= opts.forcedCylinderSchurBlockMaximumInverseGain && ...
                relativeTrialResidual <= 1.0e-6
            factor = trialFactor;
            regularization = 0;
            inverseGain = trialGain;
            return;
        end
        lastMessage = sprintf([ ...
            'hinted nullspace border gave inverse gain %.3e and ', ...
            'augmented residual %.3e'],trialGain,relativeTrialResidual);
    catch factorError
        lastMessage = factorError.message;
    end
end

% First try the exact unshifted equations. If their primitive-variable
% pressure block is rank deficient, complete only that nullspace with a
% bordered minimum-norm constraint. A uniform diagonal shift perturbs every
% momentum, incompressibility, boundary, and interface equation and is kept
% solely as the final seed fallback.
useColumnQr = strcmp(opts.forcedCylinderSchurFactorMethod,'columnQr') && ...
    ~contains(description,'reduced temporal Schur system');
if useColumnQr
    try
        factor = construct_column_equilibrated_qr_factor( ...
            matrix,nullityHint);
        regularization = 0;
        inverseGain = NaN;
        return;
    catch factorError
        lastMessage = factorError.message;
    end
end
if opts.forcedCylinderSchurUseEquilibratedBlocks
    if sprank(matrix) == matrixSize
        try
            trialFactor = construct_equilibrated_lu_factor(matrix);
            factor = trialFactor;
            regularization = 0;
            inverseGain = NaN;
            return;
        catch factorError
            lastMessage = factorError.message;
        end
    else
        lastMessage = ...
            'the unshifted matrix is structurally rank deficient';
    end
end
try
    [trialFactor,trialGain,relativeTrialResidual] = ...
        construct_lu_factor(matrix,probe,probeNorm);
    if trialGain <= opts.forcedCylinderSchurBlockMaximumInverseGain && ...
            relativeTrialResidual <= 1.0e-6
        factor = trialFactor;
        regularization = 0;
        inverseGain = trialGain;
        return;
    end
    lastMessage = sprintf( ...
        'unshifted LU gave inverse gain %.3e and residual %.3e', ...
        trialGain,relativeTrialResidual);
catch factorError
    lastMessage = factorError.message;
end

if opts.forcedCylinderSchurUseNullspaceBordering
    try
        [trialFactor,trialGain,relativeTrialResidual] = ...
            construct_nullspace_bordered_factor( ...
            matrix,probe,probeNorm,opts,nullityHint,description);
        if trialGain <= opts.forcedCylinderSchurBlockMaximumInverseGain && ...
                relativeTrialResidual <= 1.0e-6
            factor = trialFactor;
            regularization = 0;
            inverseGain = trialGain;
            return;
        end
        lastMessage = sprintf([ ...
            'nullspace-bordered LU gave inverse gain %.3e and augmented ', ...
            'residual %.3e'],trialGain,relativeTrialResidual);
    catch factorError
        lastMessage = factorError.message;
    end
end

for shiftIndex = find(shifts > 0)
    shift = shifts(shiftIndex);
    shiftedMatrix = matrix;
    shiftedMatrix = shiftedMatrix+shift*identity;
    try
        [trialFactor,trialGain,relativeTrialResidual] = ...
            construct_lu_factor(shiftedMatrix,probe,probeNorm);
        if trialGain <= ...
                opts.forcedCylinderSchurBlockMaximumInverseGain && ...
                relativeTrialResidual <= 1.0e-6
            factor = trialFactor;
            regularization = shift;
            inverseGain = trialGain;
            return;
        end
        lastMessage = sprintf( ...
            'shift %.3e gave inverse gain %.3e and residual %.3e', ...
            shift,trialGain,relativeTrialResidual);
    catch factorError
        lastMessage = factorError.message;
    end
end
error('vi_cylinder_wnl_forced_schur_solve:TemporalFactorization', ...
    ['Could not construct a finite factor for %s after %d ', ...
     'attempt(s): %s'],description,numel(shifts),lastMessage);
end

function factor = construct_column_equilibrated_qr_factor( ...
        matrix,nullityHint)
if isempty(nullityHint)
    nullityHint = 0;
end
matrixSize = size(matrix,1);
columnNorm = full(sqrt(sum(abs(matrix).^2,1))).';
largest = max(columnNorm);
columnScale = ones(size(columnNorm));
if ~isempty(largest) && largest > eps
    floorValue = max(largest*1.0e-14,eps);
    positive = columnNorm > 0;
    columnScale(positive) = ...
        1./max(columnNorm(positive),floorValue);
end
scaledMatrix = matrix*spdiags( ...
    columnScale,0,matrixSize,matrixSize);
factor = struct('method','column-equilibrated-qr', ...
    'decomposition',decomposition(scaledMatrix,'qr'), ...
    'numberOfStateRows',matrixSize,'nullity',nullityHint, ...
    'columnScale',columnScale);
end

function factor = construct_equilibrated_lu_factor(matrix)
matrixSize = size(matrix,1);
rowNorm = full(sqrt(sum(abs(matrix).^2,2)));
rowScale = 1./max(rowNorm,eps);
rowScaled = spdiags(rowScale,0,matrixSize,matrixSize)*matrix;
columnNorm = full(sqrt(sum(abs(rowScaled).^2,1))).';
columnScale = 1./max(columnNorm,eps);
scaledMatrix = rowScaled*spdiags( ...
    columnScale,0,matrixSize,matrixSize);
warningStates = suppress_factor_warnings();
cleanup = onCleanup(@() restore_warning_states(warningStates));
scaledDecomposition = decomposition(scaledMatrix,'lu');
clear cleanup;
factor = struct('method','equilibrated-lu','decomposition', ...
    scaledDecomposition, ...
    'numberOfStateRows',size(matrix,1),'nullity',0, ...
    'rowScale',rowScale,'columnScale',columnScale);
end

function [factor,inverseGain,relativeResidual] = ...
        construct_lu_factor(matrix,probe,probeNorm)
factor = struct('method','lu','decomposition', ...
    decomposition(matrix,'lu'),'numberOfStateRows',size(matrix,1), ...
    'nullity',0);
trialSolution = factor_solve_safely(factor,probe);
solutionNorm = sqrt(sum(abs(trialSolution).^2,1));
inverseGain = max(solutionNorm./probeNorm);
trialResidual = matrix*trialSolution-probe;
relativeResidual = max( ...
    sqrt(sum(abs(trialResidual).^2,1))./probeNorm);
if any(~isfinite(real(trialSolution(:)))) || ...
        any(~isfinite(imag(trialSolution(:))))
    inverseGain = Inf;
    relativeResidual = Inf;
end
end

function [factor,inverseGain,relativeResidual] = ...
        construct_fixed_range_border_factor( ...
        matrix,probe,probeNorm,opts,basis,description)
matrixSize = size(matrix,1);
rangeBorder = sparse(basis.left);
rightNull = sparse(basis.right);
nullity = basis.nullity;
if size(rangeBorder,2) ~= nullity || size(rightNull,2) ~= nullity
    error('vi_cylinder_wnl_forced_schur_solve:KnownNullspaceSize', ...
        'The certified nullspace size is inconsistent for %s.',description);
end
nullTolerance = max( ...
    100*opts.forcedCylinderSchurNullTolerance,1.0e-9);
rightResidual = norm(matrix*rightNull,'fro') / ...
    max(norm(matrix,1)*norm(rightNull,'fro'),eps);
if rightResidual > nullTolerance
    error('vi_cylinder_wnl_forced_schur_solve:KnownNullspaceResidual', ...
        ['The certified right pressure-nullspace is not invariant for ', ...
         '%s: residual %.3e.'],description,rightResidual);
end
% The left nullspace varies with temporal frequency because U'*B need not
% vanish. Its reference value is nevertheless a valid fixed range-
% completion border whenever the augmented factor is nonsingular. Any
% nonzero border multiplier appears in the original physical residual,
% which remains the final acceptance gate.
bordered = [matrix,rangeBorder;rightNull',sparse(nullity,nullity)];
factor = struct('method','fixed-range-bordered-lu', ...
    'decomposition',decomposition(bordered,'lu'), ...
    'numberOfStateRows',matrixSize,'nullity',nullity);
rightHandSide = [probe;complex(zeros(nullity,size(probe,2)))];
warningStates = suppress_factor_warnings();
cleanup = onCleanup(@() restore_warning_states(warningStates));
trialAugmented = factor.decomposition\rightHandSide;
clear cleanup;
trialSolution = trialAugmented(1:matrixSize,:);
solutionNorm = sqrt(sum(abs(trialSolution).^2,1));
inverseGain = max(solutionNorm./probeNorm);
trialResidual = bordered*trialAugmented-rightHandSide;
relativeResidual = max( ...
    sqrt(sum(abs(trialResidual).^2,1))./probeNorm);
if any(~isfinite(real(trialAugmented(:)))) || ...
        any(~isfinite(imag(trialAugmented(:))))
    inverseGain = Inf;
    relativeResidual = Inf;
end
end

function [factor,inverseGain,relativeResidual] = ...
        construct_nullspace_bordered_factor( ...
        matrix,probe,probeNorm,opts,nullityHint,description)
matrixSize = size(matrix,1);
if isempty(nullityHint)
    numberRequested = min(opts.forcedCylinderSchurMaximumNullity, ...
        max(matrixSize-1,1));
else
    numberRequested = min(nullityHint,max(matrixSize-1,1));
end
[leftVectors,singularValues,rightVectors] = ...
    svds(matrix,numberRequested,'smallest');
[singularValues,order] = sort(real(diag(singularValues)),'ascend');
leftVectors = leftVectors(:,order);
rightVectors = rightVectors(:,order);
threshold = opts.forcedCylinderSchurNullTolerance*max(norm(matrix,1),1);
nullity = sum(singularValues <= threshold);
if nullity == 0
    error('vi_cylinder_wnl_forced_schur_solve:NoTemporalNullspace', ...
        ['The exact LU failed for %s, but no singular value fell below ', ...
         'the nullspace threshold %.3e.'],description,threshold);
end
if isempty(nullityHint) && nullity == numberRequested && ...
        numberRequested < matrixSize-1
    error('vi_cylinder_wnl_forced_schur_solve:TemporalNullspaceLimit', ...
        ['At least %d temporal null vectors were detected; increase ', ...
         'forcedCylinderSchurMaximumNullity.'],numberRequested);
end
leftNull = sparse(leftVectors(:,1:nullity));
rightNull = sparse(rightVectors(:,1:nullity));
bordered = [matrix,leftNull;rightNull',sparse(nullity,nullity)];
factor = struct('method','nullspace-bordered-lu','decomposition', ...
    decomposition(bordered,'lu'),'numberOfStateRows',matrixSize, ...
    'nullity',nullity);
rightHandSide = [probe;complex(zeros(nullity,size(probe,2)))];
warningStates = suppress_factor_warnings();
cleanup = onCleanup(@() restore_warning_states(warningStates));
trialAugmented = factor.decomposition\rightHandSide;
clear cleanup;
trialSolution = trialAugmented(1:matrixSize,:);
solutionNorm = sqrt(sum(abs(trialSolution).^2,1));
inverseGain = max(solutionNorm./probeNorm);
trialResidual = bordered*trialAugmented-rightHandSide;
relativeResidual = max( ...
    sqrt(sum(abs(trialResidual).^2,1))./probeNorm);
if any(~isfinite(real(trialAugmented(:)))) || ...
        any(~isfinite(imag(trialAugmented(:))))
    inverseGain = Inf;
    relativeResidual = Inf;
end
end

function solution = refined_factor_solve( ...
        factor,matrix,rightHandSide,numberOfSteps)
solution = factor_solve_safely(factor,rightHandSide);
if isstruct(factor) && ...
        strcmp(factor.method,'column-equilibrated-qr')
    % QR already returns the least-squares solution in the selected rank
    % subspace. Re-solving its orthogonal residual cannot improve that
    % objective and only repeats an expensive sparse backsolve.
    numberOfSteps = 0;
end
for refinementIndex = 1:numberOfSteps
    residual = rightHandSide-matrix*solution;
    oldNorm = norm(residual,'fro');
    if oldNorm <= 64*eps*max(norm(rightHandSide,'fro'),1)
        break;
    end
    correction = factor_solve_safely(factor,residual);
    trial = solution+correction;
    newNorm = norm(rightHandSide-matrix*trial,'fro');
    if ~isfinite(newNorm) || newNorm >= oldNorm*(1-64*eps)
        break;
    end
    solution = trial;
end
if any(~isfinite(real(solution(:)))) || ...
        any(~isfinite(imag(solution(:))))
    error('vi_cylinder_wnl_forced_schur_solve:NonfiniteFactorSolve', ...
        'A temporal-block factor solve returned a nonfinite value.');
end
end

function solution = factor_solve_safely(factor,rightHandSide)
warningStates = suppress_factor_warnings();
cleanup = onCleanup(@() restore_warning_states(warningStates));
if isstruct(factor)
    if strcmp(factor.method,'column-equilibrated-qr')
        solution = factor.columnScale.* ...
            (factor.decomposition\rightHandSide);
    elseif strcmp(factor.method,'equilibrated-lu')
        solution = factor.columnScale.* ...
            (factor.decomposition\(factor.rowScale.*rightHandSide));
    elseif factor.nullity > 0
        augmentedRightHandSide = [rightHandSide;complex(zeros( ...
            factor.nullity,size(rightHandSide,2)))];
        augmentedSolution = factor.decomposition\augmentedRightHandSide;
        solution = augmentedSolution(1:factor.numberOfStateRows,:);
    else
        solution = factor.decomposition\rightHandSide;
    end
else
    solution = factor\rightHandSide;
end
clear cleanup;
end

function warningStates = suppress_factor_warnings()
warningIds = {'MATLAB:singularMatrix', ...
    'MATLAB:nearlySingularMatrix', ...
    'MATLAB:illConditionedMatrix', ...
    'MATLAB:rankDeficientMatrix'};
warningStates = repmat(struct('identifier','','state',''), ...
    numel(warningIds),1);
for warningIndex = 1:numel(warningIds)
    warningStates(warningIndex) = warning('off', ...
        warningIds{warningIndex});
end
end

function restore_warning_states(warningStates)
for warningIndex = 1:numel(warningStates)
    warning(warningStates(warningIndex).state, ...
        warningStates(warningIndex).identifier);
end
end

function basis = pressure_nullspace_basis(matrix,mass,spec,opts)
basis = struct('left',complex(zeros(size(matrix,1),0)), ...
    'right',complex(zeros(size(matrix,1),0)), ...
    'nullity',0,'certified',false,'rightResidual',NaN, ...
    'massResidual',NaN,'leftResidual',NaN,'leftMassResidual',NaN);
if spec.m ~= 0 || ~opts.forcedCylinderSchurUseNullspaceBordering
    return;
end
numberRequested = min(opts.forcedCylinderSchurMaximumNullity, ...
    size(matrix,1)-1);
try
    [leftVectors,singularValues,rightVectors] = ...
        svds(matrix,numberRequested,'smallest');
catch
    return;
end
[singularValues,order] = sort(real(diag(singularValues)),'ascend');
leftVectors = leftVectors(:,order);
rightVectors = rightVectors(:,order);
threshold = opts.forcedCylinderSchurNullTolerance* ...
    max(norm(matrix,1),1);
nullity = sum(singularValues <= threshold);
if nullity == 0 || ...
        (nullity == numberRequested && numberRequested < size(matrix,1)-1)
    return;
end
leftNull = leftVectors(:,1:nullity);
rightNull = rightVectors(:,1:nullity);
massColumnNorm = full(sqrt(sum(abs(mass).^2,1))).';
algebraicColumns = massColumnNorm <= ...
    max(max(massColumnNorm)*1.0e-13,eps);
projectedRight = complex(zeros(size(rightNull)));
projectedRight(algebraicColumns,:) = rightNull(algebraicColumns,:);
[projectedRight,~] = qr(projectedRight,0);
rawResidual = norm(matrix*rightNull,'fro') / ...
    max(norm(matrix,1)*norm(rightNull,'fro'),eps);
projectedResidual = norm(matrix*projectedRight,'fro') / ...
    max(norm(matrix,1)*norm(projectedRight,'fro'),eps);
if projectedResidual <= max(100*rawResidual, ...
        10*opts.forcedCylinderSchurNullTolerance)
    rightNull = projectedRight;
end
basis.left = leftNull;
basis.right = rightNull;
basis.nullity = nullity;
basis.rightResidual = norm(matrix*rightNull,'fro') / ...
    max(norm(matrix,1)*norm(rightNull,'fro'),eps);
basis.massResidual = norm(mass*rightNull,'fro') / ...
    max(norm(mass,1)*norm(rightNull,'fro'),eps);
basis.leftResidual = norm(leftNull'*matrix,'fro') / ...
    max(norm(matrix,1)*norm(leftNull,'fro'),eps);
basis.leftMassResidual = norm(leftNull'*mass,'fro') / ...
    max(norm(mass,1)*norm(leftNull,'fro'),eps);
certificationTolerance = max( ...
    100*opts.forcedCylinderSchurNullTolerance,1.0e-10);
basis.certified = isfinite(basis.rightResidual) && ...
    isfinite(basis.massResidual) && isfinite(basis.leftResidual) && ...
    basis.rightResidual <= certificationTolerance && ...
    basis.massResidual <= certificationTolerance && ...
    basis.leftResidual <= certificationTolerance;
end

function leftNull = harmonic_left_compatibility_basis( ...
        factor,matrix,opts)
leftNull = complex(zeros(size(matrix,1),0));
if ~opts.forcedCylinderSchurUsePressureCompatibilityProjection || ...
        ~isstruct(factor) || factor.nullity <= 0 || ...
        ~strcmp(factor.method,'fixed-range-bordered-lu')
    return;
end
matrixSize = size(matrix,1);
nullity = factor.nullity;
rightHandSide = [complex(zeros(matrixSize,nullity));speye(nullity)];
warningStates = suppress_factor_warnings();
cleanup = onCleanup(@() restore_warning_states(warningStates));
augmentedSolution = factor.decomposition'\rightHandSide;
clear cleanup;
candidate = augmentedSolution(1:matrixSize,:);
if any(~isfinite(real(candidate(:)))) || ...
        any(~isfinite(imag(candidate(:))))
    return;
end
[candidate,~] = qr(candidate,0);
relativeResidual = norm(candidate'*matrix,'fro') / ...
    max(norm(matrix,1)*norm(candidate,'fro'),eps);
certificationTolerance = max( ...
    100*opts.forcedCylinderSchurNullTolerance,1.0e-9);
if isfinite(relativeResidual) && ...
        relativeResidual <= certificationTolerance
    leftNull = candidate;
end
end

function result = pressure_compatible_residual( ...
        residual,diagonals,basis,leftCompatibilityBasis,forcingNorm,opts)
result = struct('applied',false, ...
    'relativeResidual',norm(residual(:))/max(forcingNorm,eps), ...
    'compatibilityRelativeResidual',0,'analyzedHarmonics',0, ...
    'maximumNullResidual',NaN);
if ~opts.forcedCylinderSchurUsePressureCompatibilityProjection || ...
        ~basis.certified || basis.nullity == 0
    return;
end
if numel(leftCompatibilityBasis) ~= size(residual,2)
    return;
end
numberOfHarmonics = size(residual,2);
analysisFloor = opts.forcedSolveResidualTolerance*forcingNorm / ...
    max(100*sqrt(numberOfHarmonics),1);
projected = residual;
maximumNullResidual = 0;
analyzed = 0;
try
    for harmonicIndex = 1:numberOfHarmonics
        if norm(residual(:,harmonicIndex)) <= analysisFloor
            continue;
        end
        diagonal = diagonals{harmonicIndex};
        rightResidual = norm(diagonal*basis.right,'fro') / ...
            max(norm(diagonal,1)*norm(basis.right,'fro'),eps);
        if rightResidual > max( ...
                100*opts.forcedCylinderSchurNullTolerance,1.0e-9)
            return;
        end
        leftNull = leftCompatibilityBasis{harmonicIndex};
        if size(leftNull,2) ~= basis.nullity
            return;
        end
        relativeNullResidual = norm(leftNull'*diagonal,'fro') / ...
            max(norm(diagonal,1)*norm(leftNull,'fro'),eps);
        maximumNullResidual = max( ...
            maximumNullResidual,relativeNullResidual);
        if relativeNullResidual > max( ...
                100*opts.forcedCylinderSchurNullTolerance,1.0e-9)
            return;
        end
        projected(:,harmonicIndex) = residual(:,harmonicIndex)- ...
            leftNull*(leftNull'*residual(:,harmonicIndex));
        analyzed = analyzed+1;
    end
catch
    return;
end
removed = residual-projected;
result.applied = analyzed > 0;
result.relativeResidual = norm(projected(:))/max(forcingNorm,eps);
result.compatibilityRelativeResidual = ...
    norm(removed(:))/max(forcingNorm,eps);
result.analyzedHarmonics = analyzed;
result.maximumNullResidual = maximumNullResidual;
end

function range = finite_range(values)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    range = [NaN,NaN];
else
    range = [min(finiteValues),max(finiteValues)];
end
end

function value = maximum_finite(values)
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    value = NaN;
else
    value = max(finiteValues);
end
end

function [reduced,rightHandSide,previousPositions,nextPositions] = ...
        assemble_reduced_system( ...
        spec,activeColumns,particular,previousResponse,nextResponse)
numberOfHarmonics = numel(spec.n);
numberOfActiveColumns = numel(activeColumns);
reducedSize = numberOfActiveColumns*numberOfHarmonics;
reduced = speye(reducedSize);
rightHandSide = reshape(particular(activeColumns,:),reducedSize,1);
previousPositions = zeros(numberOfHarmonics,1);
nextPositions = zeros(numberOfHarmonics,1);
for harmonicIndex = 1:numberOfHarmonics
    rowIndices = reduced_indices( ...
        harmonicIndex,numberOfActiveColumns);
    previousPosition = find( ...
        spec.n == spec.n(harmonicIndex)-1,1);
    nextPosition = find(spec.n == spec.n(harmonicIndex)+1,1);
    if ~isempty(previousPosition)
        previousPositions(harmonicIndex) = previousPosition;
        columnIndices = reduced_indices( ...
            previousPosition,numberOfActiveColumns);
        reduced(rowIndices,columnIndices) = ...
            reduced(rowIndices,columnIndices) + ...
            previousResponse{harmonicIndex}(activeColumns,:);
    end
    if ~isempty(nextPosition)
        nextPositions(harmonicIndex) = nextPosition;
        columnIndices = reduced_indices( ...
            nextPosition,numberOfActiveColumns);
        reduced(rowIndices,columnIndices) = ...
            reduced(rowIndices,columnIndices) + ...
            nextResponse{harmonicIndex}(activeColumns,:);
    end
end
end

function residual = complete_temporal_residual( ...
        blocks,omega,spec,lambda,state,forcing)
numberOfHarmonics = numel(spec.n);
residual = complex(zeros(size(state)));
previousCoupling = -blocks.Lplus;
nextCoupling = -blocks.Lminus;
for harmonicIndex = 1:numberOfHarmonics
    diagonal = temporal_diagonal( ...
        blocks,omega,spec,lambda,harmonicIndex);
    value = diagonal*state(:,harmonicIndex)-forcing(:,harmonicIndex);
    previousPosition = find( ...
        spec.n == spec.n(harmonicIndex)-1,1);
    nextPosition = find(spec.n == spec.n(harmonicIndex)+1,1);
    if ~isempty(previousPosition)
        value = value+previousCoupling*state(:,previousPosition);
    end
    if ~isempty(nextPosition)
        value = value+nextCoupling*state(:,nextPosition);
    end
    residual(:,harmonicIndex) = value;
end
end

function indices = reduced_indices(position,blockSize)
indices = (position-1)*blockSize+(1:blockSize);
end
