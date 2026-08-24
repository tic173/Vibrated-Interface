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
    diagonalFactors = cell(numberOfHarmonics,1);
    previousResponse = cell(numberOfHarmonics,1);
    nextResponse = cell(numberOfHarmonics,1);
    maximumResponseResidual = 0;
    blockRegularization = zeros(numberOfHarmonics,1);
    blockInverseGain = nan(numberOfHarmonics,1);
    reducedRegularization = NaN;
    reducedInverseGain = NaN;

    factorClock = tic;
    for harmonicIndex = 1:numberOfHarmonics
        diagonal = temporal_diagonal( ...
            blocks,omega,spec,lambda,harmonicIndex);
        [diagonalFactors{harmonicIndex}, ...
            blockRegularization(harmonicIndex), ...
            blockInverseGain(harmonicIndex)] = ...
            factor_temporal_system(diagonal,opts,sprintf( ...
            'temporal block %d',harmonicIndex));
        if numberOfActiveColumns > 0
            responseRightHandSide = [ ...
                previousCoupling(:,activeColumns), ...
                nextCoupling(:,activeColumns)];
            response = refined_factor_solve( ...
                diagonalFactors{harmonicIndex},diagonal, ...
                responseRightHandSide, ...
                opts.forcedCylinderSchurBlockRefinementSteps);
            previousResponse{harmonicIndex} = ...
                response(:,1:numberOfActiveColumns);
            nextResponse{harmonicIndex} = ...
                response(:,numberOfActiveColumns+1:end);
            responseResidual = norm( ...
                diagonal*response-responseRightHandSide,'fro') / ...
                max(norm(responseRightHandSide,'fro'),eps);
            maximumResponseResidual = max( ...
                maximumResponseResidual,responseResidual);
        else
            previousResponse{harmonicIndex} = ...
                complex(zeros(ndof,0));
            nextResponse{harmonicIndex} = complex(zeros(ndof,0));
        end
    end
    information.factorSeconds = toc(factorClock);
    information.maximumResponseResidual = maximumResponseResidual;
    information.blockRegularizationRange = ...
        finite_range(blockRegularization);
    information.maximumBlockInverseGain = maximum_finite(blockInverseGain);

    particularClock = tic;
    forcingByHarmonic = reshape(forcing,ndof,numberOfHarmonics);
    particular = complex(zeros(ndof,numberOfHarmonics));
    for harmonicIndex = 1:numberOfHarmonics
        diagonal = temporal_diagonal( ...
            blocks,omega,spec,lambda,harmonicIndex);
        particular(:,harmonicIndex) = refined_factor_solve( ...
            diagonalFactors{harmonicIndex},diagonal, ...
            forcingByHarmonic(:,harmonicIndex), ...
            opts.forcedCylinderSchurBlockRefinementSteps);
    end
    information.particularSeconds = toc(particularClock);

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
            'reduced temporal Schur system');
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

    vector = stateByHarmonic(:);
    residual = complete_temporal_residual( ...
        blocks,omega,spec,lambda,stateByHarmonic,forcingByHarmonic);
    information.relativeResidual = norm(residual(:)) / ...
        max(norm(forcing),eps);
    information.available = all(isfinite(real(vector))) && ...
        all(isfinite(imag(vector))) && ...
        isfinite(information.relativeResidual);
    information.passedPhysicalGate = information.available && ...
        information.relativeResidual <= ...
        opts.forcedSolveResidualTolerance;
    information.totalSeconds = toc(solverClock);
    if information.passedPhysicalGate
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
    'reducedResidual',NaN,'maximumResponseResidual',NaN, ...
    'blockRegularizationRange',[NaN,NaN], ...
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
        factor_temporal_system(matrix,opts,description)
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
for shiftIndex = 1:numel(shifts)
    shift = shifts(shiftIndex);
    shiftedMatrix = matrix;
    if shift > 0
        shiftedMatrix = shiftedMatrix+shift*identity;
    end
    try
        trialFactor = decomposition(shiftedMatrix,'lu');
        trialSolution = factor_solve_safely(trialFactor,probe);
        solutionNorm = sqrt(sum(abs(trialSolution).^2,1));
        trialGain = max(solutionNorm./probeNorm);
        trialResidual = shiftedMatrix*trialSolution-probe;
        relativeTrialResidual = max( ...
            sqrt(sum(abs(trialResidual).^2,1))./probeNorm);
        finiteTrial = all(isfinite(real(trialSolution(:)))) && ...
            all(isfinite(imag(trialSolution(:))));
        if finiteTrial && ...
                trialGain <= ...
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

function solution = refined_factor_solve( ...
        factor,matrix,rightHandSide,numberOfSteps)
solution = factor_solve_safely(factor,rightHandSide);
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
solution = factor\rightHandSide;
clear cleanup;
end

function warningStates = suppress_factor_warnings()
warningIds = {'MATLAB:singularMatrix', ...
    'MATLAB:nearlySingularMatrix', ...
    'MATLAB:illConditionedMatrix'};
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
