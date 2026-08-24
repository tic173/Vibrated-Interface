function solution = wnl_solve_forced(model, spec, forcingCoeff, neutralModes, userOpts)
%WNL_SOLVE_FORCED Solve a second-order Floquet block, with bordering.
%
% For retained modes in the output block, the system is
%
%   [ A       Bslow*Phi ] [q]      [f]
%   [ Y'*Bslow    0     ] [lambda] [0].
%
% lambda reports resonant forcing. A nonzero Y'*f means that the associated
% quadratic amplitude term must be retained.

if nargin < 4 || isempty(neutralModes)
    neutralModes = {};
end
if nargin < 5
    userOpts = struct();
end
opts = wnl_options(userOpts);
forcedSolveWallClock = tic;
block = model.block(spec);
forcing = forcingCoeff(:);
n = size(block.A, 1);
if numel(forcing) ~= n
    error('wnl_solve_forced:ForcingSize', ...
        'Forcing has %d entries; block requires %d.', numel(forcing), n);
end

matches = wnl_matching_modes(neutralModes, spec);
projections = complex(zeros(numel(matches), 1));
for j = 1:numel(matches)
    projections(j) = matches{j}.left' * forcing;
end

rankDiagnostics = disabled_rank_information( ...
    'rank-aware solve was not required');
modelSolveDiagnostics = disabled_model_solve_information( ...
    'no model-specific forced solver was used');
if isempty(matches)
    hasModelSolve = isfield(block, 'solve') && ...
        isa(block.solve, 'function_handle');
    if hasModelSolve
        modelCandidate = block.solve(forcing, spec, opts);
        structuredModelSolve = isstruct(modelCandidate) && ...
            isfield(modelCandidate,'vector');
        if structuredModelSolve
            q = modelCandidate.vector(:);
            if isfield(modelCandidate,'diagnostics') && ...
                    isstruct(modelCandidate.diagnostics)
                modelSolveDiagnostics = modelCandidate.diagnostics;
            else
                modelSolveDiagnostics = disabled_model_solve_information( ...
                    'model solver returned no diagnostics');
                modelSolveDiagnostics.attempted = true;
            end
        else
            q = modelCandidate(:);
            modelSolveDiagnostics = disabled_model_solve_information( ...
                'legacy model solver returned a vector');
            modelSolveDiagnostics.attempted = true;
            modelSolveDiagnostics.available = true;
        end
        if numel(q) ~= n
            error('wnl_solve_forced:ModelSolveSize', ...
                ['The model-specific solver returned %d rather than ', ...
                 '%d entries.'],numel(q),n);
        end
        [q, gmresDiagnostics] = block_gmres_solution( ...
            block.A, forcing, q, spec, opts);
        % A structured model solve may declare itself unavailable (for
        % example, a singular temporal diagonal in the cylinder Schur
        % elimination). In that case retain the established global fallback.
        if structuredModelSolve && ...
                (~isfield(modelSolveDiagnostics,'available') || ...
                ~modelSolveDiagnostics.available || ...
                ~isfield(modelSolveDiagnostics,'passedPhysicalGate') || ...
                ~modelSolveDiagnostics.passedPhysicalGate)
            currentResidual = norm(forcing-block.A*q) / ...
                max(norm(forcing),eps);
            if currentResidual > opts.forcedSolveResidualTolerance && ...
                    opts.forcedUseRankAwareMinimumNorm
                gmresDiagnostics.minimumNormFallbackUsed = true;
                [minimumNormCandidate,rankDiagnostics] = ...
                    wnl_rank_aware_minimum_norm( ...
                    block.A,forcing,block.Bslow,opts);
                minimumNormResidual = norm( ...
                    forcing-block.A*minimumNormCandidate) / ...
                    max(norm(forcing),eps);
                gmresDiagnostics.minimumNormFallbackResidual = ...
                    minimumNormResidual;
                if isfinite(minimumNormResidual) && ...
                        minimumNormResidual < currentResidual
                    q = minimumNormCandidate;
                end
            end
        end
    else
        % The temporal-block GMRES path is optional. Models without a
        % dedicated forced solver use the rank-aware, column-equilibrated
        % minimum-norm solve below; unlike row equilibration, column
        % equilibration does not change the physical least-squares objective.
        q = complex(zeros(n,1));
        [q, gmresDiagnostics] = block_gmres_solution( ...
            block.A, forcing, q, spec, opts);
        currentResidual = norm(forcing-block.A*q) / ...
            max(norm(forcing),eps);
        if currentResidual > opts.forcedSolveResidualTolerance && ...
                opts.forcedUseRankAwareMinimumNorm
            gmresDiagnostics.minimumNormFallbackUsed = true;
            [minimumNormCandidate,rankDiagnostics] = ...
                wnl_rank_aware_minimum_norm( ...
                block.A,forcing,block.Bslow,opts);
            minimumNormResidual = norm( ...
                forcing-block.A*minimumNormCandidate) / ...
                max(norm(forcing),eps);
            gmresDiagnostics.minimumNormFallbackResidual = ...
                minimumNormResidual;
            if isfinite(minimumNormResidual) && ...
                    minimumNormResidual < currentResidual
                q = minimumNormCandidate;
            end
        elseif currentResidual > opts.forcedSolveResidualTolerance
            gmresDiagnostics.minimumNormFallbackUsed = true;
            minimumNormCandidate = minimum_norm_solve( ...
                block.A,forcing,opts);
            minimumNormResidual = norm( ...
                forcing-block.A*minimumNormCandidate) / ...
                max(norm(forcing),eps);
            gmresDiagnostics.minimumNormFallbackResidual = ...
                minimumNormResidual;
            rankDiagnostics = disabled_rank_information( ...
                ['rank-aware solve disabled; legacy minimum-norm ', ...
                 'candidate used']);
            rankDiagnostics.relativeResidual = minimumNormResidual;
            if isfinite(minimumNormResidual) && ...
                    minimumNormResidual < currentResidual
                q = minimumNormCandidate;
            end
        else
            rankDiagnostics = disabled_rank_information( ...
                'the optional GMRES seed already passed the physical gate');
        end
    end
    [q, solveDiagnostics, algebraicCompletion] = ...
        refine_forced_solution(block.A, block.Bslow, forcing, q, opts);
    lambda = complex(zeros(0, 1));
else
    phi = complex(zeros(n, numel(matches)));
    y = complex(zeros(n, numel(matches)));
    for j = 1:numel(matches)
        phi(:, j) = matches{j}.vector;
        y(:, j) = matches{j}.left;
    end
    topRight = block.Bslow * phi;
    bottomLeft = y' * block.Bslow;
    bordered = [block.A, topRight; ...
        bottomLeft, sparse(numel(matches), numel(matches))];
    rhs = [forcing; complex(zeros(numel(matches), 1))];
    borderedSolution = minimum_norm_solve(bordered, rhs, opts);
    gmresDiagnostics = struct('attempted', false, ...
        'accepted', false, 'flag', NaN, 'iterations', 0, ...
        'relativeResidual', NaN, ...
        'improvedSeed',false,'passedPhysicalGate',false, ...
        'minimumNormFallbackUsed',true, ...
        'minimumNormFallbackResidual',NaN);
    rankDiagnostics = disabled_rank_information( ...
        'bordered resonant solve uses its projection constraints');
    [borderedSolution, solveDiagnostics] = refine_bordered_solution( ...
        bordered, rhs, borderedSolution, opts);
    q = borderedSolution(1:n);
    lambda = borderedSolution(n + 1:end);
    algebraicCompletion = disabled_completion_information( ...
        'bordered resonant solve: algebraic completion not applied');
end

solution = struct();
solution.spec = spec;
solution.field = wnl_make_field(spec, q);
solution.vector = q;
solution.forcing = forcing;
solution.projections = projections;
solution.lambda = lambda;
solution.quadraticResonance = any(abs(projections) > ...
    opts.resonanceTolerance * max(1.0, norm(forcing)));
solution.unborderedResidual = norm(block.A * q - forcing);
solution.equationResidual = solution.unborderedResidual;
solution.constraintResidual = 0.0;
if ~isempty(matches)
    solution.equationResidual = norm(block.A * q + ...
        block.Bslow * phi * lambda - forcing);
    solution.constraintResidual = norm(y' * (block.Bslow * q));
end
solution.forcingNorm = norm(forcing);
solution.relativeEquationResidual = solution.equationResidual / ...
    max(solution.forcingNorm, eps);
solution.relativeConstraintResidual = solution.constraintResidual / ...
    max(norm(block.Bslow*q), 1.0);
solution.valid = isfinite(solution.relativeEquationResidual) && ...
    isfinite(solution.relativeConstraintResidual) && ...
    solution.relativeEquationResidual <= ...
        opts.forcedSolveResidualTolerance && ...
    solution.relativeConstraintResidual <= ...
        opts.forcedSolveConstraintTolerance;
exploratoryTolerance = forced_exploratory_tolerance(opts);
solution.exploratoryResidualTolerance = exploratoryTolerance;
solution.exploratoryUsable = ...
    isfinite(solution.relativeEquationResidual) && ...
    isfinite(solution.relativeConstraintResidual) && ...
    solution.relativeEquationResidual <= exploratoryTolerance && ...
    solution.relativeConstraintResidual <= ...
        opts.forcedSolveConstraintTolerance;
solution.solveDiagnostics = solveDiagnostics;
solution.gmresDiagnostics = gmresDiagnostics;
solution.rankAwareDiagnostics = rankDiagnostics;
solution.modelSolveDiagnostics = modelSolveDiagnostics;
solution.algebraicCompletion = algebraicCompletion;
solution.fullStateNorm = norm(q);
solution.descriptorStateNorm = norm(block.Bslow*q);
solution.fullToDescriptorNormRatio = solution.fullStateNorm / ...
    max(solution.descriptorStateNorm,eps);
solution.residualReport = [];
solution.solveSeconds = toc(forcedSolveWallClock);

if opts.verbose
    fprintf(['WNL forced solve %s: equation residual %.3e ', ...
        '(forcing-relative %.3e), constraint residual %.3e, ', ...
        'valid = %d, refinement flag/iterations/restarts = %d/%d/%d\n'], ...
        spec.label, solution.equationResidual, ...
        solution.relativeEquationResidual, solution.constraintResidual, ...
        solution.valid, solveDiagnostics.flag, ...
        solveDiagnostics.totalIterations, solveDiagnostics.restarts);
    fprintf('  forced solve wall-clock time = %.3f s\n', ...
        solution.solveSeconds);
    if ~solution.valid && solution.exploratoryUsable
        fprintf(['  exploratory forced-field ceiling passed: %.3e ', ...
            '<= %.3e (strict validity remains false)\n'], ...
            solution.relativeEquationResidual,exploratoryTolerance);
    end
    if modelSolveDiagnostics.attempted
        fprintf(['  model forced solver: %s; available/gate = %d/%d; ', ...
            'physical residual %.3e\n'], ...
            model_solve_text(modelSolveDiagnostics,'method','unspecified'), ...
            model_solve_logical(modelSolveDiagnostics,'available'), ...
            model_solve_logical(modelSolveDiagnostics, ...
                'passedPhysicalGate'), ...
            model_solve_numeric(modelSolveDiagnostics, ...
                'relativeResidual'));
        if isfield(modelSolveDiagnostics,'reducedDimension')
            fprintf(['  temporal Schur dimension/active columns = %d/%d; ', ...
                'block-response/reduced residual %.3e/%.3e; ', ...
                'factor/total time %.3f/%.3f s\n'], ...
                modelSolveDiagnostics.reducedDimension, ...
                modelSolveDiagnostics.numberOfActiveColumns, ...
                modelSolveDiagnostics.maximumResponseResidual, ...
                modelSolveDiagnostics.reducedResidual, ...
                modelSolveDiagnostics.factorSeconds, ...
                modelSolveDiagnostics.totalSeconds);
        end
        if isfield(modelSolveDiagnostics,'blockRegularizationRange')
            fprintf(['  temporal Schur regularization range ', ...
                '[%.3e, %.3e]; max inverse gain %.3e; ', ...
                'reduced shift %.3e\n'], ...
                modelSolveDiagnostics.blockRegularizationRange(1), ...
                modelSolveDiagnostics.blockRegularizationRange(2), ...
                modelSolveDiagnostics.maximumBlockInverseGain, ...
                modelSolveDiagnostics.reducedRegularization);
        end
        if isfield(modelSolveDiagnostics,'stopReason')
            fprintf('  model forced-solver stop: %s\n', ...
                modelSolveDiagnostics.stopReason);
        end
    end
    if gmresDiagnostics.attempted
        fprintf(['  scaled block-GMRES flag/iterations/physical ', ...
            'residual/improved-seed/gate ', ...
            '= %d/%d/%.3e/%d/%d\n'], gmresDiagnostics.flag, ...
            gmresDiagnostics.iterations, ...
            gmresDiagnostics.relativeResidual, ...
            gmresDiagnostics.improvedSeed, ...
            gmresDiagnostics.passedPhysicalGate);
        if isfield(gmresDiagnostics,'scaledRelativeResidual')
            fprintf(['  scaled block-GMRES residual %.3e; initial ', ...
                'physical residual %.3e; factor time %.3f s; ', ...
                'Krylov time %.3f s\n'], ...
                gmresDiagnostics.scaledRelativeResidual, ...
                gmresDiagnostics.initialRelativeResidual, ...
                gmresDiagnostics.factorSeconds, ...
                gmresDiagnostics.gmresSeconds);
        end
        if isfield(gmresDiagnostics,'blockRegularization') && ...
                ~isempty(gmresDiagnostics.blockRegularization)
            fprintf(['  forced block regularization range ', ...
                '[%.3e, %.3e], maximum probe inverse gain %.3e\n'], ...
                min(gmresDiagnostics.blockRegularization), ...
                max(gmresDiagnostics.blockRegularization), ...
                max(gmresDiagnostics.blockInverseGain));
        end
    end
    if rankDiagnostics.attempted
        if isnan(rankDiagnostics.selectedTolerance)
            toleranceText = 'default';
        else
            toleranceText = sprintf('%.3e', ...
                rankDiagnostics.selectedTolerance);
        end
        fprintf(['  rank-aware minimum norm: attempts %d, selected ', ...
            '%d, rank tolerance %s, physical residual %.3e, gate=%d\n'], ...
            rankDiagnostics.numberOfAttempts, ...
            rankDiagnostics.selectedAttempt,toleranceText, ...
            rankDiagnostics.relativeResidual, ...
            rankDiagnostics.passedPhysicalGate);
        fprintf(['  rank-aware state/descriptor ratio %.3e; ', ...
            'column-scale range [%.3e, %.3e]\n'], ...
            rankDiagnostics.fullToDescriptorNormRatio, ...
            rankDiagnostics.columnScaleRange(1), ...
            rankDiagnostics.columnScaleRange(2));
    end
    if algebraicCompletion.attempted
        fprintf(['  forced DAE completion available/accepted = %d/%d, ', ...
            'residual %.3e -> %.3e, state norm %.3e -> %.3e\n'], ...
            algebraicCompletion.available, ...
            algebraicCompletion.accepted, ...
            algebraicCompletion.initialResidualNorm, ...
            algebraicCompletion.finalResidualNorm, ...
            algebraicCompletion.initialFullNorm, ...
            algebraicCompletion.finalFullNorm);
        if isfield(algebraicCompletion,'acceptance')
            fprintf(['  forced DAE seed selection: residual-improved=%d, ', ...
                'substantial-compaction=%d, norm ratio %.3e, %s\n'], ...
                algebraicCompletion.acceptance.residualImproved, ...
                algebraicCompletion.acceptance.substantialCompaction, ...
                algebraicCompletion.acceptance.fullNormRatio, ...
                algebraicCompletion.acceptance.reason);
        end
    end
    fprintf(['  forced state norm/descriptor norm/ratio = ', ...
        '%.3e / %.3e / %.3e\n'],solution.fullStateNorm, ...
        solution.descriptorStateNorm,solution.fullToDescriptorNormRatio);
    if isfield(solveDiagnostics,'pressureSafeRefinement') && ...
            solveDiagnostics.pressureSafeRefinement && ...
            ~isempty(solveDiagnostics.correctionEquationResiduals)
        fprintf(['  pressure-safe correction equation residuals: ', ...
            'min %.3e, final %.3e; trial completion attempts = %d\n'], ...
            min(solveDiagnostics.correctionEquationResiduals), ...
            solveDiagnostics.correctionEquationResiduals(end), ...
            solveDiagnostics.trialCompletionCount);
    end
    if ~solution.valid && isfield(model,'residualLayout')
        solution.residualReport = wnl_forced_residual_report( ...
            solution,block,model.residualLayout);
    end
    if solution.quadraticResonance
        fprintf('  Quadratic resonant projections:');
        fprintf(' %.3e', abs(projections));
        fprintf('\n');
    end
end

function [x, information] = block_gmres_solution(A, b, x, spec, opts)
information = struct('attempted', false, 'accepted', false, ...
    'improvedSeed',false,'passedPhysicalGate',false, ...
    'minimumNormFallbackUsed',false, ...
    'minimumNormFallbackResidual',NaN, ...
    'flag', NaN, 'iterations', 0, 'relativeResidual', NaN, ...
    'scaledRelativeResidual',NaN,'initialRelativeResidual',NaN, ...
    'factorSeconds',0,'gmresSeconds',0, ...
    'blockRegularization',[],'blockInverseGain',[], ...
    'cycles',0,'stagnated',false);
numberOfHarmonics = numel(spec.n);
if ~opts.forcedUseBlockGmres || numberOfHarmonics <= 1 || ...
        size(A,1) ~= spec.ndof*numberOfHarmonics || ...
        size(A,2) ~= size(A,1)
    return;
end
try
    validate_forced_block_gmres_options(opts);
    initialResidual = norm(b-A*x)/max(norm(b),eps);
    information.initialRelativeResidual = initialResidual;
    if initialResidual <= opts.forcedSolveResidualTolerance
        information.relativeResidual = initialResidual;
        information.passedPhysicalGate = true;
        information.skipReason = ...
            'initial field already passed the physical residual gate';
        return;
    end
    information.attempted = true;

    [scaledA,rowScaledB,rowScale] = ...
        equilibrate_forced_equations(A,b);
    [scaledA,columnScale] = ...
        equilibrate_forced_columns(scaledA);
    scaledInitial = x./columnScale;
    if any(~isfinite(real(scaledInitial))) || ...
            any(~isfinite(imag(scaledInitial)))
        scaledInitial = complex(zeros(size(x)));
    end

    factorClock = tic;
    [factors,blockRegularization,blockInverseGain] = ...
        factor_forced_temporal_blocks(scaledA,spec.ndof, ...
            numberOfHarmonics,opts);
    information.factorSeconds = toc(factorClock);
    information.blockRegularization = blockRegularization;
    information.blockInverseGain = blockInverseGain;
    information.rowScaleRange = positive_forced_range(rowScale);
    information.columnScaleRange = positive_forced_range(columnScale);
    preconditioner = @(value) apply_block_inverse( ...
        value,factors,spec.ndof,numberOfHarmonics);

    validateattributes(opts.forcedGmresStagnationCycles,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(opts.forcedGmresMinimumCycleImprovement, ...
        {'numeric'},{'scalar','real','nonnegative','<',1,'finite'});

    gmresClock = tic;
    bestScaled = scaledInitial;
    bestCandidate = x;
    bestPhysicalResidual = initialResidual;
    bestScaledResidual = norm(rowScaledB-scaledA*bestScaled) / ...
        max(norm(rowScaledB),eps);
    totalIterations = 0;
    stagnantCycles = 0;
    lastFlag = 1;
    residualHistory = cell(opts.forcedGmresMaxCycles,1);
    for cycle = 1:opts.forcedGmresMaxCycles
        previousBest = bestPhysicalResidual;
        [trialScaled,lastFlag,reportedResidual,iteration,history] = ...
            gmres(scaledA,rowScaledB,opts.forcedGmresRestart, ...
            opts.forcedSolveRefinementTolerance,1,preconditioner,[], ...
            bestScaled);
        totalIterations = totalIterations + ...
            gmres_iteration_count(iteration,opts.forcedGmresRestart);
        residualHistory{cycle} = history;
        trial = columnScale.*trialScaled;
        trialFinite = all(isfinite(real(trial))) && ...
            all(isfinite(imag(trial)));
        if trialFinite
            trialPhysicalResidual = norm(b-A*trial)/max(norm(b),eps);
        else
            trialPhysicalResidual = Inf;
        end
        if trialPhysicalResidual < bestPhysicalResidual*(1-64*eps)
            bestCandidate = trial;
            bestScaled = trialScaled;
            bestPhysicalResidual = trialPhysicalResidual;
            bestScaledResidual = reportedResidual;
        end
        relativeImprovement = (previousBest-bestPhysicalResidual) / ...
            max(previousBest,eps);
        if relativeImprovement < ...
                opts.forcedGmresMinimumCycleImprovement
            stagnantCycles = stagnantCycles+1;
        else
            stagnantCycles = 0;
        end
        information.cycles = cycle;
        if bestPhysicalResidual <= opts.forcedSolveResidualTolerance
            break;
        end
        if stagnantCycles >= opts.forcedGmresStagnationCycles
            information.stagnated = true;
            break;
        end
    end
    information.gmresSeconds = toc(gmresClock);
    information.flag = lastFlag;
    information.iterations = totalIterations;
    information.relativeResidual = bestPhysicalResidual;
    information.scaledRelativeResidual = bestScaledResidual;
    information.reportedRelativeResidual = bestScaledResidual;
    information.residualHistory = residualHistory(1:information.cycles);
    information.improvedSeed = isfinite(bestPhysicalResidual) && ...
        bestPhysicalResidual < initialResidual*(1-64*eps);
    information.passedPhysicalGate = isfinite(bestPhysicalResidual) && ...
        bestPhysicalResidual <= opts.forcedSolveResidualTolerance;
    % Retain the former field for callers written against V39.  It means
    % only that the candidate improved the seed, not that it passed.
    information.accepted = information.improvedSeed;
    if information.improvedSeed
        x = bestCandidate;
    end
catch solverError
    information.errorIdentifier = solverError.identifier;
    information.errorMessage = solverError.message;
end

function count = gmres_iteration_count(iteration,restart)
if numel(iteration) > 1
    count = max(iteration(1)-1,0)*restart+iteration(2);
else
    count = iteration;
end
end
end

function value = apply_block_inverse(rhs,factors,ndof,numberOfHarmonics)
value = complex(zeros(size(rhs)));
for harmonicIndex = 1:numberOfHarmonics
    indices = (harmonicIndex-1)*ndof+(1:ndof);
    value(indices,:) = factors{harmonicIndex}\rhs(indices,:);
end
if any(~isfinite(real(value(:)))) || ...
        any(~isfinite(imag(value(:))))
    error('wnl_solve_forced:BlockPreconditionerNonfinite', ...
        ['The regularized temporal-block inverse returned a ', ...
         'nonfinite value.']);
end
end

function [factors,regularization,inverseGain] = ...
        factor_forced_temporal_blocks(A,ndof,numberOfHarmonics,opts)
factors = cell(numberOfHarmonics,1);
regularization = zeros(numberOfHarmonics,1);
inverseGain = inf(numberOfHarmonics,1);
identity = speye(ndof);
probeIndex = (1:ndof).';
probe = [complex(ones(ndof,1)),exp(1i*sqrt(2)*probeIndex)];
probeNorm = sqrt(sum(abs(probe).^2,1));
for harmonicIndex = 1:numberOfHarmonics
    indices = (harmonicIndex-1)*ndof+(1:ndof);
    diagonalBlock = A(indices,indices);
    oneNorm = norm(diagonalBlock,1);
    infinityNorm = norm(diagonalBlock,inf);
    blockScale = max(sqrt(max(oneNorm,0)*max(infinityNorm,0)),1.0);
    lastMessage = '';
    for attempt = 1:opts.forcedBlockRegularizationAttempts
        shift = opts.forcedBlockRegularization*blockScale * ...
            opts.forcedBlockRegularizationGrowth^(attempt-1);
        shiftedBlock = diagonalBlock+shift*identity;
        try
            trialFactor = decomposition(shiftedBlock,'lu');
            trialSolution = trialFactor\probe;
            solutionNorm = sqrt(sum(abs(trialSolution).^2,1));
            trialGain = max(solutionNorm./probeNorm);
            trialResidual = shiftedBlock*trialSolution-probe;
            relativeTrialResidual = max( ...
                sqrt(sum(abs(trialResidual).^2,1))./probeNorm);
            finiteTrial = all(isfinite(real(trialSolution(:)))) && ...
                all(isfinite(imag(trialSolution(:))));
            if finiteTrial && ...
                    trialGain <= opts.forcedBlockMaximumInverseGain && ...
                    relativeTrialResidual <= 1.0e-5
                factors{harmonicIndex} = trialFactor;
                regularization(harmonicIndex) = shift;
                inverseGain(harmonicIndex) = trialGain;
                break;
            end
            lastMessage = sprintf( ...
                'inverse gain %.3e, factor residual %.3e', ...
                trialGain,relativeTrialResidual);
        catch factorError
            lastMessage = factorError.message;
        end
    end
    if isempty(factors{harmonicIndex})
        error('wnl_solve_forced:BlockPreconditionerFactorization', ...
            ['Could not construct a finite regularized inverse for ', ...
             'forced temporal block %d after %d attempts (%s).'], ...
            harmonicIndex,opts.forcedBlockRegularizationAttempts, ...
            lastMessage);
    end
end
end

function validate_forced_block_gmres_options(opts)
validateattributes(opts.forcedGmresRestart,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opts.forcedGmresMaxCycles,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opts.forcedBlockRegularization,{'numeric'}, ...
    {'scalar','real','positive','finite'});
validateattributes(opts.forcedBlockRegularizationGrowth,{'numeric'}, ...
    {'scalar','real','>',1,'finite'});
validateattributes(opts.forcedBlockRegularizationAttempts,{'numeric'}, ...
    {'scalar','integer','positive'});
validateattributes(opts.forcedBlockMaximumInverseGain,{'numeric'}, ...
    {'scalar','real','positive','finite'});
end

function range = positive_forced_range(values)
positive = values(isfinite(values) & values > 0);
if isempty(positive)
    range = [NaN,NaN];
else
    range = [min(positive),max(positive)];
end
end

function [x, diagnostics, completion] = refine_forced_solution( ...
        A, Bslow, b, x, opts)
% Pressure-safe iterative refinement for a nonresonant forced field.
%
% The initial minimum-norm solution can have an enormous zero-mass
% pressure/gauge component.  First compact that algebraic component while
% keeping all descriptor-supported variables fixed.  Subsequent correction
% solves use physical-coordinate Tikhonov regularization and retain the
% candidate with the smallest UNEQUILIBRATED equation residual.
scale = max(norm(b), eps);
[acceptanceTolerance,refinementTolerance] = ...
    forced_refinement_tolerances(opts);
initialRelativeResidual = norm(b-A*x)/scale;
if initialRelativeResidual <= acceptanceTolerance
    exploratoryEarlyStop = initialRelativeResidual > ...
        opts.forcedSolveResidualTolerance;
    diagnostics = struct('flag',0,'totalIterations',0,'restarts',0, ...
        'relativeResidual',initialRelativeResidual, ...
        'rawRefinement',struct('flag',0,'totalIterations',0, ...
            'restarts',0,'relativeResidual',initialRelativeResidual, ...
            'pressureSafeRefinement',false), ...
        'correctionEquationResiduals',zeros(0,1), ...
        'acceptedStepLengths',zeros(0,1), ...
        'trialCompletionCount',0,'pressureSafeRestarts',0, ...
        'pressureSafeRefinement',false, ...
        'exploratoryEarlyStop',exploratoryEarlyStop);
    if exploratoryEarlyStop
        stopReason = ['initial field passed the exploratory residual ', ...
            'ceiling; strict forced-field validity remains false'];
    else
        stopReason = 'initial field passed the forced physical residual gate';
    end
    completion = disabled_completion_information(stopReason);
    return;
end
rawOpts = opts;
rawOpts.forcedSolveRefinementTolerance = refinementTolerance;
[x,rawDiagnostics] = refine_bordered_solution(A,b,x,rawOpts);
bestX = x;
bestRelativeResidual = norm(b-A*x)/scale;
currentX = x;
currentRelativeResidual = bestRelativeResidual;
completion = disabled_completion_information( ...
    'raw residual correction passed the forced-field gate');
if bestRelativeResidual <= acceptanceTolerance
    diagnostics = rawDiagnostics;
    diagnostics.rawRefinement = rawDiagnostics;
    diagnostics.correctionEquationResiduals = zeros(0,1);
    diagnostics.acceptedStepLengths = zeros(0,1);
    diagnostics.trialCompletionCount = 0;
    diagnostics.pressureSafeRefinement = false;
    diagnostics.exploratoryEarlyStop = bestRelativeResidual > ...
        opts.forcedSolveResidualTolerance;
    return;
end
completion.stopReason = 'forced algebraic completion disabled';
if opts.forcedUseAlgebraicCompletion
    [completed, completion] = complete_forced_state( ...
        A,Bslow,currentX,b,opts);
    [completion.accepted,completion.acceptance] = ...
        accept_forced_completion_seed(completion,scale,opts);
    if completion.accepted
        currentX = completed;
        currentRelativeResidual = norm(b-A*currentX)/scale;
        if currentRelativeResidual < bestRelativeResidual
            bestX = currentX;
            bestRelativeResidual = currentRelativeResidual;
        end
    end
end

totalIterations = rawDiagnostics.totalIterations;
lastFlag = rawDiagnostics.flag;
restarts = 0;
correctionEquationResiduals = nan(opts.forcedSolveMaxRestarts,1);
acceptedStepLengths = nan(opts.forcedSolveMaxRestarts,1);
trialCompletionCount = 0;
for attempt = 1:opts.forcedSolveMaxRestarts
    if bestRelativeResidual <= refinementTolerance
        break;
    end
    residual = b-A*currentX;
    [correction, flag, iterationCount, correctionInformation] = ...
        regularized_forced_correction(A,residual,opts);
    totalIterations = totalIterations + iterationCount;
    correctionEquationResiduals(attempt) = ...
        correctionInformation.equationResidual;
    lastFlag = flag;
    restarts = attempt;
    if ~correctionInformation.available
        break;
    end

    [candidate, candidateRelativeResidual, stepLength] = ...
        best_forced_line_search(A,b,currentX,correction, ...
        currentRelativeResidual,scale,opts.forcedLineSearchMaxCuts);

    % The full correction is the most physically useful state to complete:
    % its dynamic variables solve the equilibrated equation, while the
    % completion recomputes only algebraic variables.  Compare it with the
    % best raw line-search trial using the original equation residual.
    if opts.forcedUseAlgebraicCompletion
        fullTrial = currentX+correction;
        [completedTrial, trialCompletion] = complete_forced_state( ...
            A,Bslow,fullTrial,b,opts);
        trialCompletionCount = trialCompletionCount+1;
        completedRelativeResidual = norm(b-A*completedTrial)/scale;
        if trialCompletion.available && trialCompletion.valid && ...
                isfinite(completedRelativeResidual) && ...
                completedRelativeResidual < candidateRelativeResidual
            candidate = completedTrial;
            candidateRelativeResidual = completedRelativeResidual;
            stepLength = 1.0;
        end
    end

    improvedCurrent = candidateRelativeResidual < ...
        currentRelativeResidual*(1-64*eps);
    if ~improvedCurrent
        break;
    end
    currentX = candidate;
    currentRelativeResidual = candidateRelativeResidual;
    acceptedStepLengths(attempt) = stepLength;
    if currentRelativeResidual < bestRelativeResidual
        bestX = currentX;
        bestRelativeResidual = currentRelativeResidual;
    end
end
x = bestX;
diagnostics = struct('flag', lastFlag, ...
    'totalIterations', totalIterations, ...
    'restarts', rawDiagnostics.restarts+restarts, ...
    'relativeResidual', bestRelativeResidual, ...
    'correctionEquationResiduals', ...
        correctionEquationResiduals(1:restarts), ...
    'acceptedStepLengths',acceptedStepLengths(1:restarts), ...
    'trialCompletionCount',trialCompletionCount, ...
    'rawRefinement',rawDiagnostics, ...
    'pressureSafeRestarts',restarts, ...
    'pressureSafeRefinement',true, ...
    'exploratoryEarlyStop',bestRelativeResidual > ...
        opts.forcedSolveResidualTolerance && ...
        bestRelativeResidual <= acceptanceTolerance);
end

function [completed, information] = complete_forced_state( ...
        A,Bslow,state,rightHandSide,opts)
completionOpts = opts;
completionOpts.eigenpairAlgebraicCompletionRegularization = ...
    opts.forcedAlgebraicCompletionRegularization;
completionOpts.eigenpairAlgebraicCompletionSolveTolerance = ...
    opts.forcedAlgebraicCompletionSolveTolerance;
completionOpts.eigenpairAlgebraicCompletionEquationTolerance = ...
    opts.forcedAlgebraicCompletionEquationTolerance;
completionOpts.eigenpairAlgebraicCompletionMaxIterations = ...
    opts.forcedAlgebraicCompletionMaxIterations;
completionOpts.eigenpairAlgebraicCompletionMaxRestarts = ...
    opts.forcedAlgebraicCompletionMaxRestarts;
[completed,information] = wnl_complete_algebraic_state( ...
    A,Bslow,state,[],completionOpts,rightHandSide);
information.accepted = false;
end

function [correction, flag, iterationCount, information] = ...
        regularized_forced_correction(A,residual,opts)
[scaledA,scaledResidual] = equilibrate_forced_equations(A,residual);
numberOfUnknowns = size(A,2);
regularization = opts.forcedCorrectionRegularization;
if regularization > 0
    objectiveA = [scaledA; ...
        sqrt(regularization)*speye(numberOfUnknowns)];
    objectiveB = [scaledResidual; ...
        complex(zeros(numberOfUnknowns,1))];
else
    objectiveA = scaledA;
    objectiveB = scaledResidual;
end
[objectiveA,columnScale] = equilibrate_forced_columns(objectiveA);
information = struct('available',false,'equationResidual',Inf, ...
    'reportedRelativeResidual',NaN,'regularization',regularization);
correction = complex(zeros(numberOfUnknowns,1));
flag = 2;
iterationCount = 0;
if opts.forcedSolveMaxIterations <= 0
    return;
end
try
    [scaledCorrection,flag,reportedResidual,iteration] = lsqr( ...
        objectiveA,objectiveB,opts.forcedSolveRefinementTolerance, ...
        opts.forcedSolveMaxIterations);
    correction = columnScale.*scaledCorrection;
    iterationCount = forced_iteration_count(iteration);
    information.available = all(isfinite(real(correction))) && ...
        all(isfinite(imag(correction)));
    information.equationResidual = norm(A*correction-residual) / ...
        max(norm(residual),eps);
    information.reportedRelativeResidual = reportedResidual;
catch solverError
    information.errorIdentifier = solverError.identifier;
    information.errorMessage = solverError.message;
end
end

function [candidate,bestResidual,stepLength] = best_forced_line_search( ...
        A,b,state,correction,currentResidual,scale,maximumCuts)
candidate = state;
bestResidual = currentResidual;
stepLength = 0;
for cut = 0:maximumCuts
    trialStep = 2^(-cut);
    trial = state+trialStep*correction;
    trialResidual = norm(b-A*trial)/scale;
    if isfinite(trialResidual) && trialResidual < bestResidual
        candidate = trial;
        bestResidual = trialResidual;
        stepLength = trialStep;
    end
end
end

function [scaledA,scaledB,rowScale] = equilibrate_forced_equations(A,b)
rowNorm = sqrt(full(sum(abs(A).^2,2)));
largest = max(rowNorm);
if isempty(largest) || largest <= eps
    rowScale = ones(size(rowNorm));
else
    floorValue = max(largest*1.0e-12,eps);
    rowScale = ones(size(rowNorm));
    positive = rowNorm > 0;
    rowScale(positive) = 1./max(rowNorm(positive),floorValue);
end
scaledA = spdiags(rowScale,0,size(A,1),size(A,1))*A;
scaledB = rowScale.*b;
end

function [scaledA,columnScale] = equilibrate_forced_columns(A)
columnNorm = sqrt(full(sum(abs(A).^2,1))).';
largest = max(columnNorm);
columnScale = ones(size(columnNorm));
if isempty(largest) || largest <= eps
    % Leave structurally zero columns unchanged.  Assigning them the
    % reciprocal of a tiny floor creates an artificial O(1/eps) scale and
    % can destroy a descriptor preconditioner without changing A.
else
    floorValue = max(largest*1.0e-12,eps);
    positive = columnNorm > 0;
    columnScale(positive) = ...
        1./max(columnNorm(positive),floorValue);
end
scaledA = A*spdiags(columnScale,0,size(A,2),size(A,2));
end

function count = forced_iteration_count(iteration)
if numel(iteration) > 1
    count = iteration(1)*iteration(2);
else
    count = iteration;
end
end

function information = disabled_completion_information(reason)
information = struct('attempted',false,'available',false, ...
    'valid',false,'accepted',false,'stopReason',reason, ...
    'initialResidualNorm',NaN,'finalResidualNorm',NaN, ...
    'initialFullNorm',NaN,'finalFullNorm',NaN);
end

function information = disabled_rank_information(reason)
information = struct('attempted',false,'available',false, ...
    'numberOfAttempts',0,'selectedAttempt',0, ...
    'selectedTolerance',NaN,'selectedMethod','', ...
    'relativeResidual',NaN,'passedPhysicalGate',false, ...
    'fullToDescriptorNormRatio',NaN, ...
    'attemptTolerances',zeros(0,1),'attemptMethods',{{}}, ...
    'attemptResiduals',zeros(0,1),'attemptFullNorms',zeros(0,1), ...
    'attemptDescriptorNorms',zeros(0,1), ...
    'attemptFullToDescriptorRatios',zeros(0,1), ...
    'attemptFinite',false(0,1),'columnScaleRange',[NaN,NaN], ...
    'stopReason',reason);
end

function information = disabled_model_solve_information(reason)
information = struct('attempted',false,'available',false, ...
    'method','none','passedPhysicalGate',false, ...
    'relativeResidual',NaN,'stopReason',reason);
end

function value = model_solve_logical(information,fieldName)
value = isfield(information,fieldName) && ...
    logical(information.(fieldName));
end

function value = model_solve_numeric(information,fieldName)
if isfield(information,fieldName) && ...
        isnumeric(information.(fieldName)) && ...
        isscalar(information.(fieldName))
    value = information.(fieldName);
else
    value = NaN;
end
end

function value = model_solve_text(information,fieldName,defaultValue)
if isfield(information,fieldName) && ...
        (ischar(information.(fieldName)) || ...
        (isstring(information.(fieldName)) && ...
        isscalar(information.(fieldName))))
    value = char(information.(fieldName));
else
    value = defaultValue;
end
end

function [accepted,diagnostic] = accept_forced_completion_seed( ...
        completion,forcingNorm,opts)
% Do not abandon a better raw field for a DAE completion that merely moves
% residual among algebraic rows. A temporarily worse completion is a valid
% refinement seed only when it substantially compacts a pressure/gauge-
% inflated state. This distinction is essential for ordinary O(10^2)
% second-order fields, where completion may otherwise increase the residual
% by two or three orders of magnitude without reducing the state norm.
validateattributes( ...
    opts.forcedAlgebraicCompletionMaximumFullNormRatio,{'numeric'}, ...
    {'scalar','real','positive','<=',1,'finite'});

initialResidual = completion.initialResidualNorm;
finalResidual = completion.finalResidualNorm;
initialFullNorm = completion.initialFullNorm;
finalFullNorm = completion.finalFullNorm;
fullNormRatio = finalFullNorm/max(initialFullNorm,eps);
maximumSeedResidual = max( ...
    opts.forcedAlgebraicCompletionMaximumResidualGrowth* ...
        initialResidual,forcingNorm);
residualBounded = isfinite(finalResidual) && ...
    finalResidual <= maximumSeedResidual;
residualImproved = isfinite(finalResidual) && ...
    finalResidual <= initialResidual*(1+64*eps);
substantialCompaction = isfinite(fullNormRatio) && ...
    fullNormRatio <= ...
        opts.forcedAlgebraicCompletionMaximumFullNormRatio;
finiteInformation = completion.available && completion.valid && ...
    isfinite(initialResidual) && isfinite(initialFullNorm) && ...
    isfinite(finalFullNorm);
accepted = finiteInformation && residualBounded && ...
    (residualImproved || substantialCompaction);

if ~finiteInformation
    reason = 'rejected: completion is unavailable, invalid, or nonfinite';
elseif ~residualBounded
    reason = 'rejected: completed residual exceeds the seed-growth bound';
elseif residualImproved
    reason = 'accepted: completion does not increase the physical residual';
elseif substantialCompaction
    reason = 'accepted: pressure/gauge state was substantially compacted';
else
    reason = ['rejected: residual increased without substantial ', ...
        'pressure/gauge compaction'];
end
diagnostic = struct('residualBounded',residualBounded, ...
    'residualImproved',residualImproved, ...
    'substantialCompaction',substantialCompaction, ...
    'fullNormRatio',fullNormRatio, ...
    'maximumFullNormRatio', ...
        opts.forcedAlgebraicCompletionMaximumFullNormRatio, ...
    'maximumSeedResidual',maximumSeedResidual,'reason',reason);
end

function tolerance = forced_exploratory_tolerance(opts)
tolerance = opts.forcedSolveResidualTolerance;
if isempty(opts.forcedExploratoryResidualTolerance)
    return;
end
validateattributes(opts.forcedExploratoryResidualTolerance,{'numeric'}, ...
    {'scalar','real','finite','positive','>=', ...
    opts.forcedSolveResidualTolerance});
tolerance = opts.forcedExploratoryResidualTolerance;
end

function [acceptanceTolerance,refinementTolerance] = ...
        forced_refinement_tolerances(opts)
acceptanceTolerance = opts.forcedSolveResidualTolerance;
refinementTolerance = opts.forcedSolveRefinementTolerance;
% A looser ceiling is a computational early-stop only when the caller has
% explicitly elected to retain unconverged fields.  Strict runs continue to
% refine toward the original acceptance criteria.
if opts.stopOnUnconvergedForcedSolve || ...
        isempty(opts.forcedExploratoryResidualTolerance)
    return;
end
acceptanceTolerance = forced_exploratory_tolerance(opts);
refinementTolerance = max(refinementTolerance,acceptanceTolerance);
end

function [x, diagnostics] = refine_bordered_solution(A, b, x, opts)
% Resonant bordered systems retain the original correction solve.  Their
% extra projection coordinates are not primitive-variable DAE unknowns and
% therefore must not enter the algebraic-completion partition.
scale = max(norm(b), eps);
bestX = x;
bestRelativeResidual = norm(b-A*x)/scale;
totalIterations = 0;
lastFlag = double(bestRelativeResidual > ...
    opts.forcedSolveRefinementTolerance);
restarts = 0;
for attempt = 1:opts.forcedSolveMaxRestarts
    if bestRelativeResidual <= opts.forcedSolveRefinementTolerance
        break;
    end
    residual = b-A*x;
    [correction, flag, ~, iteration] = lsqr(A, residual, ...
        opts.forcedSolveRefinementTolerance, ...
        opts.forcedSolveMaxIterations);
    x = x + correction;
    totalIterations = totalIterations+forced_iteration_count(iteration);
    relativeResidual = norm(b-A*x)/scale;
    improved = relativeResidual < bestRelativeResidual*(1-64*eps);
    if improved
        bestX = x;
        bestRelativeResidual = relativeResidual;
    else
        x = bestX;
    end
    lastFlag = flag;
    restarts = attempt;
    if flag == 0 || ~improved
        break;
    end
end
x = bestX;
diagnostics = struct('flag', lastFlag, ...
    'totalIterations', totalIterations, 'restarts', restarts, ...
    'relativeResidual', bestRelativeResidual, ...
    'pressureSafeRefinement',false);
end
end

function x = minimum_norm_solve(A, b, opts)
n = max(size(A));
if n <= opts.fullSvdMax
    [u, s, v] = svd(full(A), 'econ');
    singularValues = diag(s);
    if isempty(singularValues)
        x = zeros(size(A, 2), 1);
        return;
    end
    keep = singularValues > opts.solveTolerance * max(singularValues);
    if ~any(keep)
        x = complex(zeros(size(A, 2), 1));
    else
        x = v(:, keep) * ((u(:, keep)' * b) ./ singularValues(keep));
    end
    return;
end

if exist('lsqminnorm', 'file') == 2
    x = lsqminnorm(A, b, opts.solveTolerance);
else
    x = A \ b;
end
end
