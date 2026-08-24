function [x, information] = wnl_rank_aware_minimum_norm( ...
        A,b,Bslow,userOpts)
%WNL_RANK_AWARE_MINIMUM_NORM Solve a forced block without hiding rank loss.
%
% The primitive-variable cylinder operator mixes velocity, interface,
% pressure, gauge, and constraint columns with very different magnitudes.
% Passing the generic iterative-solver tolerance directly to LSQMINNORM as
% an absolute QR-rank tolerance can therefore discard weak but physically
% necessary directions.  This routine first equilibrates COLUMNS (which
% preserves the physical least-squares norm), then tries progressively
% smaller explicit rank tolerances.  Every candidate is selected using the
% original, unequilibrated residual norm ||A*x-b||/||b||.
%
% No projection of b is made and the physical acceptance tolerance is not
% relaxed.  If no candidate reaches that tolerance, the best physical
% candidate is returned and the caller must keep the forced field invalid.

if nargin < 3 || isempty(Bslow)
    Bslow = speye(size(A,2));
end
if nargin < 4
    userOpts = struct();
end
opts = wnl_options(userOpts);

validateattributes(opts.forcedRankToleranceFactors,{'numeric'}, ...
    {'vector','real','positive','finite'});
validateattributes(opts.forcedColumnEquilibrationFloorRatio,{'numeric'}, ...
    {'scalar','real','positive','<=',1,'finite'});

scale = max(norm(b),eps);
[equilibratedA,columnScale] = equilibrate_columns( ...
    A,opts.forcedColumnEquilibrationFloorRatio);
tolerances = opts.solveTolerance * ...
    opts.forcedRankToleranceFactors(:);
tolerances = unique(tolerances,'stable');
maximumAttempts = numel(tolerances) + ...
    double(opts.forcedTryDefaultRankTolerance);

residuals = nan(maximumAttempts,1);
fullNorms = nan(maximumAttempts,1);
descriptorNorms = nan(maximumAttempts,1);
ratios = nan(maximumAttempts,1);
usedTolerances = nan(maximumAttempts,1);
methods = cell(maximumAttempts,1);
candidateFinite = false(maximumAttempts,1);

x = complex(zeros(size(A,2),1));
bestResidual = Inf;
bestRatio = Inf;
selectedAttempt = 0;
attempt = 0;

useDenseSvd = max(size(equilibratedA)) <= opts.fullSvdMax;
if useDenseSvd
    [u,s,v] = svd(full(equilibratedA),'econ');
    singularValues = diag(s);
end

for toleranceIndex = 1:numel(tolerances)
    attempt = attempt+1;
    tolerance = tolerances(toleranceIndex);
    if useDenseSvd
        scaledCandidate = truncated_svd_solution( ...
            u,singularValues,v,b,tolerance);
        method = 'column-equilibrated SVD';
    elseif exist('lsqminnorm','file') == 2
        scaledCandidate = lsqminnorm(equilibratedA,b,tolerance);
        method = 'column-equilibrated lsqminnorm';
    else
        scaledCandidate = equilibratedA\b;
        method = 'column-equilibrated backslash';
    end
    candidate = columnScale.*scaledCandidate;
    [x,bestResidual,bestRatio,selectedAttempt,passed] = ...
        consider_candidate(A,b,Bslow,candidate,attempt,x, ...
        bestResidual,bestRatio,selectedAttempt,scale,opts);
    residuals(attempt) = passed.residual;
    fullNorms(attempt) = passed.fullNorm;
    descriptorNorms(attempt) = passed.descriptorNorm;
    ratios(attempt) = passed.ratio;
    candidateFinite(attempt) = passed.finite;
    usedTolerances(attempt) = tolerance;
    methods{attempt} = method;
    if passed.residual <= opts.forcedSolveResidualTolerance
        break;
    end
    if ~useDenseSvd && exist('lsqminnorm','file') ~= 2
        break;
    end
end

if bestResidual > opts.forcedSolveResidualTolerance && ...
        opts.forcedTryDefaultRankTolerance && ...
        exist('lsqminnorm','file') == 2
    attempt = attempt+1;
    scaledCandidate = lsqminnorm(equilibratedA,b);
    candidate = columnScale.*scaledCandidate;
    [x,bestResidual,bestRatio,selectedAttempt,passed] = ...
        consider_candidate(A,b,Bslow,candidate,attempt,x, ...
        bestResidual,bestRatio,selectedAttempt,scale,opts);
    residuals(attempt) = passed.residual;
    fullNorms(attempt) = passed.fullNorm;
    descriptorNorms(attempt) = passed.descriptorNorm;
    ratios(attempt) = passed.ratio;
    candidateFinite(attempt) = passed.finite;
    usedTolerances(attempt) = NaN;
    methods{attempt} = 'column-equilibrated lsqminnorm default rank';
end

information = struct();
information.attempted = attempt > 0;
information.available = selectedAttempt > 0;
information.numberOfAttempts = attempt;
information.selectedAttempt = selectedAttempt;
information.selectedTolerance = NaN;
information.selectedMethod = '';
if selectedAttempt > 0
    information.selectedTolerance = usedTolerances(selectedAttempt);
    information.selectedMethod = methods{selectedAttempt};
end
information.relativeResidual = bestResidual;
information.passedPhysicalGate = isfinite(bestResidual) && ...
    bestResidual <= opts.forcedSolveResidualTolerance;
information.fullToDescriptorNormRatio = bestRatio;
information.attemptTolerances = usedTolerances(1:attempt);
information.attemptMethods = methods(1:attempt);
information.attemptResiduals = residuals(1:attempt);
information.attemptFullNorms = fullNorms(1:attempt);
information.attemptDescriptorNorms = descriptorNorms(1:attempt);
information.attemptFullToDescriptorRatios = ratios(1:attempt);
information.attemptFinite = candidateFinite(1:attempt);
information.columnScaleRange = positive_range(columnScale);
if information.passedPhysicalGate
    information.stopReason = ...
        'the original unequilibrated forced-field residual gate was met';
elseif information.available
    information.stopReason = ...
        'rank sweep exhausted; returning the smallest physical residual';
else
    information.stopReason = ...
        'no finite rank-aware least-squares candidate was available';
end
end

function scaledCandidate = truncated_svd_solution( ...
        u,singularValues,v,b,tolerance)
if isempty(singularValues)
    scaledCandidate = complex(zeros(size(v,1),1));
    return;
end
keep = singularValues > tolerance;
if ~any(keep)
    scaledCandidate = complex(zeros(size(v,1),1));
else
    scaledCandidate = v(:,keep) * ...
        ((u(:,keep)'*b)./singularValues(keep));
end
end

function [bestX,bestResidual,bestRatio,selectedAttempt,diagnostic] = ...
        consider_candidate(A,b,Bslow,candidate,attempt,bestX, ...
        bestResidual,bestRatio,selectedAttempt,scale,opts)
finiteCandidate = all(isfinite(real(candidate))) && ...
    all(isfinite(imag(candidate)));
if finiteCandidate
    residual = norm(A*candidate-b)/scale;
    fullNorm = norm(candidate);
    descriptorNorm = norm(Bslow*candidate);
    ratio = fullNorm/max(descriptorNorm,eps);
else
    residual = Inf;
    fullNorm = Inf;
    descriptorNorm = NaN;
    ratio = Inf;
end

strictlyBetter = isfinite(residual) && ...
    residual < bestResidual*(1-64*eps);
sameResidual = isfinite(residual) && isfinite(bestResidual) && ...
    abs(residual-bestResidual) <= ...
        64*eps*max([1,residual,bestResidual]);
smallerPhysicalState = sameResidual && ratio < bestRatio;
if finiteCandidate && (selectedAttempt == 0 || ...
        strictlyBetter || smallerPhysicalState)
    bestX = candidate;
    bestResidual = residual;
    bestRatio = ratio;
    selectedAttempt = attempt;
end

diagnostic = struct('finite',finiteCandidate,'residual',residual, ...
    'fullNorm',fullNorm,'descriptorNorm',descriptorNorm,'ratio',ratio, ...
    'physicalGate',residual <= opts.forcedSolveResidualTolerance);
end

function [scaledA,columnScale] = equilibrate_columns(A,floorRatio)
columnNorm = sqrt(full(sum(abs(A).^2,1))).';
largest = max(columnNorm);
columnScale = ones(size(columnNorm));
if ~isempty(largest) && largest > eps
    floorValue = max(largest*floorRatio,eps);
    positive = columnNorm > 0;
    columnScale(positive) = ...
        1./max(columnNorm(positive),floorValue);
end
scaledA = A*spdiags(columnScale,0,size(A,2),size(A,2));
end

function range = positive_range(values)
positive = values(isfinite(values) & values > 0);
if isempty(positive)
    range = [NaN,NaN];
else
    range = [min(positive),max(positive)];
end
end
