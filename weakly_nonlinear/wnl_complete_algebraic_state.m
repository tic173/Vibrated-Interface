function [completedState, information] = wnl_complete_algebraic_state( ...
        A, Bslow, state, preservedDofs, userOpts, rightHandSide)
%WNL_COMPLETE_ALGEBRAIC_STATE Minimum-norm DAE algebraic completion.
%
% [Q,INFO] = WNL_COMPLETE_ALGEBRAIC_STATE(A,BSLOW,Q0,KEEP,OPTS)
% [Q,INFO] = WNL_COMPLETE_ALGEBRAIC_STATE(A,BSLOW,Q0,KEEP,OPTS,F)
% leaves every descriptor-supported variable and every index in KEEP fixed.
% It recomputes the remaining pressure, gauge, boundary, and constraint
% variables from
%
%   min_z || A_f*q_f + A_a*z - F ||_2^2 + alpha*||z||_2^2.
%
% F defaults to zero.  The optional nonzero right-hand side is used by the
% second-order forced-field solver; the original zero-right-hand-side form
% is the homogeneous eigenpair completion.
%
% The regularization is inserted in physical coordinates before column
% equilibration. Thus scaling can accelerate LSQR without making a very
% large pressure/null-space component artificially cheap.

if nargin < 4 || isempty(preservedDofs)
    preservedDofs = zeros(0,1);
end
if nargin < 5
    userOpts = struct();
end
if nargin < 6 || isempty(rightHandSide)
    rightHandSide = complex(zeros(size(A,1),1));
end
opts = wnl_options(userOpts);

if size(A,1) ~= size(A,2) || ~isequal(size(Bslow),size(A))
    error('wnl_complete_algebraic_state:MatrixSize', ...
        'A and Bslow must be square matrices of the same size.');
end
state = state(:);
rightHandSide = rightHandSide(:);
numberOfUnknowns = size(A,2);
if numel(state) ~= numberOfUnknowns
    error('wnl_complete_algebraic_state:StateSize', ...
        'The state vector has the wrong number of entries.');
end
if numel(rightHandSide) ~= size(A,1)
    error('wnl_complete_algebraic_state:RightHandSideSize', ...
        'The right-hand side has the wrong number of entries.');
end
preservedDofs = unique(preservedDofs(:));
if any(preservedDofs < 1) || any(preservedDofs > numberOfUnknowns) || ...
        any(preservedDofs ~= round(preservedDofs))
    error('wnl_complete_algebraic_state:PreservedDofs', ...
        'preservedDofs must contain valid integer state indices.');
end

validateattributes(opts.eigenpairAlgebraicCompletionRegularization, ...
    {'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(opts.eigenpairAlgebraicCompletionSolveTolerance, ...
    {'numeric'},{'scalar','real','positive','finite'});
validateattributes(opts.eigenpairAlgebraicCompletionEquationTolerance, ...
    {'numeric'},{'scalar','real','positive','<',1,'finite'});
validateattributes(opts.eigenpairAlgebraicCompletionMaxIterations, ...
    {'numeric'},{'scalar','integer','nonnegative'});
validateattributes(opts.eigenpairAlgebraicCompletionMaxRestarts, ...
    {'numeric'},{'scalar','integer','nonnegative'});

supportMagnitude = sqrt(full(sum(abs(Bslow).^2,1))).';
maximumSupport = max(supportMagnitude);
if isempty(maximumSupport) || ~isfinite(maximumSupport) || ...
        maximumSupport <= 0
    descriptorSupport = false(numberOfUnknowns,1);
    supportThreshold = NaN;
else
    supportThreshold = 32*eps(maximumSupport);
    descriptorSupport = supportMagnitude > supportThreshold;
end
fixed = descriptorSupport;
fixed(preservedDofs) = true;
fixedDofs = find(fixed);
algebraicDofs = find(~fixed);

initialResidualNorm = norm(A*state-rightHandSide);
initialFullNorm = norm(state);
initialDescriptorNorm = norm(Bslow*state);
completedState = state;
information = base_information();
information.supportThreshold = supportThreshold;
information.numberOfFixedDofs = numel(fixedDofs);
information.numberOfAlgebraicDofs = numel(algebraicDofs);
information.initialResidualNorm = initialResidualNorm;
information.finalResidualNorm = initialResidualNorm;
information.initialFullNorm = initialFullNorm;
information.finalFullNorm = initialFullNorm;
information.initialDescriptorNorm = initialDescriptorNorm;
information.finalDescriptorNorm = initialDescriptorNorm;

if isempty(algebraicDofs) || ...
        opts.eigenpairAlgebraicCompletionMaxIterations == 0
    information.stopReason = 'no enabled algebraic completion solve';
    return;
end

information.attempted = true;
algebraicMatrix = A(:,algebraicDofs);
algebraicRightHandSide = rightHandSide- ...
    A(:,fixedDofs)*state(fixedDofs);
[scaledMatrix,scaledRightHandSide,rowScale] = ...
    equilibrate_equations(algebraicMatrix,algebraicRightHandSide);

regularization = opts.eigenpairAlgebraicCompletionRegularization;
numberOfAlgebraicDofs = numel(algebraicDofs);
if regularization > 0
    objectiveMatrix = [scaledMatrix; ...
        sqrt(regularization)*speye(numberOfAlgebraicDofs)];
    objectiveRightHandSide = [scaledRightHandSide; ...
        complex(zeros(numberOfAlgebraicDofs,1))];
else
    objectiveMatrix = scaledMatrix;
    objectiveRightHandSide = scaledRightHandSide;
end
[objectiveMatrix,columnScale] = equilibrate_columns(objectiveMatrix);

scaledGuess = complex(zeros(numberOfAlgebraicDofs,1));
maximumAttempts = opts.eigenpairAlgebraicCompletionMaxRestarts+1;
bestState = [];
bestResidualNorm = Inf;
bestEquationResidual = Inf;
bestFlag = NaN;
bestRelativeResidual = NaN;
bestAttempt = NaN;
totalIterations = 0;
attemptResidualNorms = nan(maximumAttempts,1);
attemptEquationResiduals = nan(maximumAttempts,1);
attemptFlags = nan(maximumAttempts,1);
attemptIterations = nan(maximumAttempts,1);

for attempt = 1:maximumAttempts
    [scaledCandidate,flag,relativeResidual,iterations] = lsqr( ...
        objectiveMatrix,objectiveRightHandSide, ...
        opts.eigenpairAlgebraicCompletionSolveTolerance, ...
        opts.eigenpairAlgebraicCompletionMaxIterations,[],[], ...
        scaledGuess);
    totalIterations = totalIterations+iterations;
    algebraicCandidate = columnScale.*scaledCandidate;
    candidateState = state;
    candidateState(algebraicDofs) = algebraicCandidate;
    candidateResidualNorm = norm(A*candidateState-rightHandSide);
    candidateEquationResidual = candidateResidualNorm / ...
        max(norm(algebraicRightHandSide),eps);

    attemptResidualNorms(attempt) = candidateResidualNorm;
    attemptEquationResiduals(attempt) = candidateEquationResidual;
    attemptFlags(attempt) = flag;
    attemptIterations(attempt) = iterations;
    finiteCandidate = all(isfinite(real(candidateState))) && ...
        all(isfinite(imag(candidateState))) && ...
        isfinite(candidateResidualNorm);
    if finiteCandidate && candidateResidualNorm < bestResidualNorm
        bestState = candidateState;
        bestResidualNorm = candidateResidualNorm;
        bestEquationResidual = candidateEquationResidual;
        bestFlag = flag;
        bestRelativeResidual = relativeResidual;
        bestAttempt = attempt;
    end
    scaledGuess = scaledCandidate;
    if flag == 0 || candidateEquationResidual <= ...
            opts.eigenpairAlgebraicCompletionEquationTolerance
        break;
    end
end

numberOfAttempts = find(~isnan(attemptFlags),1,'last');
if isempty(numberOfAttempts)
    numberOfAttempts = 0;
end
information.numberOfAttempts = numberOfAttempts;
information.totalIterations = totalIterations;
information.attemptResidualNorms = ...
    attemptResidualNorms(1:numberOfAttempts);
information.attemptEquationResiduals = ...
    attemptEquationResiduals(1:numberOfAttempts);
information.attemptFlags = attemptFlags(1:numberOfAttempts);
information.attemptIterations = attemptIterations(1:numberOfAttempts);
information.rowScaleRange = positive_range(rowScale);
information.columnScaleRange = positive_range(columnScale);
information.regularization = regularization;

if isempty(bestState)
    information.stopReason = 'no finite algebraic completion candidate';
    return;
end

completedState = bestState;
information.available = true;
information.selectedAttempt = bestAttempt;
information.flag = bestFlag;
information.relativeResidual = bestRelativeResidual;
information.equationResidual = bestEquationResidual;
information.finalResidualNorm = bestResidualNorm;
information.finalFullNorm = norm(completedState);
information.finalDescriptorNorm = norm(Bslow*completedState);
information.fixedStateChange = norm( ...
    completedState(fixedDofs)-state(fixedDofs));
information.residualReductionFactor = bestResidualNorm / ...
    max(initialResidualNorm,eps);
information.fullNormReductionFactor = information.finalFullNorm / ...
    max(initialFullNorm,eps);
information.descriptorChange = norm( ...
    Bslow*(completedState-state));
information.valid = information.fixedStateChange <= ...
    100*eps(max(norm(state(fixedDofs)),1)) && ...
    information.descriptorChange <= ...
    100*eps(max(initialDescriptorNorm,1));
if information.valid
    information.stopReason = 'minimum-norm algebraic completion available';
else
    information.stopReason = ...
        'algebraic completion changed a preserved descriptor variable';
end
end

function information = base_information()
information = struct('attempted',false,'available',false,'valid',false, ...
    'stopReason','','numberOfFixedDofs',0,'numberOfAlgebraicDofs',0, ...
    'numberOfAttempts',0,'totalIterations',0,'selectedAttempt',NaN, ...
    'flag',NaN,'relativeResidual',NaN,'equationResidual',Inf, ...
    'regularization',NaN,'initialResidualNorm',NaN, ...
    'finalResidualNorm',NaN,'initialFullNorm',NaN, ...
    'finalFullNorm',NaN,'initialDescriptorNorm',NaN, ...
    'finalDescriptorNorm',NaN,'fixedStateChange',NaN, ...
    'descriptorChange',NaN,'residualReductionFactor',NaN, ...
    'fullNormReductionFactor',NaN,'supportThreshold',NaN, ...
    'rowScaleRange',[NaN,NaN],'columnScaleRange',[NaN,NaN], ...
    'attemptResidualNorms',zeros(0,1), ...
    'attemptEquationResiduals',zeros(0,1), ...
    'attemptFlags',zeros(0,1),'attemptIterations',zeros(0,1));
end

function [scaledA,scaledB,rowScale] = equilibrate_equations(A,b)
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

function [scaledA,columnScale] = equilibrate_columns(A)
columnNorm = sqrt(full(sum(abs(A).^2,1))).';
largest = max(columnNorm);
if isempty(largest) || largest <= eps
    columnScale = ones(size(columnNorm));
else
    floorValue = max(largest*1.0e-12,eps);
    columnScale = 1./max(columnNorm,floorValue);
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
