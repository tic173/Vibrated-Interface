function diagnostic = wnl_reference_eigenvalue_consistency(mode,userOpts)
%WNL_REFERENCE_EIGENVALUE_CONSISTENCY Compare reduced and full exponents.
%
% The prescribed-interface seed carries the reduced-solver exponent. A
% full-state Newton correction is allowed and expected, but a large shift
% means that the two discretizations do not yet represent the same linear
% operating-point mode. Passing A(lambda)*phi=0 is then an algebraic result,
% not sufficient evidence for a quantitative linear/WNL comparison.

if nargin < 2
    userOpts = struct();
end
opts = wnl_options(userOpts);
diagnostic = struct('applicable',false,'consistent',true, ...
    'initialLambda',NaN,'finalLambda',wnl_spec_lambda(mode.spec), ...
    'absoluteMismatch',0,'relativeMismatch',0, ...
    'allowedMismatch',Inf,'relativeTolerance', ...
        opts.maximumReducedFullEigenvalueRelativeMismatch, ...
    'absoluteTolerance', ...
        opts.maximumReducedFullEigenvalueAbsoluteMismatch);
initialLambda = [];
if isfield(mode,'tracking') && ...
        isfield(mode.tracking,'eigenpairRefinement')
    refinement = mode.tracking.eigenpairRefinement;
    if isstruct(refinement) && isfield(refinement,'attempted') && ...
            refinement.attempted && isfield(refinement,'initialLambda')
        initialLambda = refinement.initialLambda;
    end
end
% A cached full mode bypasses the expensive Newton reconstruction, but its
% current spec still carries the reduced operating-point reference used to
% establish branch agreement.  Preserve the same reduced/full gate instead
% of treating a reused mode as automatically consistent.
if isempty(initialLambda) && isfield(mode.spec,'reducedReferenceLambda') && ...
        ~isempty(mode.spec.reducedReferenceLambda)
    initialLambda = mode.spec.reducedReferenceLambda;
end
if isempty(initialLambda)
    return;
end
validateattributes(opts.maximumReducedFullEigenvalueRelativeMismatch, ...
    {'numeric'},{'scalar','real','nonnegative','finite'});
validateattributes(opts.maximumReducedFullEigenvalueAbsoluteMismatch, ...
    {'numeric'},{'scalar','real','nonnegative','finite'});
finalLambda = wnl_spec_lambda(mode.spec);
scale = max([abs(initialLambda),abs(finalLambda), ...
    opts.maximumReducedFullEigenvalueAbsoluteMismatch,eps]);
absoluteMismatch = abs(finalLambda-initialLambda);
allowedMismatch = max( ...
    opts.maximumReducedFullEigenvalueAbsoluteMismatch, ...
    opts.maximumReducedFullEigenvalueRelativeMismatch*scale);
diagnostic.applicable = true;
diagnostic.initialLambda = initialLambda;
diagnostic.finalLambda = finalLambda;
diagnostic.absoluteMismatch = absoluteMismatch;
diagnostic.relativeMismatch = absoluteMismatch/scale;
diagnostic.allowedMismatch = allowedMismatch;
diagnostic.consistent = absoluteMismatch <= allowedMismatch;
end
