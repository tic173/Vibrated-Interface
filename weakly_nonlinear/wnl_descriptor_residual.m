function [value, details] = wnl_descriptor_residual( ...
        A, Bslow, vector, side)
%WNL_DESCRIPTOR_RESIDUAL Scale a DAE eigenvector residual physically.
%
% VALUE = WNL_DESCRIPTOR_RESIDUAL(A,BSLOW,VECTOR,'direct') returns
%
%   ||A*VECTOR|| / (||A||_scale ||VECTOR(dynamic columns)||).
%
% The dynamic columns are the columns on which BSLOW is nonzero. Pressure,
% gauge, and other algebraic variables have zero BSLOW columns and therefore
% cannot enlarge the denominator. This is essential for primitive-variable
% incompressible pencils: ||A*q||/(||A||*||q||) can be made arbitrarily small
% by adding a huge algebraic pressure/null component even though the physical
% velocity/interface part of q has not converged.
%
% With SIDE='adjoint', the analogous metric is
%
%   ||A'*VECTOR|| / (||A||_scale ||VECTOR(dynamic rows)||),
%
% where the dynamic rows are the nonzero rows of BSLOW. DETAILS reports both
% the physical and full norms so algebraic inflation is visible in saved
% diagnostics.

if nargin < 4
    side = 'direct';
end
if ~(ischar(side) || (isstring(side) && isscalar(side)))
    error('wnl_descriptor_residual:SideType', ...
        'side must be ''direct'' or ''adjoint''.');
end
side = lower(char(side));
if size(A, 1) ~= size(A, 2) || ...
        ~isequal(size(Bslow), size(A))
    error('wnl_descriptor_residual:MatrixSize', ...
        'A and Bslow must be square matrices of the same size.');
end
vector = vector(:);

switch side
    case 'direct'
        if numel(vector) ~= size(A, 2)
            error('wnl_descriptor_residual:VectorSize', ...
                'The direct vector has the wrong number of entries.');
        end
        residualVector = A*vector;
        supportMagnitude = sqrt(full(sum(abs(Bslow).^2, 1))).';
        descriptorComponentNorm = norm(Bslow*vector);
        supportKind = 'nonzero Bslow columns';
    case 'adjoint'
        if numel(vector) ~= size(A, 1)
            error('wnl_descriptor_residual:VectorSize', ...
                'The adjoint vector has the wrong number of entries.');
        end
        residualVector = A'*vector;
        supportMagnitude = sqrt(full(sum(abs(Bslow).^2, 2)));
        descriptorComponentNorm = norm(Bslow'*vector);
        supportKind = 'nonzero Bslow rows';
    otherwise
        error('wnl_descriptor_residual:SideValue', ...
            'side must be ''direct'' or ''adjoint''.');
end

maximumSupport = max(supportMagnitude);
if isempty(maximumSupport) || ~isfinite(maximumSupport) || ...
        maximumSupport <= 0
    support = false(size(supportMagnitude));
    supportThreshold = NaN;
else
    % Structural zeros in Bslow are exact. The small tolerance only removes
    % roundoff-level fill while retaining light-fluid mass coefficients.
    supportThreshold = 32*eps(maximumSupport);
    support = supportMagnitude > supportThreshold;
end

physicalVectorNorm = norm(vector(support));
fullVectorNorm = norm(vector);
equationResidualNorm = norm(residualVector);
operatorScale = max(1.0, deterministic_matrix_scale(A));
if ~isfinite(physicalVectorNorm) || physicalVectorNorm <= eps || ...
        ~isfinite(equationResidualNorm)
    effectiveResidual = Inf;
    value = Inf;
else
    effectiveResidual = equationResidualNorm/physicalVectorNorm;
    value = effectiveResidual/operatorScale;
end

details = struct();
details.side = side;
details.supportKind = supportKind;
details.supportThreshold = supportThreshold;
details.numberOfPhysicalEntries = nnz(support);
details.numberOfAlgebraicEntries = numel(support)-nnz(support);
details.physicalVectorNorm = physicalVectorNorm;
details.fullVectorNorm = fullVectorNorm;
details.fullToPhysicalNormRatio = fullVectorNorm/max(physicalVectorNorm,eps);
details.descriptorComponentNorm = descriptorComponentNorm;
details.equationResidualNorm = equationResidualNorm;
details.operatorScale = operatorScale;
details.effectiveResidual = effectiveResidual;
details.scaledResidual = value;
end

function scale = deterministic_matrix_scale(A)
oneNorm = norm(A, 1);
infinityNorm = norm(A, inf);
scale = max(sqrt(max(oneNorm, 0)*max(infinityNorm, 0)), eps);
end
