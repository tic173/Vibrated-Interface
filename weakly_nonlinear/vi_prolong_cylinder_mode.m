function [seed,diagnostics] = vi_prolong_cylinder_mode( ...
        mode,sourceMetadata,targetR,targetVerticalGrid)
%VI_PROLONG_CYLINDER_MODE Interpolate a saved full Floquet mode in r and z.
%
% Multidomain vertical grids are interpolated element by element so the two
% copies of every C0 interface node remain distinct. The returned vector is
% only a recovery seed; WNL_COMPUTE_MODE rechecks it on the target operator.

if ~isstruct(mode) || ~isfield(mode,'field') || ...
        ~isfield(mode.field,'coeff')
    error('vi_prolong_cylinder_mode:BadMode', ...
        'mode.field.coeff is required.');
end
if ~isstruct(sourceMetadata) || ~isfield(sourceMetadata,'layout') || ...
        ~isfield(sourceMetadata,'radialGrid') || ...
        ~isfield(sourceMetadata,'verticalGrid')
    error('vi_prolong_cylinder_mode:BadMetadata', ...
        'Source layout, radial grid, and vertical grids are required.');
end
if ~isstruct(targetVerticalGrid) || ...
        ~isfield(targetVerticalGrid,'lower') || ...
        ~isfield(targetVerticalGrid,'upper')
    error('vi_prolong_cylinder_mode:BadTargetVerticalGrid', ...
        'Target lower and upper vertical grids are required.');
end

sourceLayout = sourceMetadata.layout;
sourceR = sourceMetadata.radialGrid.r(:);
targetR = targetR(:);
sourceNr = sourceLayout.nr;
targetNr = numel(targetR);
validate_grid(sourceR,sourceNr,'source radial');
validate_grid(targetR,targetNr,'target radial');

sourceCoeff = mode.field.coeff;
if size(sourceCoeff,1) ~= sourceLayout.ndof
    error('vi_prolong_cylinder_mode:CoefficientSize', ...
        'The saved coefficient rows do not match the source layout.');
end
numberOfHarmonics = size(sourceCoeff,2);
sourceLower = sourceMetadata.verticalGrid.lower;
sourceUpper = sourceMetadata.verticalGrid.upper;
targetLower = targetVerticalGrid.lower;
targetUpper = targetVerticalGrid.upper;
targetLayout = cylinder_layout(targetNr, ...
    numel(targetLower.z),numel(targetUpper.z));
radialInterpolation = barycentric_matrix(sourceR,targetR);
lowerInterpolation = vertical_interpolation(sourceLower,targetLower);
upperInterpolation = vertical_interpolation(sourceUpper,targetUpper);
targetCoeff = complex(zeros(targetLayout.ndof,numberOfHarmonics));

variables = {'ur','ut','w','p'};
for variableIndex = 1:numel(variables)
    name = variables{variableIndex};
    targetCoeff(targetLayout.d.(name),:) = interpolate_volume( ...
        sourceCoeff(sourceLayout.d.(name),:),radialInterpolation, ...
        lowerInterpolation,sourceNr,targetNr,sourceLayout.nzD, ...
        targetLayout.nzD,numberOfHarmonics);
    targetCoeff(targetLayout.l.(name),:) = interpolate_volume( ...
        sourceCoeff(sourceLayout.l.(name),:),radialInterpolation, ...
        upperInterpolation,sourceNr,targetNr,sourceLayout.nzL, ...
        targetLayout.nzL,numberOfHarmonics);
end
targetCoeff(targetLayout.zeta,:) = radialInterpolation* ...
    sourceCoeff(sourceLayout.zeta,:);

seed = targetCoeff(:);
diagnostics.sourceNr = sourceNr;
diagnostics.targetNr = targetNr;
diagnostics.sourceNzLower = sourceLayout.nzD;
diagnostics.targetNzLower = targetLayout.nzD;
diagnostics.sourceNzUpper = sourceLayout.nzL;
diagnostics.targetNzUpper = targetLayout.nzL;
diagnostics.numberOfHarmonics = numberOfHarmonics;
diagnostics.sourceNdof = sourceLayout.ndof;
diagnostics.targetNdof = targetLayout.ndof;
diagnostics.radialConstantResidual = constant_residual( ...
    radialInterpolation);
diagnostics.lowerConstantResidual = constant_residual( ...
    lowerInterpolation);
diagnostics.upperConstantResidual = constant_residual( ...
    upperInterpolation);
diagnostics.finite = all(isfinite(real(seed))) && ...
    all(isfinite(imag(seed)));
end

function values = interpolate_volume(coefficients,radialInterpolation, ...
        verticalInterpolation,sourceNr,targetNr,sourceNz,targetNz,nh)
sourceValues = reshape(coefficients,sourceNr,sourceNz,nh);
targetValues = complex(zeros(targetNr,targetNz,nh));
for harmonicIndex = 1:nh
    targetValues(:,:,harmonicIndex) = radialInterpolation* ...
        sourceValues(:,:,harmonicIndex)*verticalInterpolation.';
end
values = reshape(targetValues,targetNr*targetNz,nh);
end

function matrix = vertical_interpolation(sourceGrid,targetGrid)
sourceZ = sourceGrid.z(:);
targetZ = targetGrid.z(:);
if isfield(sourceGrid,'elementIndices') && ...
        isfield(targetGrid,'elementIndices') && ...
        numel(sourceGrid.elementIndices) == ...
        numel(targetGrid.elementIndices)
    matrix = zeros(numel(targetZ),numel(sourceZ));
    for elementIndex = 1:numel(sourceGrid.elementIndices)
        sourceIndex = sourceGrid.elementIndices{elementIndex};
        targetIndex = targetGrid.elementIndices{elementIndex};
        matrix(targetIndex,sourceIndex) = barycentric_matrix( ...
            sourceZ(sourceIndex),targetZ(targetIndex));
    end
else
    validate_grid(sourceZ,numel(sourceZ),'source vertical');
    validate_grid(targetZ,numel(targetZ),'target vertical');
    matrix = barycentric_matrix(sourceZ,targetZ);
end
end

function layout = cylinder_layout(nr,nzD,nzL)
nD = nr*nzD;
nL = nr*nzL;
next = 0;
layout.d.ur = next+(1:nD); next = next+nD;
layout.d.ut = next+(1:nD); next = next+nD;
layout.d.w = next+(1:nD); next = next+nD;
layout.d.p = next+(1:nD); next = next+nD;
layout.l.ur = next+(1:nL); next = next+nL;
layout.l.ut = next+(1:nL); next = next+nL;
layout.l.w = next+(1:nL); next = next+nL;
layout.l.p = next+(1:nL); next = next+nL;
layout.zeta = next+(1:nr); next = next+nr;
layout.ndof = next;
layout.nr = nr;
layout.nzD = nzD;
layout.nzL = nzL;
end

function matrix = barycentric_matrix(source,target)
numberOfSourceNodes = numel(source);
weights = ones(numberOfSourceNodes,1);
for index = 1:numberOfSourceNodes
    difference = source(index)-source([1:index-1,index+1:end]);
    weights(index) = 1/prod(difference);
end
matrix = zeros(numel(target),numberOfSourceNodes);
scale = max(max(abs(source)),1);
for index = 1:numel(target)
    [distance,match] = min(abs(target(index)-source));
    if distance <= 64*eps(scale)
        matrix(index,match) = 1;
    else
        row = weights./(target(index)-source);
        matrix(index,:) = (row/sum(row)).';
    end
end
end

function validate_grid(nodes,expectedSize,label)
if numel(nodes) ~= expectedSize || expectedSize < 2 || ...
        any(~isfinite(nodes)) || numel(unique(nodes)) ~= expectedSize
    error('vi_prolong_cylinder_mode:BadGrid', ...
        'The %s grid must contain distinct finite nodes.',label);
end
end

function value = constant_residual(matrix)
value = norm(matrix*ones(size(matrix,2),1)- ...
    ones(size(matrix,1),1),inf);
end
