function [seed,diagnostics] = vi_prolong_cylinder_mode_radially( ...
        mode,sourceMetadata,targetR)
%VI_PROLONG_CYLINDER_MODE_RADIALLY Interpolate a saved full Floquet mode.
%
% Every primitive variable and interface coefficient is interpolated in r;
% vertical and temporal indices are unchanged. The returned vector is only
% a recovery seed: WNL_COMPUTE_MODE must recheck/refine the target operator.

if ~isstruct(mode) || ~isfield(mode,'field') || ...
        ~isfield(mode.field,'coeff')
    error('vi_prolong_cylinder_mode_radially:BadMode', ...
        'mode.field.coeff is required.');
end
if ~isstruct(sourceMetadata) || ~isfield(sourceMetadata,'layout') || ...
        ~isfield(sourceMetadata,'radialGrid') || ...
        ~isfield(sourceMetadata.radialGrid,'r')
    error('vi_prolong_cylinder_mode_radially:BadMetadata', ...
        'Source layout and radial-grid nodes are required.');
end

sourceLayout = sourceMetadata.layout;
sourceR = sourceMetadata.radialGrid.r(:);
targetR = targetR(:);
sourceNr = sourceLayout.nr;
targetNr = numel(targetR);
if numel(sourceR) ~= sourceNr
    error('vi_prolong_cylinder_mode_radially:SourceGridSize', ...
        'The source radial nodes and layout have inconsistent sizes.');
end
if targetNr < 2 || any(~isfinite(targetR)) || ...
        numel(unique(targetR)) ~= targetNr
    error('vi_prolong_cylinder_mode_radially:TargetGrid', ...
        'targetR must contain at least two distinct finite nodes.');
end

sourceCoeff = mode.field.coeff;
if size(sourceCoeff,1) ~= sourceLayout.ndof
    error('vi_prolong_cylinder_mode_radially:CoefficientSize', ...
        'The saved coefficient rows do not match the source layout.');
end
numberOfHarmonics = size(sourceCoeff,2);
targetLayout = cylinder_layout( ...
    targetNr,sourceLayout.nzD,sourceLayout.nzL);
interpolation = barycentric_matrix(sourceR,targetR);
targetCoeff = complex(zeros(targetLayout.ndof,numberOfHarmonics));

variables = {'ur','ut','w','p'};
for variableIndex = 1:numel(variables)
    name = variables{variableIndex};
    targetCoeff(targetLayout.d.(name),:) = interpolate_block( ...
        sourceCoeff(sourceLayout.d.(name),:),interpolation, ...
        sourceNr,targetNr,sourceLayout.nzD,numberOfHarmonics);
    targetCoeff(targetLayout.l.(name),:) = interpolate_block( ...
        sourceCoeff(sourceLayout.l.(name),:),interpolation, ...
        sourceNr,targetNr,sourceLayout.nzL,numberOfHarmonics);
end
targetCoeff(targetLayout.zeta,:) = interpolation* ...
    sourceCoeff(sourceLayout.zeta,:);

seed = targetCoeff(:);
diagnostics.sourceNr = sourceNr;
diagnostics.targetNr = targetNr;
diagnostics.numberOfHarmonics = numberOfHarmonics;
diagnostics.sourceNdof = sourceLayout.ndof;
diagnostics.targetNdof = targetLayout.ndof;
diagnostics.interpolationConstantResidual = norm( ...
    interpolation*ones(sourceNr,1)-ones(targetNr,1),inf);
diagnostics.finite = all(isfinite(real(seed))) && ...
    all(isfinite(imag(seed)));
end

function values = interpolate_block( ...
        coefficients,interpolation,sourceNr,targetNr,nz,nh)
sourceValues = reshape(coefficients,sourceNr,nz,nh);
values = interpolation*reshape(sourceValues,sourceNr,[]);
values = reshape(values,targetNr*nz,nh);
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
