function report = wnl_forced_residual_report(solution, block, rowInformation, numberToPrint)
%WNL_FORCED_RESIDUAL_REPORT Locate an unconverged forced equation residual.
if nargin < 4 || isempty(numberToPrint)
    numberToPrint = 5;
end
metadata = [];
if isstruct(rowInformation) && isfield(rowInformation,'layout')
    metadata = rowInformation;
    layout = rowInformation.layout;
else
    layout = rowInformation;
end
residual = block.A*solution.vector-solution.forcing;
if ~isempty(solution.lambda)
    % Bordered corrections are already included in equationResidual but the
    % associated neutral vectors are not stored here. Report unbordered rows.
    report.note = 'unbordered residual shown';
else
    report.note = 'complete equation residual shown';
end
names = {'dense ur','dense utheta','dense w','dense p', ...
    'light ur','light utheta','light w','light p','zeta'};
local = {layout.d.ur,layout.d.ut,layout.d.w,layout.d.p, ...
    layout.l.ur,layout.l.ut,layout.l.w,layout.l.p,layout.zeta};
offsets = (0:numel(solution.spec.n)-1)*solution.spec.ndof;
families = repmat(struct('name','','norm',0,'fraction',0, ...
    'peakFrequency',0),numel(names),1);
total = max(norm(residual),eps);
for familyIndex = 1:numel(names)
    indices = reshape(local{familyIndex}(:)+offsets,[],1);
    byHarmonic = zeros(numel(offsets),1);
    for harmonicIndex = 1:numel(offsets)
        rows = local{familyIndex}(:)+offsets(harmonicIndex);
        byHarmonic(harmonicIndex) = norm(residual(rows));
    end
    [~,peak] = max(byHarmonic);
    families(familyIndex).name = names{familyIndex};
    families(familyIndex).norm = norm(residual(indices));
    families(familyIndex).fraction = families(familyIndex).norm/total;
    families(familyIndex).peakFrequency = ...
        solution.spec.n(peak)+solution.spec.s;
end
[~,order] = sort([families.norm],'descend');
numberToPrint = min(numberToPrint,numel(order));
fprintf('  Largest forced residual row families for %s\n',solution.spec.label);
for rank = 1:numberToPrint
    family = families(order(rank));
    fprintf('    %.3e  %5.1f%%  n+s=%+.3g  %s rows\n', ...
        family.norm,100*family.fraction,family.peakFrequency,family.name);
end
report.total = norm(residual);
report.families = families;
if ~isempty(metadata)
    report.largestRows = print_largest_rows( ...
        residual,solution.spec,layout,metadata,numberToPrint);
else
    report.largestRows = struct([]);
end
end

function entries = print_largest_rows( ...
        residual,spec,layout,metadata,numberToPrint)
[~,order] = sort(abs(residual),'descend');
numberToPrint = min(numberToPrint,numel(order));
entries = repmat(struct('magnitude',0,'frequency',0,'localRow',0, ...
    'description',''),numberToPrint,1);
fprintf('  Largest individual forced-residual equations\n');
for rank = 1:numberToPrint
    globalRow = order(rank);
    harmonicPosition = floor((globalRow-1)/spec.ndof)+1;
    localRow = globalRow-(harmonicPosition-1)*spec.ndof;
    description = describe_cylinder_row( ...
        localRow,layout,metadata,spec.m);
    entries(rank).magnitude = abs(residual(globalRow));
    entries(rank).frequency = ...
        spec.n(harmonicPosition)+spec.s;
    entries(rank).localRow = localRow;
    entries(rank).description = description;
    fprintf('    %.3e  n+s=%+.3g  local row %d: %s\n', ...
        entries(rank).magnitude,entries(rank).frequency, ...
        localRow,description);
end
end

function description = describe_cylinder_row(localRow,layout,metadata,m)
families = {'ur','ut','w','p'};
layerNames = {'d','l'};
fluidNames = {'dense','light'};
for layerIndex = 1:2
    layerName = layerNames{layerIndex};
    layer = layout.(layerName);
    for familyIndex = 1:numel(families)
        familyName = families{familyIndex};
        position = find(layer.(familyName) == localRow,1);
        if isempty(position)
            continue;
        end
        ir = mod(position-1,layout.nr)+1;
        iz = floor((position-1)/layout.nr)+1;
        role = cylinder_row_role( ...
            layerName,familyName,ir,iz,layout,metadata,m);
        description = sprintf('%s %s, ir=%d, iz=%d (%s)', ...
            fluidNames{layerIndex},familyName,ir,iz,role);
        return;
    end
end
zetaPosition = find(layout.zeta == localRow,1);
if ~isempty(zetaPosition)
    if zetaPosition == 1
        role = 'axis regularity';
    elseif zetaPosition == layout.nr
        role = 'contact line';
    elseif m == 0 && zetaPosition == ...
            max(2,min(layout.nr-1,ceil(layout.nr/2)))
        role = 'volume constraint';
    else
        role = 'kinematic condition';
    end
    description = sprintf('zeta, ir=%d (%s)',zetaPosition,role);
    return;
end
description = 'unclassified cylinder row';
end

function role = cylinder_row_role( ...
        layerName,familyName,ir,iz,layout,metadata,m)
if strcmp(layerName,'d')
    nz = layout.nzD;
    physicalWallIz = 1;
    physicalInterfaceIz = nz;
    if isfield(metadata,'verticalGrid') && ...
            isfield(metadata.verticalGrid,'lower')
        interfacePairs = metadata.verticalGrid.lower.interfacePairs;
    else
        interfacePairs = zeros(0,2);
    end
else
    nz = layout.nzL;
    physicalWallIz = nz;
    physicalInterfaceIz = 1;
    if isfield(metadata,'verticalGrid') && ...
            isfield(metadata.verticalGrid,'upper')
        interfacePairs = metadata.verticalGrid.upper.interfacePairs;
    else
        interfacePairs = zeros(0,2);
    end
end
if ir == 1
    role = 'axis regularity';
    return;
elseif ir == layout.nr && ~strcmp(familyName,'p')
    role = 'sidewall condition';
    return;
elseif ir > 1 && ir < layout.nr && ...
        iz == physicalWallIz && ~strcmp(familyName,'p')
    role = 'horizontal-wall condition';
    return;
elseif ir > 1 && ir < layout.nr && iz == physicalInterfaceIz
    if strcmp(layerName,'d') && ~strcmp(familyName,'p')
        role = 'interface velocity continuity';
        return;
    elseif strcmp(layerName,'l') && ~strcmp(familyName,'p')
        role = 'interface traction';
        return;
    end
end
leftCopy = any(interfacePairs(:,1) == iz);
rightCopy = any(interfacePairs(:,2) == iz);
if ir > 1 && ir < layout.nr && leftCopy
    role = 'vertical-element C0 matching';
elseif ir > 1 && ir < layout.nr && rightCopy && ...
        ~strcmp(familyName,'p')
    role = 'vertical-element C1 matching';
elseif strcmp(familyName,'p')
    if m == 0
        role = 'incompressibility or pressure gauge';
    else
        role = 'incompressibility';
    end
else
    role = 'bulk momentum';
end
end
