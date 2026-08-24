function grid = vi_cylinder_vertical_grid(parameters, layerName)
%VI_CYLINDER_VERTICAL_GRID Single- or multi-domain cylinder z grid.
%
% Multi-domain defaults create thin Chebyshev elements at both the physical
% interface and rigid wall. Their width is based on the layer's oscillatory
% Stokes thickness. Explicit breaks/point counts override the automatic
% construction.

layerName = validatestring(layerName, {'lower', 'upper'});
numerics = parameters.numerics;
rhoUpperOverLower = (1-parameters.At)/(1+parameters.At);
if strcmp(layerName, 'lower')
    interval = [-1, 0];
    kinematicCoefficient = parameters.C;
    singlePointField = 'NzLower';
    breakField = 'lowerBreaks';
    pointField = 'lowerPoints';
else
    interval = [0, 1];
    kinematicCoefficient = ...
        parameters.eta*parameters.C/rhoUpperOverLower;
    singlePointField = 'NzUpper';
    breakField = 'upperBreaks';
    pointField = 'upperPoints';
end
stokesThickness = sqrt(2*kinematicCoefficient/parameters.omegaStar);

verticalOptions = struct();
verticalOptions.type = 'single';
verticalOptions.boundaryLayerWidthFactor = 4;
verticalOptions.maximumBoundaryLayerFraction = 0.20;
verticalOptions.pointsPerBoundaryLayer = 11;
verticalOptions.pointsInBulk = 15;
verticalOptions.lowerBreaks = [];
verticalOptions.lowerPoints = [];
verticalOptions.upperBreaks = [];
verticalOptions.upperPoints = [];
if isfield(numerics, 'verticalGrid') && ...
        ~isempty(numerics.verticalGrid)
    supplied = numerics.verticalGrid;
    names = fieldnames(supplied);
    for nameIndex = 1:numel(names)
        verticalOptions.(names{nameIndex}) = ...
            supplied.(names{nameIndex});
    end
end
gridType = validatestring(verticalOptions.type, ...
    {'single', 'multidomain'});

if strcmp(gridType, 'single')
    if ~isfield(numerics, singlePointField)
        error('vi_cylinder_vertical_grid:MissingSinglePointCount', ...
            ['numerics.%s is required when ', ...
             'verticalGrid.type=''single''.'], singlePointField);
    end
    singlePointCount = numerics.(singlePointField);
    validateattributes(singlePointCount, {'numeric'}, ...
        {'scalar', 'integer', '>=', 7});
    breaks = interval;
    points = singlePointCount;
    automatic = false;
else
    explicitBreaks = verticalOptions.(breakField);
    explicitPoints = verticalOptions.(pointField);
    if xor(isempty(explicitBreaks), isempty(explicitPoints))
        error('vi_cylinder_vertical_grid:IncompleteExplicitGrid', ...
            ['verticalGrid.%s and verticalGrid.%s must either both be ', ...
             'empty or both be supplied.'], breakField, pointField);
    end
    if ~isempty(explicitBreaks)
        breaks = explicitBreaks;
        points = explicitPoints;
        automatic = false;
    else
        validateattributes(verticalOptions.boundaryLayerWidthFactor, ...
            {'numeric'}, {'scalar', 'real', 'positive', 'finite'});
        validateattributes(verticalOptions.maximumBoundaryLayerFraction, ...
            {'numeric'}, {'scalar', 'real', 'positive', '<', 0.5, 'finite'});
        validateattributes(verticalOptions.pointsPerBoundaryLayer, ...
            {'numeric'}, {'scalar', 'integer', '>=', 3});
        validateattributes(verticalOptions.pointsInBulk, ...
            {'numeric'}, {'scalar', 'integer', '>=', 3});
        layerDepth = diff(interval);
        boundaryWidth = min( ...
            verticalOptions.boundaryLayerWidthFactor*stokesThickness, ...
            verticalOptions.maximumBoundaryLayerFraction*layerDepth);
        if strcmp(layerName, 'lower')
            breaks = [interval(1), interval(1)+boundaryWidth, ...
                interval(2)-boundaryWidth, interval(2)];
        else
            breaks = [interval(1), interval(1)+boundaryWidth, ...
                interval(2)-boundaryWidth, interval(2)];
        end
        points = [verticalOptions.pointsPerBoundaryLayer, ...
            verticalOptions.pointsInBulk, ...
            verticalOptions.pointsPerBoundaryLayer];
        automatic = true;
    end
end

breaks = breaks(:).';
if abs(breaks(1)-interval(1)) > 1.0e-12 || ...
        abs(breaks(end)-interval(2)) > 1.0e-12
    error('vi_cylinder_vertical_grid:WrongLayerEndpoints', ...
        'The %s grid must span exactly [%g,%g].', ...
        layerName, interval(1), interval(2));
end
grid = vi_multidomain_chebyshev_grid(breaks, points);
grid.type = gridType;
grid.layerName = layerName;
grid.automatic = automatic;
grid.stokesThickness = stokesThickness;
grid.kinematicCoefficient = kinematicCoefficient;
grid.wallCoordinate = interval(1);
grid.interfaceCoordinate = interval(2);
if strcmp(layerName, 'upper')
    grid.wallCoordinate = interval(2);
    grid.interfaceCoordinate = interval(1);
end

wallDistance = abs(grid.z-grid.wallCoordinate);
interfaceDistance = abs(grid.z-grid.interfaceCoordinate);
positiveWallDistance = wallDistance(wallDistance > 0);
positiveInterfaceDistance = interfaceDistance(interfaceDistance > 0);
grid.firstWallSpacing = min(positiveWallDistance);
grid.firstInterfaceSpacing = min(positiveInterfaceDistance);
grid.pointsInsideWallStokesLayer = ...
    nnz(wallDistance <= stokesThickness);
grid.pointsInsideInterfaceStokesLayer = ...
    nnz(interfaceDistance <= stokesThickness);
end
