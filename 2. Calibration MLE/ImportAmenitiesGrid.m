%% Script for the import of the amenity data at the SP level

% Import of the amenity files at the grid level for the extrapolation
importfile([path_nedum,'grid_amenities.csv']);

% Airport cones
airportCone = airport_cone;
airportCone(airport_cone==55) = 1;
airportCone(airport_cone==60) = 1;
airportCone(airport_cone==65) = 1;
airportCone(airport_cone==70) = 1;
airportCone(airport_cone==75) = 1;

%Distance to RDP housing
if distanceRDP ~= 2
    matrixDistance = ((repmat(grid.xCoord,length(grid.xCoord),1) - repmat(grid.xCoord',1,length(grid.xCoord))).^2 ... 
                        + (repmat(grid.yCoord,length(grid.xCoord),1) - repmat(grid.yCoord',1,length(grid.xCoord))).^2)...
                        < distanceRDP ^ 2;
    gridDistanceRDP = (double(land.numberPropertiesRDP > 5) * (matrixDistance)')' > 1;
    % save(strcat('.', slash, '0. Precalculated inputs', slash, 'gridDistanceRDP'), 'gridDistanceRDP')
else 
    load(strcat('.', slash, '0. Precalculated inputs', slash, 'gridDistanceRDP'))
end

% Output as a table

tableAmenitiesGrid = table(distance_distr_parks < 2, ...
        distance_ocean < 2, distance_ocean > 2 & distance_ocean < 4,...
        distance_world_herit < 2, distance_world_herit > 2 & distance_world_herit < 4, ...
        distance_urban_herit < 2, distance_UCT < 2,...
        airportCone, slope > 1 & slope < 5, slope > 5, ...
        distance_train < 2, distance_protected_envir < 2, ..., 
        distance_protected_envir > 2 & distance_protected_envir < 4,...
        gridDistanceRDP, distance_power_station < 2, distance_biosphere_reserve < 2);
tableAmenitiesGrid.Properties.VariableNames = {'distance_distr_parks' 'distance_ocean' 'distance_ocean_2_4' 'distance_world_herit' 'distance_world_herit_2_4' 'distance_urban_herit' 'distance_UCT' 'airport_cone2' 'slope_1_5' 'slope_5' 'distance_train' 'distance_protected_envir' 'distance_protected_envir_2_4' 'RDP_proximity' 'distance_power_station' 'distance_biosphere_reserve'};


