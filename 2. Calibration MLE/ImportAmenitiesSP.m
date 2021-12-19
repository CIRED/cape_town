%% Script for the import of the amenity data at the SP level

% Import of the amenity files at the SP level (for the regression)
importfile([path_nedum,'SP_amenities.csv']);

% Airport cones
airportCone = airport_cone;
airportCone(airport_cone==55) = 1;
airportCone(airport_cone==60) = 1;
airportCone(airport_cone==65) = 1;    
airportCone(airport_cone==70) = 1;
airportCone(airport_cone==75) = 1;
    
% Distance to RDP houses
distanceRDP = 2;
if distanceRDP ~= 2 
    matrixDistance = ((repmat(grille.xCoord, length(data.spX), 1) - repmat(data.spX, 1, length(grille.coord_horiz))).^2 ... 
                        + (repmat(grille.yCoord, length(data.spY), 1) - repmat(data.spY, 1, length(grille.yCoord))).^2)...
                        < distanceRDP ^ 2;
    SP_distance_RDP = (double(land.numberPropertiesRDP > 5) * (matrixDistance)')' > 1;
else
    load(strcat('.', slash, '0. Precalculated inputs', slash, 'SPdistanceRDP'))
end

tableAmenities = table(SP_CODE, distance_distr_parks < 2, ...
        distance_ocean < 2, distance_ocean > 2 & distance_ocean < 4,...
        distance_world_herit < 2, distance_world_herit > 2 & distance_world_herit < 4, ...
        distance_urban_herit < 2, distance_UCT < 2,...
        airportCone, slope > 1 & slope < 5, slope > 5, ...
        distance_train < 2, distance_protected_envir < 2, ..., 
        distance_protected_envir > 2 & distance_protected_envir < 4,...
        SP_distance_RDP, distance_power_station < 2, distance_biosphere_reserve < 2);
tableAmenities.Properties.VariableNames = {'SP_CODE' 'distance_distr_parks' 'distance_ocean' 'distance_ocean_2_4' 'distance_world_herit' 'distance_world_herit_2_4' 'distance_urban_herit' 'distance_UCT' 'airport_cone2' 'slope_1_5' 'slope_5' 'distance_train' 'distance_protected_envir' 'distance_protected_envir_2_4' 'RDP_proximity' 'distance_power_station' ,'distance_biosphere_reserve'};

% reorder tableAmenities to match order from data.XXX
[~,idx] = ismember(data.spCode, SP_CODE);
tableAmenities = tableAmenities(idx,:);

