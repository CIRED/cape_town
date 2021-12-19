function ExportDataForMaps(grid, param, land, data, poly, t, simulation, yearUrbanized, filename, scenarioName)
% Function that creates the export as csv that can then be imported in R to
% display simulation maps

finalState.dwellingSize(:,:) = simulation.dwellingSize(length(t),:, :);
finalState.rent(:,:) = simulation.rent(length(t),:, :);
finalState.households(:,:,:) = simulation.households(length(t),:, :, :);
finalState.householdsHousingType(:,:) = simulation.householdsHousingType(length(t),:,:);
finalState.householdsCenter(:,:) = simulation.householdsCenter(length(t), :, :);
finalState.housingSupply(:,:) = simulation.housingSupply(length(t),:, :);

% Extent of urban area for different dates
people(:,:) = sum(simulation.householdsHousingType, 2);
yearUrbanizedMat = (people > 10) .* repmat(simulation.T + param.yearBegin, 1, length(grid.distanceCBD));
yearUrbanizedMat(yearUrbanizedMat == 0) = NaN;
yearUrbanized = min(yearUrbanizedMat,[],1, 'omitnan');


statFinalState.householdsIncomeGroup = finalState.householdsCenter;

table_export = table(grid.ID, ...
                     grid.distanceCBD', ...
                     sum(statFinalState.householdsIncomeGroup .* data.householdSizeIncomeGroup')', ...
                     sum(finalState.householdsHousingType(:,:))', ...
                     (finalState.householdsHousingType(1,:) + finalState.householdsHousingType(4,:))', ...
                     finalState.householdsHousingType(2,:)', ...
                     finalState.householdsHousingType(3,:)', ...
                     finalState.rent(1,:)', ...
                     finalState.dwellingSize(1,:)', ...
                     sum(land.coeffLand)', ...
                     yearUrbanized');
table_export.Properties.VariableNames = {'ID_grid' ...
                                         'distance_center'...
                                         'total_people' ...
                                         'total_households' ...
                                         'total_formal' ...
                                         'total_backyard' ...
                                         'total_settlement' ...
                                         'rent_formal' ...
                                         'dwelling_size_formal' ...
                                         'coeff_land_total'...
                                         'year_urbanized'};

writetable(table_export, strcat(filename, scenarioName, '.csv'))

end
