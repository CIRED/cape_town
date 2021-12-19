% Copyright or ?? or Copr. Ecole des Ponts ParisTech 2012
% contributor(s) : Vincent Viguie
% 
% viguie@centre-cired.fr
% 
% This software is a computer program whose purpose is to analyze city growth over time.
% 
% This software is governed by the CeCILL license under French law and
% abiding by the rules of distribution of free software.  You can  use, 
% modify and/ or redistribute the software under the terms of the CeCILL 
% license as circulated by CEA, CNRS and INRIA at the following URL
% "http://www.cecill.info". 
% 
% As a counterpart to the access to the source code and  rights to copy,
% modify and redistribute granted by the license, users are provided only
% with a limited warranty  and the software's author,  the holder of the
% economic rights,  and the successive licensors  have only  limited
% liability. 
% 
% In this respect, the user's attention is drawn to the risks associated
% with loading,  using,  modifying and/or developing or reproducing the
% software by the user in light of its specific status of free software,
% that may mean  that it is complicated to manipulate,  and  that  also
% therefore means  that it is reserved for developers  and  experienced
% professionals having in-depth computer knowledge. Users are therefore
% encouraged to load and test the software's suitability as regards their
% requirements in conditions enabling the security of their systems and/or 
% data to be ensured and,  more generally, to use and operate it in the 
% same conditions as regards security. 
% 
% The fact that you are presently reading this means that you have had
% knowledge of the CeCILL license and that you accept its terms.
% 




function statDynamics = ComputeFinalStatisticsCapeTown(macro,option,trans,land,grid,param,poly,simulation)
% Estimation of different variables for dynamic evolution

%% Useful inputs

interestRateSimulation = zeros(length(simulation.T), 1);
coeffLandSimulation = zeros(length(simulation.T),size(land.coeffLand,1), size(land.coeffLand,2));
for i = 1:length(simulation.T)
    interestRateSimulation(i) = InterpolateInterestRateEvolution(macro, simulation.T(i));
    coeffLandSimulation(i,:,:) = InterpolateLandCoefficientEvolution(land, option, param, simulation.T(i));
end

statDynamics.households(:,:) = sum(simulation.householdsHousingType,2);

radiusCBD = grid.distanceCBD < 6;
rentFormal(:,:) = simulation.rent(:,1,:);
housingFormal(:,:) = simulation.housingSupply(:,1,:);
peopleFormal(:,:) = simulation.householdsHousingType(:,1,:) ;
peopleRDP(:,:) = simulation.householdsHousingType(:,4,:)  ;
peopleBackyard(:,:) = simulation.householdsHousingType(:,2,:)  ;
peopleSettlement(:,:) = simulation.householdsHousingType(:,3,:) ;
dwellingSizeFormal(:,:) = simulation.dwellingSize(:,1,:);


%% Prices

statDynamics.avgRentFormalCBD = nanmean(rentFormal(:,radiusCBD),2);
statDynamics.avgPricePerLandCBD = nanmean(...
             rentFormal(:,radiusCBD),2)./...
             (param.depreciationRate + ppval(macro.splineInterestRate, simulation.T));
     
%% Urban footprint

statDynamics.urbanFootprint = sum(statDynamics.households > param.limitPeopleCityEdge, 2) .* grid.sizeSquare^2; % in sqkm

%% Average prices

whichAveragePrice(:) = (sum(simulation.householdsHousingType(1,:,:) ,2) > 10);
statDynamics.averagePricePerFormalDwelling = nansum(...
    rentFormal(:,whichAveragePrice).*dwellingSizeFormal(:,whichAveragePrice)./repmat(param.depreciationRate + ppval(macro.splineInterestRate, simulation.T),1,size(rentFormal(:,whichAveragePrice),2))./100 ...
    ,2) ./ sum(whichAveragePrice);
statDynamics.weightedAveragePricePerFormalDwelling = nansum(...
    peopleFormal.*rentFormal.*dwellingSizeFormal./repmat(param.depreciationRate + ppval(macro.splineInterestRate, simulation.T),1,size(rentFormal,2))./100 ...
    ,2)./sum(peopleFormal,2);

%% Affordability

statDynamics.averageIncome = ppval(macro.splineIncome, simulation.T) .* ppval(macro.splineInflation, simulation.T)./ppval(macro.splineInflation, 2011 - param.yearBegin);
statDynamics.affordabilityIndex = statDynamics.averagePricePerFormalDwelling ./ statDynamics.averageIncome;
statDynamics.weightedAffordabilityIndex = statDynamics.weightedAveragePricePerFormalDwelling ./ statDynamics.averageIncome;

%% Total by housing types

statDynamics.totalBackyard(:) = sum(simulation.householdsHousingType(:,2,:), 3) ;
statDynamics.totalSettlement(:) = sum(simulation.householdsHousingType(:,3,:), 3) ;
statDynamics.totalPrivate(:) = sum(simulation.householdsHousingType(:,1,:),3) ;
statDynamics.totalRDP(:) = sum(simulation.householdsHousingType(:,4,:) ,3) ;
statDynamics.totalInformal(:) = statDynamics.totalBackyard + statDynamics.totalSettlement;
statDynamics.totalFormalIncludingRDP(:) = statDynamics.totalPrivate + statDynamics.totalRDP;
statDynamics.totalPopulation = statDynamics.totalFormalIncludingRDP + statDynamics.totalInformal;

% poor in formal 
for i = 1:length(simulation.T)
    statDynamics.numberPoorInFormal(i) = sum(sum(simulation.households(i, 1, poly.incomeGroup(i,:) < 3, :)));
end


%% Average distance to the CBD

statDynamics.averageDistanceToCBD = sum(statDynamics.households(:,:).*repmat(grid.distanceCBD, length(simulation.T),1),2) ./ sum(statDynamics.households(:,:),2);

householdsCenterPoor(:,:) = sum(simulation.householdsCenter(:,1:2,:), 2);
householdsCenterRich(:,:) = sum(simulation.householdsCenter(:,3:4,:), 2);
statDynamics.averageDistanceToCBDPoor = sum(householdsCenterPoor .*repmat(grid.distanceCBD, length(simulation.T),1), 2) ./ sum(householdsCenterPoor,2);
statDynamics.averageDistanceToCBDRich = sum(householdsCenterRich .*repmat(grid.distanceCBD, length(simulation.T),1), 2) ./ sum(householdsCenterRich,2);

statDynamics.averageDistanceToCBD_RDP = sum(peopleRDP .*repmat(grid.distanceCBD, length(simulation.T),1), 2) ./ sum(peopleRDP,2);
statDynamics.averageDistanceToCBDBackyard = sum(peopleBackyard .*repmat(grid.distanceCBD, length(simulation.T),1), 2) ./ sum(peopleBackyard,2);


end