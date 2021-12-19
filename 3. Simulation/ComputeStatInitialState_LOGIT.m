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




function output = ComputeStatInitialState_LOGIT(trans,land,grid,macro,param,option,initialState,poly, yearEquilibrium)
% calcule divers chifres sur un etat d'equilibre polycentrique

carto2 = @(x) scatter(grid.xCoord,grid.yCoord,150,x,'.');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mean values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

statInitial.households = initialState.householdsHousingType ; %.* land.coeffLand;
statInitial.population = sum(sum(statInitial.households));
population = statInitial.population;

% Function to estimate means
ComputeMean = @(input) ComputeMeanFunction(grid, param, ones(size(sum(initialState.householdsHousingType,1))), sum(initialState.householdsHousingType,1), sum(land.coeffLand,1), population, input);

% Monetary cost of transport
monetaryCost = ppval(macro.splineFuelCost,yearEquilibrium) + param.Tax;

%coeffLand
statInitial.coeffLand = land.coeffLand;

% Average income
income = InterpolateIncomeEvolution(macro,param,option,grid, poly,yearEquilibrium);
statInitial.householdsIncomeGroup = initialState.householdsCenter;
[~, statInitial.incomeGroup] = max(initialState.householdsCenter);
statInitial.incomeGroup(sum(initialState.householdsCenter) == 0) = NaN;

[~, whichYear] = min(abs((trans.yearTransport - param.yearBegin) - yearEquilibrium));
statInitial.averageIncome = sum(trans.averageIncome(:,:,whichYear) .* initialState.householdsCenter, 1)./sum(initialState.householdsCenter,1);

% Modal distribution
for i = 1:trans.numberModes
    modalSharePerX = permute(nansum(permute(trans.modalShares(:, :, i, :, whichYear), [1 2 4 3 5]) .* trans.ODflows(:,:,:,whichYear), 1) ./ nansum(trans.ODflows(:,:,:,whichYear), 1), [3 2 1]);
    statInitial.fractionMode(i) = sum(nansum(modalSharePerX .* initialState.householdsCenter)) ./ sum(nansum(initialState.householdsCenter));
end

% Mean distance by mode for a direct trip 
% statInitial.distance = ComputeMean(ComputeMeanPolycentric(initialState.householdsCenter, sum(trans.whichModes(:,:,:,1).*repmat(trans.distanceOutput(1,:,:,:),size(trans.whichModes,1),1),3)));

% Mean residence-workplace distance
% statInitial.commutingDistance = sum(initialState.householdsCenter .* trans.distanceOutput(:,:,1), 1) ./ sum(initialState.householdsCenter, 1);

% Total urbanized area
statInitial.surfaceUrbanized = ComputeTotalSumCity(grid,param,sum((initialState.housingSupply > param.minDensityUrban).*land.coeffLand,1)./land.maxCoeffLand);

% Total floor space built (all type of housing)
statInitial.totalFloorSpace = ComputeTotalSumCity(grid,param,sum(initialState.housingSupply.*land.coeffLand,1)./land.maxCoeffLand);
% Total floor space built for formal housing
statInitial.totalFloorSpaceFormalPrivate = ComputeTotalSumCity(grid,param, initialState.housingSupply(1,:).*(initialState.housingSupply(1,:)>param.minDensityUrban).*land.coeffLand(1,:)./land.maxCoeffLand);

output = statInitial;

end

function output = ComputeMeanPolycentric(people,input)
% Compute mean for a variable, weighted by number of people

output = nansum(input.*people,1)./nansum(people,1);

end

function output = ComputeMeanFunction(grid,param,limit,people,coeffLand_ici,population,input)
% Compute mean for a variable, weighted by number of people, by households

filter = (~isnan(people))&(~isnan(input));
if length(grid.sizeSquare)==1
   output = sum(input(filter).*limit(filter).*people(filter).*coeffLand_ici(filter).*grid.sizeSquare.*grid.sizeSquare,2)./population;
else
   output = sum(input(filter).*limit(filter).*people(filter).*coeffLand_ici(filter).*grid.sizeSquare(filter).*grid.sizeSquare(filter),2)./population;
end

end

function output = ComputeTotalSumCity(grid,param,entree)
% Calculate the sum for the whole city
% careful, does not integrate coeffLand

output = sum(entree.*grid.sizeSquare.*grid.sizeSquare,2);

end
