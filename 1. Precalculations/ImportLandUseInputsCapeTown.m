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




function output = ImportLandUseInputsCapeTown(grid, option, param, macro, data)
% import land use datasets
% To be used with the grid_500 from June 2018 on, with constraints from EA data

disp ('importing land use...');


global path_nedum

land.maxCoeffLand = param.maxCoeffLand;
land.maxCoeffLandBackyard = param.maxCoeffLandBackyard;
land.maxCoeffLandSettlement = param.maxCoeffLandSettlement;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import of land_use
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

areaPixel = grid.sizeSquare^2 .* 1000000;

% Land Cover Data from our estimation (see R code for details)
importfile([path_nedum,'grid_NEDUM_Cape_Town_500.csv'])

land.urbanized = urban'/areaPixel;
land.informal = informal'/areaPixel; 
land.coeffLandNoUrbanEdge = (unconstrained_out' + unconstrained_UE')/areaPixel;
land.coeffLandUrbanEdge = unconstrained_UE'/areaPixel;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Import of RDP houses data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we estimate the number of RDP/BNG dwellings and the area available
% for backyarding in each subplace using GV2012 data
land.numberPropertiesRDP = data.gridCountRDPfromGV;
numberPropertiesRDP2000 = data.gridCountRDPfromGV .* (1 - grid.distanceCBD ./ max(grid.distanceCBD(data.gridCountRDPfromGV > 0)));
land.areaRDP = data.gridAreaRDPfromGV .* param.sizeRDP./(param.sizeBackyard + param.sizeRDP)./areaPixel;
land.areaBackyard = data.gridAreaRDPfromGV .* param.sizeBackyard./(param.sizeBackyard + param.sizeRDP)./areaPixel;

land.coeffLandBackyard = min(land.urbanized, land.areaBackyard);

method = 'linear';

if (option.futureConstructionRDP == 1)
    
% Scenario with future construction
% param.backyardSizeFuture will indicate whether backyarding is possible
% in future RDP/BNG settlements

    % Evolution for the Backyard area
    importfile([path_nedum, 'grid_new_RDP_projects.csv'])
    
    % Given the pace, at what years will ST and LT projects by completed
    yearBeginRDP = 2015;
    yearRDP = (yearBeginRDP:2100) - param.yearBegin;
    numberRDP = ppval(macro.splineRDP, yearRDP);
    [~, yearShortTerm] = min(abs(sum(total_yield_DU_ST) - (numberRDP - numberRDP(1))));
    [~, yearLongTerm] = min(abs(sum(total_yield_DU_LT + total_yield_DU_ST) - (numberRDP - numberRDP(1))));

    areaRDPShortTerm = min(area_ST, (param.sizeBackyard + param.sizeRDP).* total_yield_DU_ST);
    areaRDPLongTerm = min(min(area_ST + area_LT, (param.sizeBackyard + param.sizeRDP).* (total_yield_DU_ST + total_yield_DU_LT)),areaPixel);

    % Share of pixel for RDP houses and backyards in ST and LT
    areaBackyardShortTerm = land.areaBackyard + max(areaRDPShortTerm' - total_yield_DU_ST' .* param.sizeRDP, 0)./areaPixel;
    areaRDPShortTerm = land.areaRDP + min(total_yield_DU_ST' .* param.sizeRDP, area_ST')./areaPixel;
    areaBackyardShortTerm = min(areaBackyardShortTerm, param.maxCoeffLand - areaRDPShortTerm);
    areaBackyardLongTerm = land.areaBackyard + max(areaRDPLongTerm' - (total_yield_DU_LT' + total_yield_DU_ST') .* param.sizeRDP, 0)./areaPixel;
    areaRDPLongTerm = land.areaRDP + min((total_yield_DU_LT' + total_yield_DU_ST') .* param.sizeRDP, areaRDPLongTerm')./areaPixel;
    areaBackyardLongTerm = min(areaBackyardLongTerm, param.maxCoeffLand - areaRDPLongTerm);

    yearDataInformal = [2000 - param.yearBegin; yearBeginRDP - param.yearBegin; yearShortTerm; yearLongTerm]';
    land.splineLandBackyard = interp1(yearDataInformal,  [land.areaBackyard; land.areaBackyard; areaBackyardShortTerm; areaBackyardLongTerm], method, 'pp');
    land.splineLandRDP = interp1(yearDataInformal,  [land.areaRDP; land.areaRDP; areaRDPShortTerm; areaRDPLongTerm], method, 'pp');
    land.splineEstimateRDP = interp1(yearDataInformal, [numberPropertiesRDP2000; land.numberPropertiesRDP; land.numberPropertiesRDP + total_yield_DU_ST'; land.numberPropertiesRDP + total_yield_DU_ST' + total_yield_DU_LT'], method, 'pp');

elseif option.futureConstructionRDP == 0
% Scenario with no future construction of RDP

    yearDataInformal = [1990; 2040]' - param.year_begin;
    land.spline_land_backyard = interp1(yearDataInformal,  [land.areaBackyard; land.areaBackyard], method, 'pp');
    land.spline_land_RDP = interp1(yearDataInformal,  [land.areaRDP; land.areaRDP], method, 'pp');
    land.spline_estimate_RDP = interp1(yearDataInformal, [land.numberPropertiesRDP; land.numberPropertiesRDP], method, 'pp');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import of informal settlement land use from Census (to be removed)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% land.informal = data_courbe.informal_settlement_grid .* param.size_shack./1000000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeff_land for each housing type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Privately developed housing
land.coeffLandPrivateUrbanEdge = (land.coeffLandUrbanEdge - land.informal - min(land.areaRDP + land.areaBackyard, land.urbanized)).*land.maxCoeffLand;
land.coeffLandPrivateNoUrbanEdge = (land.coeffLandNoUrbanEdge - land.informal - min(land.areaRDP + land.areaBackyard, land.urbanized)).*land.maxCoeffLand;
land.coeffLandPrivateUrbanEdge(land.coeffLandPrivateUrbanEdge<0) = 0;
land.coeffLandPrivateNoUrbanEdge(land.coeffLandPrivateNoUrbanEdge<0) = 0;

% Evolution of constraints
if option.urbanEdge == 0
    yearConstraints = [1990; param.yearUrbanEdge - 1; param.yearUrbanEdge; 2040]' - param.yearBegin;
    land.splineLandConstraints = interp1(yearConstraints, [land.coeffLandUrbanEdge; land.coeffLandUrbanEdge; land.coeffLandNoUrbanEdge; land.coeffLandNoUrbanEdge], method, 'pp');
else
    yearConstraints = [1990;2040]' - param.yearBegin;
    land.splineLandConstraints = interp1(yearConstraints, [land.coeffLandUrbanEdge; land.coeffLandUrbanEdge], method, 'pp');
end

% For the initial state
land.coeffLandPrivate = land.coeffLandPrivateUrbanEdge;

% For the initial state
% Backyard and settlements
land.coeffLandBackyard = land.coeffLandBackyard .* land.maxCoeffLandBackyard;
land.coeffLandBackyard(land.coeffLandBackyard < 0) = 0;
land.coeffLandSettlement = land.informal .* land.maxCoeffLandSettlement;
land.coeffLandRDP = ones(size(land.coeffLandPrivate));


% For the initial state
land.coeffLand = [land.coeffLandPrivate; land.coeffLandBackyard; land.coeffLandSettlement; land.coeffLandRDP];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building limit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

centerRegulation = (grid.distanceCBD <= param.historicRadius);
outsideRegulation = (grid.distanceCBD > param.historicRadius);

%The variable that will be used in practice:
land.housingLimit = param.limitHeightCenter*1000000*centerRegulation + param.limitHeightOut*1000000*outsideRegulation;
land.housingLimitPolicy = land.housingLimit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output = land;


end

