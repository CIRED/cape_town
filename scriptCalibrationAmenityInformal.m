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

clear;
clc;
close all;
tic();
disp ('**************** NEDUM-Cape-Town - Polycentric Version - Formal and Informal housing ****************');

global slash
slash = '/'; % To be changed when on Mac or Windows
% Load function addresses
pathProgram;
% Load data address
pathData;

%% Parameters and options

%%% Import Parameters %%%
param = ChoiceParameters;

%%% Set options for the simulation %%% 
option.polycentric = 1;
option.futureConstructionRDP = 1;
option.logit = 1;
option.taxOutUrbanEdge = 0;
option.doubleStoreyShacks = 0;
option.urbanEdge = 0; %1 means we keep the urban edge
option.ownInitializationSolver = 0;
option.adjustHousingSupply = 1;
option.loadTransportTime = 1;
option.sortEmploymentCenters = 1;

%%% Parameters for the scenario %%%
param.yearUrbanEdge = 2016; % in case option.urban_edge = 0, the year the constraint is removed
param.taxOutUrbanEdge = 10000;
param.coeffDoubleStorey = 0.02; % Random for now

%%% Years for the simulation %%%
param.yearBegin = 2011;
t = 0:1:6; % to go up to 2017

%%% Import inputs for the simulation %%%
ImportInputsData;

%%%  Calibrated parameters %%% 
load('./0. Precalculated inputs/calibratedAmenities.mat')
land.amenities = amenities ./ mean(amenities);
load('./0. Precalculated inputs/calibratedUtility_beta.mat')
load('./0. Precalculated inputs/calibratedUtility_q0.mat')
param.beta = calibratedUtility_beta;
param.alpha = 1 - param.beta; 
param.basic_q = calibratedUtility_q0;
load('./0. Precalculated inputs/calibratedHousing_b.mat')
load('./0. Precalculated inputs/calibratedHousing_lambda.mat')
param.coeff_b = coeff_b;
param.coeff_a = 1 - param.coeff_b;
param.coeff_grandA = coeff_lambda;
load('./0. Precalculated inputs/calibratedParamAmenities.mat')
param.amenityBackyard = calibratedParamAmenities(1);
param.amenitySettlement = calibratedParamAmenities(2);

%%% Define minimum lot-size %%% 
param.miniLotSize = min(data.spDwellingSize((data.spInformalSettlement + data.spInformalBackyard)./data.spTotalDwellings < 0.1));
% We transform the rent per unit of land to a rent per unit of floor area:
param.agriculturalRent2011 = param.agriculturalRent2011.^(param.coeff_a) .* (param.depreciationRate + InterpolateInterestRateEvolution(macro, 0)) ./ (param.coeff_grandA .* param.coeff_b .^ param.coeff_b);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% List of parameters
listAmenityBackyard = 0.6:0.032:0.9; % listAmenityBackyard = 0.5:0.033:0.8;
listAmenitySettlement = 0.6:0.034:0.9; % listAmenitySettlement = 0.5:0.033:0.8;
listParam = combvec(listAmenityBackyard, listAmenitySettlement);
housingTypeTotal = zeros(4, size(listParam, 2));

Uo_init = poly.averageIncome(1,:) .* 0.2;
option.ownInitializationSolver = 1;

param.maxIteration = 2 * param.maxIteration;

sumHousingTypes = @(initialState) sum(initialState.householdsHousingType, 2);

for i = 1:size(listParam, 2)

    param.amenityBackyard = listParam(1,i);
    param.amenitySettlement = listParam(2,i);
    initialState = RunEquilibriumSolverNEDUM(t(1),trans,option,land,grid,macro,param,poly,Uo_init,param.housing_in);
    statInitialState = ComputeStatInitialState(trans,land,grid,macro,param,option,initialState,poly,t(1));
    housingTypeTotal(:,i) = sumHousingTypes(initialState);
    Uo_init = initialState.utility;
    
end

save 'calibration housing types'

disp('*** End of simulations for chosen parameters ***');

%% Test

housingTypeData = [sum(data.gridFormal(data.limitCapeTown)) - sum(initialState.householdsHousingType(4,:)); sum(data.gridInformalBackyard(data.limitCapeTown)); ...
    sum(data.gridInformalSettlement(data.limitCapeTown)); sum(data.gridFormal(data.limitCapeTown) + data.gridInformalBackyard(data.limitCapeTown) + data.gridInformalSettlement(data.limitCapeTown))];

distanceShare = abs(housingTypeTotal(1:3,:)./sum(housingTypeTotal) - repmat(housingTypeData(1:3)./housingTypeData(4), 1, size(housingTypeTotal,2)));

% scatter(list_param(2,:), distance_share(3,:), 200, list_param(1,:), '.')

distanceShareScore = distanceShare(2,:) + distanceShare(3,:);
[~, which] = min(distanceShareScore);
calibratedParamAmenities = listParam(:,which);

save('./0. Precalculated inputs/calibratedParamAmenities.mat', 'calibratedParamAmenities')
disp('*** Amenity parameters saved ***');
