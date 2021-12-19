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
disp ('**************** NEDUM-Cape-Town - Estimation of parameters ****************');

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
t = 0:2; 

%%% Import inputs %%%
ImportInputsData;

%%% Define minimum lot-size %%% 
param.miniLotSize = min(data.spDwellingSize((data.spInformalSettlement + data.spInformalBackyard)./data.spTotalDwellings < 0.2));

%% Calibration data

% Data coordinates (SP)
xData = data.spX;
yData = data.spY;

% Data at the SP level
dataPrice = data.spPrice(2,:);
dataDwellingSize = data.spDwellingSize;

% Income classes
% dataIncomeGroup = ones(1,length(dataDwellingSize));
% for i = 1:length(dataIncomeGroup)
%     for j = 2:param.numberIncomeGroup
%         if data.sp2011AverageIncome(i) > data.thresholdIncomeDistribution(j-1) 
%             dataIncomeGroup(i) = j;
%         end
%     end
% end
[~, dataIncomeGroup] = max(data.sp2011IncomeDistributionNClass'); % Other possibility


% Income net of commuting costs, extrapolated at the SP level
% income = unique(poly.averageIncome(1,:))';
% transportCostIndex = zeros(4,length(dataPrice));
% for i = 1:4
%     transportCostIndex(i,:) = griddata(grid.xCoord, grid.yCoord, trans.generalizedCostIndex(i,:), xData, yData);
% end

%%% Import amenities at the SP level %%%
ImportAmenitiesSP
variablesRegression = {'distance_ocean', 'distance_ocean_2_4', 'slope_1_5', 'slope_5', 'airport_cone2', 'distance_distr_parks', 'distance_urban_herit', 'distance_protected_envir'};
% corrplot(table_amenities(:,which_regression))

%%  Estimation of coefficient of construction function

dataNumberFormal = (data.spTotalDwellings - data.spInformalBackyard - data.spInformalSettlement)';
dataDensity = dataNumberFormal./(data.spUnconstrainedArea .* param.maxCoeffLand ./ 1000000);
selectedDensity = (data.spUnconstrainedArea > 0.25.*1000000.*data.sp2011Area') & dataIncomeGroup > 1 & data.sp2011Distance' < 20 & dataPrice > 750 & ~data.sp2011MitchellsPlain';
housingSupply = dataDwellingSize(selectedDensity) .* dataDensity(selectedDensity) ./ (data.spUnconstrainedArea(selectedDensity));
modelConstruction = fitlm(log(dataPrice(selectedDensity)), log(housingSupply));
% modelConstruction = LinearModel.fit([log(dataPrice(selectedDensity)); log(dataDwellingSize(selectedDensity)); log(data.spUnconstrainedArea(selectedDensity) .* param.maxCoeffLand)]', log(dataNumberFormal(selectedDensity))');

coeff_b = modelConstruction.Coefficients.Estimate(2); % modelConstruction.Coefficients.Estimate(2) ./ (1 + modelConstruction.Coefficients.Estimate(2));
coeff_a = 1 - coeff_b;
% coeff_lambda = (1./coeff_b.^coeff_b) .* exp(coeff_a .* modelConstruction.Coefficients.Estimate(1));  
coeff_lambda = (1./coeff_b.^coeff_b) .* exp(modelConstruction.Coefficients.Estimate(1));  


% Correcting data for rents
dataRent = dataPrice.^(coeff_a) .* (param.depreciationRate + InterpolateInterestRateEvolution(macro, t(1))) ./ (coeff_lambda .* coeff_b .^ coeff_b);
dataRent(data.sp2011MitchellsPlain) = min(dataRent(data.sp2011MitchellsPlain), dataPrice(data.sp2011MitchellsPlain) .* (param.depreciationRate + InterpolateInterestRateEvolution(macro, t(1))) ./ data.spFormalDensityHFA(data.sp2011MitchellsPlain));
% dataRent = (param.depreciationRate + InterpolateInterestRateEvolution(macro, t(1))) .* dataPrice ./0.5;


%% Estimation of housing demand parameters

% In which areas we actually measure the likelihood
selectedSPForEstimation = (data.spInformalBackyard' + data.spInformalSettlement') ./ data.spTotalDwellings' < 0.1 &... % I remove the areas where there is informal housing, because dwelling size data is not reliable
        dataIncomeGroup > 1; 
    
% Income net of commuting costs
incomeNetOfCommuting = unique(poly.averageIncome(2,:))' - transportCostIndex;


% Rho is the coefficient for spatial autocorrelation
listRho = 0; 

% Coefficients of the model
listBeta = 0.2:0.2:0.8; %necessarily between 0 and 1
listBasicQ = 9:3:15; %necessarily between 0 and 20

% Utilities
utilityTarget = [300;1000;5000;10000];
listVariation = [0.2:0.25:2.2];
initUti2 = utilityTarget(2); 
listUti3 = utilityTarget(3) .* listVariation;
listUti4 = utilityTarget(4) .* listVariation;

[parametersScan, scoreScan, parametersAmenitiesScan, modelAmenityScan, parametersHousing, ~] = ...
    EstimateParametersByScanning(incomeNetOfCommuting, dataRent, dataDwellingSize, dataIncomeGroup, ...
    dataDensity, selectedDensity, xData, yData, selectedSPForEstimation, tableAmenities, variablesRegression, ...
    listRho, listBeta, listBasicQ, initUti2, listUti3, listUti4);

%% 

% Now run the optimization algo with identified value of the parameters
initBeta = parametersScan(1); %0.38; %0.1:0.025:0.45; %necessarily between 0 and 1
initBasicQ = min(parametersScan(2), 10.1); %12;

% Utilities
initUti3 = parametersScan(3);
initUti4 = parametersScan(4);

[parameters, scoreTot, parametersAmenities, modelAmenity, parametersHousing, selectedSPRent] = ...
    EstimateParametersByOptimization(incomeNetOfCommuting, dataRent, dataDwellingSize, dataIncomeGroup, ...
    dataDensity, selectedDensity, xData, yData, selectedSPForEstimation, tableAmenities, variablesRegression, ...
    listRho, initBeta, initBasicQ, initUti2, initUti3, initUti4);

%% Generating the map of amenities

modelAmenity

% Map of amenties
ImportAmenitiesGrid
amenities = exp(parametersAmenities(2:end)' * table2array(tableAmenitiesGrid(:,variablesRegression))');

%% Exporting and saving

utilitiesCorrected = parameters(3:end) ./ exp(parametersAmenities(1));
calibratedUtility_beta = parameters(1);
calibratedUtility_q0 = parameters(2);

save('./0. Precalculated inputs/calibratedAmenities', 'amenities')
save('./0. Precalculated inputs/calibratedUtility_beta', 'calibratedUtility_beta')
save('./0. Precalculated inputs/calibratedUtility_q0', 'calibratedUtility_q0')
save('./0. Precalculated inputs/calibratedUtilities', 'utilitiesCorrected')
save('./0. Precalculated inputs/calibratedHousing_b', 'coeff_b')
save('./0. Precalculated inputs/calibratedHousing_lambda', 'coeff_lambda')
disp('*** Parameters saved ***')
