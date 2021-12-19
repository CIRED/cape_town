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
t = 0:1; 

%%% Import inputs %%%
ImportInputsData_LOGIT;

% Transport data
yearTraffic = [t(1):sign(t(length(t))):t(length(t))];
trans = LoadTransportCostCapeTown_LOGIT(option, grid, macro, param, poly, data, yearTraffic, 'SP', 1);
disp ('*** Data and parameters imported succesfully ***')

%% Calibration data

% Data coordinates (SP)
xData = data.spX;
yData = data.spY;

% Data at the SP level
dataPrice = data.spPrice(3,:);
dataDwellingSize = data.spDwellingSize;

% Income classes
% [~, dataIncomeGroup] = max(data.sp2011IncomeDistributionNClass');
dataIncomeGroup = ones(1,length(data.sp2011AverageIncome));
for j = 1:param.numberIncomeGroup-1
    dataIncomeGroup(data.sp2011AverageIncome > data.thresholdIncomeDistribution(j)) = j+1;
end

%%% Import amenities at the SP level %%%
ImportAmenitiesSP
% variablesRegression = {'distance_ocean', 'distance_ocean_2_4', 'slope_1_5', 'slope_5', 'airport_cone2', 'distance_distr_parks', 'distance_biosphere_reserve', 'distance_world_herit', 'distance_train', 'distance_urban_herit'};
variablesRegression = {'distance_ocean', 'distance_ocean_2_4', 'slope_1_5', 'slope_5', 'airport_cone2', 'distance_distr_parks', 'distance_biosphere_reserve', 'distance_train', 'distance_urban_herit'};
% corrplot(tableAmenities(:,variablesRegression))

%%  Estimation of coefficient of construction function

dataNumberFormal = (data.spTotalDwellings - data.spInformalBackyard - data.spInformalSettlement)';
dataDensity = dataNumberFormal./(data.spUnconstrainedArea .* param.maxCoeffLand ./ 1000000);
% selectedDensity = (data.spUnconstrainedArea > 0.2.*1000000.*data.sp2011Area') & dataIncomeGroup > 1 & data.sp2011Distance' < 20 & dataPrice > quantile(dataPrice, 0.1) & ~data.sp2011MitchellsPlain';
% housingSupply = dataDwellingSize(selectedDensity) .* dataDensity(selectedDensity);% ./ (data.spUnconstrainedArea(selectedDensity));
% modelConstruction = fitlm(log(dataPrice(selectedDensity)), log(housingSupply));

% Other possibility
% selectedDensity = (data.spUnconstrainedArea > 0.3.*1000000.*data.sp2011Area') & dataIncomeGroup > 1 & data.sp2011Distance' < 30 & dataPrice > quantile(dataPrice, 0.2) & ~data.sp2011MitchellsPlain' & data.spUnconstrainedArea < quantile(data.spUnconstrainedArea, 0.8);
selectedDensity = (data.spUnconstrainedArea > 0.6.*1000000.*data.sp2011Area') & dataIncomeGroup > 1 & ~data.sp2011MitchellsPlain' & data.sp2011Distance' < 40 & dataPrice > quantile(dataPrice, 0.2) & data.spUnconstrainedArea < quantile(data.spUnconstrainedArea, 0.8);
modelConstruction = fitlm([log(dataPrice(selectedDensity)); log(param.maxCoeffLand.*data.spUnconstrainedArea(selectedDensity)); log(dataDwellingSize(selectedDensity))]', log(dataNumberFormal(selectedDensity))');

% Export outputs
coeff_b = modelConstruction.Coefficients.Estimate(2); % modelConstruction.Coefficients.Estimate(2) ./ (1 + modelConstruction.Coefficients.Estimate(2));
coeff_a = 1 - coeff_b;
% coeff_lambda = (1./coeff_b.^coeff_b) .* exp(coeff_a .* modelConstruction.Coefficients.Estimate(1));  
coeffKappa = (1./coeff_b.^coeff_b) .* exp(modelConstruction.Coefficients.Estimate(1));  

% Correcting data for rents
dataRent = dataPrice.^(coeff_a) .* (param.depreciationRate + InterpolateInterestRateEvolution(macro, t(1))) ./ (coeffKappa .* coeff_b .^ coeff_b);
% dataRent(data.sp2011MitchellsPlain) = min(dataRent(data.sp2011MitchellsPlain), dataPrice(data.sp2011MitchellsPlain) .* (param.depreciationRate + InterpolateInterestRateEvolution(macro, t(1))) ./ data.spFormalDensityHFA(data.sp2011MitchellsPlain));
% dataRent = (param.depreciationRate + InterpolateInterestRateEvolution(macro, t(1))) .* dataPrice ./0.5;

%% With a CES function: fit the function

interestRate = (param.depreciationRate + InterpolateInterestRateEvolution(macro, t(1)));

% Cobb-Douglas: 
simulHousing_CD = coeffKappa.^(1/coeff_a)...
        .*(coeff_b/interestRate).^(coeff_b/coeff_a)...
        .*(dataRent).^(coeff_b/coeff_a);

% CES: 
coeffSigma_CES = 0.75; % 0.7;
coeff_a_CES = 0.1;
coeff_b_CES = 0.9;
coeffKappa_CES = coeffKappa/31; % coeffKappa/17;
simulHousing_CES = coeffKappa_CES .* coeff_a_CES^(-coeffSigma_CES/(1-coeffSigma_CES)) .* (1 - coeff_b_CES^coeffSigma_CES .* (coeffKappa_CES .* dataRent/interestRate).^(coeffSigma_CES - 1)).^(coeffSigma_CES/(1 - coeffSigma_CES));

f1=fit(data.sp2011Distance(selectedDensity), data.spFormalDensityHFA(selectedDensity)','poly5');
f2=fit(data.sp2011Distance(~isnan(simulHousing_CD)), simulHousing_CD(~isnan(simulHousing_CD))','poly5');
f3=fit(data.sp2011Distance(imag(simulHousing_CES) == 0 & ~isnan(simulHousing_CES)), simulHousing_CES(imag(simulHousing_CES) == 0 & ~isnan(simulHousing_CES))','poly5');

close all

nameFolder = strcat('..', slash, 'results_Cape_Town', slash, 'calibration_072019');
PlotHousingFunction

%% Estimation of incomes and commuting parameters

% listLambda = 10.^[0.5:0.05:0.7];
listLambda = 10.^[0.6:0.01:0.65];
listLambda = 10.^[0.6:0.005:0.61];
% listLambda = 3.93:0.023:4;

[modalShares, incomeCenters, timeDistribution, distanceDistribution] = EstimateIncome(param, trans, poly, data, listLambda);

dataModalShares = [7.8 14.8 39.5+0.7 16 8] ./ (7.8+14.8+39.5+0.7+16+8).*100;
dataTimeDistribution = [18.3 32.7 35.0 10.5 3.4] ./ sum(18.3+32.7+35.0+10.5+3.4);
dataDistanceDistribution = [45.6174222, 18.9010734, 14.9972971, 9.6725616, 5.9425438, 2.5368754, 0.9267125, 0.3591011, 1.0464129];

%% Compute accessibility index

% bhattacharyyaModes = -log(sum(sqrt(dataModalShares'./100 .* modalShares)));
bhattacharyyaDistances = -log(sum(sqrt(dataDistanceDistribution'./100 .* distanceDistribution)));
[~, whichLambda] = min(bhattacharyyaDistances);

lambdaKeep = listLambda(whichLambda);
modalSharesKeep = modalShares(:,whichLambda);
timeDistributionKeep = timeDistribution(:,whichLambda);
distanceDistributionKeep = distanceDistribution(:,whichLambda);
incomeCentersKeep = incomeCenters(:,:,whichLambda);

save('./0. Precalculated Inputs/lambda', 'lambdaKeep')
save('./0. Precalculated Inputs/incomeCentersKeep', 'incomeCentersKeep')

[incomeNetOfCommuting, ~, ~, ~] = ComputeIncomeNetOfCommuting(param, trans, grid, poly, data, lambdaKeep, incomeCentersKeep, 'SP', 1);

Simul = distanceDistributionKeep'.*100;
Data = dataDistanceDistribution;
labels = {'0-5','5-10','10-15','15-20','20-25', '25-30', '30-35', '35-40', '> 40'};
b = bar([Data; Simul]', 1);
b(1).FaceColor = [0.3 0.3 0.3];
b(2).FaceColor = [0.3 0.94 0.3];
lgd = legend(b, {'Data', 'Model'});
lgd.FontSize = 16;
lgd.FontName = 'Arial';
set(gca,'xticklabel',labels, 'FontSize', 14, 'FontName', 'Arial');
ylabel('% of commuters', 'FontSize', 16, 'FontName','Arial');
xlabel('Residence-workplace distance (km)', 'FontSize', 16, 'FontName','Arial');
set(gcf,'units','points','position',[50,50,700,320])

%% Estimation of housing demand parameters

% In which areas we actually measure the likelihood
selectedSPForEstimation = (data.spInformalBackyard' + data.spInformalSettlement') ./ data.spTotalDwellings' < 0.1 &... % I remove the areas where there is informal housing, because dwelling size data is not reliable
        dataIncomeGroup > 1; 
    

% Rho is the coefficient for spatial autocorrelation
listRho = 0; 

% Coefficients of the model
listBeta = 0.1:0.2:0.5; %necessarily between 0 and 1
listBasicQ = 5:5:15; %necessarily between 0 and 20

% Utilities
utilityTarget = [300;1000;3000;10000];
listVariation = [0.5:0.3:2];
initUti2 = utilityTarget(2); 
listUti3 = utilityTarget(3) .* listVariation;
listUti4 = utilityTarget(4) .* listVariation;

[parametersScan, scoreScan, parametersAmenitiesScan, modelAmenityScan, parametersHousing, ~] = ...
    EstimateParametersByScanning(incomeNetOfCommuting, dataRent, dataDwellingSize, dataIncomeGroup, ...
    dataDensity, selectedDensity, xData, yData, selectedSPForEstimation, tableAmenities, variablesRegression, ...
    listRho, listBeta, listBasicQ, initUti2, listUti3, listUti4);

%
% Now run the optimization algo with identified value of the parameters
initBeta = parametersScan(1); 
initBasicQ = max(parametersScan(2), 5.1); 

% Utilities
initUti3 = parametersScan(3);
initUti4 = parametersScan(4);

[parameters, scoreTot, parametersAmenities, modelAmenity, parametersHousing, selectedSPRent] = ...
    EstimateParametersByOptimization(incomeNetOfCommuting, dataRent, dataDwellingSize, dataIncomeGroup, ...
    dataDensity, selectedDensity, xData, yData, selectedSPForEstimation, tableAmenities, variablesRegression, ...
    listRho, initBeta, initBasicQ, initUti2, initUti3, initUti4);

%% Generating the map of amenities

% Export the result of the regression of amenities
modelAmenity
% fname='./0. Precalculated inputs/Regression amenities.xls';
% writetable(cell2table([modelAmenity.CoefficientNames'...
%     num2cell(round([modelAmenity.Coefficients.Estimate, modelAmenity.Coefficients.SE, modelAmenity.Coefficients.pValue, ...
%     repmat(modelAmenity.Rsquared.Ordinary, length(modelAmenity.CoefficientNames)), repmat(modelAmenity.Rsquared.Adjusted, length(modelAmenity.CoefficientNames))], 3))]),...
%     fname,'writevariablenames',0);
save('./0. Precalculated inputs/modelAmenity', 'modelAmenity')

% Map of amenties
ImportAmenitiesGrid
amenities = exp(parametersAmenities(2:end)' * table2array(tableAmenitiesGrid(:,variablesRegression))');

%% Exporting and saving

utilitiesCorrected = parameters(3:end) ./ exp(parametersAmenities(1));
calibratedUtility_beta = parameters(1);
calibratedUtility_q0 = parameters(2);

%save('./0. Precalculated inputs/calibratedAmenities', 'amenities')
%save('./0. Precalculated inputs/calibratedUtility_beta', 'calibratedUtility_beta')
save('./0. Precalculated inputs/calibratedUtility_q0', 'calibratedUtility_q0')
save('./0. Precalculated inputs/calibratedUtilities', 'utilitiesCorrected')
save('./0. Precalculated inputs/calibratedHousing_b', 'coeff_b')
save('./0. Precalculated inputs/calibratedHousing_kappa', 'coeffKappa')
save('./0. Precalculated inputs/calibratedHousing_b_CES', 'coeff_b_CES')
save('./0. Precalculated inputs/calibratedHousing_kappa_CES', 'coeffKappa_CES')
save('./0. Precalculated inputs/calibratedHousing_sigma_CES', 'coeffSigma_CES')
disp('*** Parameters saved ***')
