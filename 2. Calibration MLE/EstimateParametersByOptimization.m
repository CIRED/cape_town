function [parameters, scoreTot, parametersAmenities, modelAmenity, parametersHousing, selectedRents] = ...
    EstimateParametersByOptimization(incomeNetOfCommuting, dataRent, dataDwellingSize, dataIncomeGroup, ...
    dataHouseholdDensity, selectedDensity, xData, yData, selectedSP, tableAmenities, variablesRegression, ...
    initRho, initBeta, initBasicQ, initUti2, initUti3, initUti4)

% Automated estimation of the parameters of NEDUM by maximizing log likelihood

% Here we minimize the log-likelihood using fminsearch 

%% Neighbor matrix

% [ weighting ] = create_weighting(x_data, y_data, 4);


%% Data as matrices, where should we regress (remove where we have no data)

% Where is which class
incomeNetOfCommuting = incomeNetOfCommuting(2:4,:); % We remove income group 1
groupLivingSpMatrix = (incomeNetOfCommuting > 0);
for i = 1:3
    groupLivingSpMatrix(i, dataIncomeGroup ~= i + 1) = false;
end
selectedTransportMatrix = (sum(groupLivingSpMatrix) == 1);
incomeNetOfCommuting(incomeNetOfCommuting < 0) = NaN;

selectedRents = ~isnan(dataRent) & selectedTransportMatrix & selectedSP;
selectedDwellingSize = ~isnan(dataDwellingSize) & ~isnan(dataRent) & selectedTransportMatrix & selectedSP;
selectedDensity = selectedDwellingSize & selectedDensity;

% For the regression of amenities
tableRegression = tableAmenities(selectedRents,:);
predictorsAmenitiesMatrix = table2array(tableRegression(:,variablesRegression));
predictorsAmenitiesMatrix = [ones(size(predictorsAmenitiesMatrix,1),1), predictorsAmenitiesMatrix];
modelAmenity = 0;

% Sort data
% dataDwellingSizeSelected = dataDwellingSize(selectedDwellingSize);
% dataRentSelected = dataRent(selectedRents);
% dataRentSelectedMatrix = repmat(dataRentSelected, size(groupLivingRentsMatrix,1),1);


%% Useful functions (precalculations for rents and dwelling sizes, likelihood function) 


% Function for dwelling sizes
% We estimate calcule_hous directly from data from rents (no extrapolation)
CalculateDwellingSize = @(beta, basic_q, incomeTemp, rentTemp) beta .* incomeTemp ./ rentTemp + (1-beta) .* basic_q;

% Log likelihood for a lognormal law
ComputeLogLikelihood = @(sigma, error)...
         nansum(- log(2*pi*sigma^2)/2 -  1./(2*sigma.^2).*(error).^2);
    


%% Optimization algorithm

% Initial value of parameters
initialVector = [initBeta; initBasicQ; initUti3; initUti4]; % So far, no spatial autocorrelation

% Function that will be minimized 
optionRegression = 0;
minusLogLikelihoodModel = @(X0) - LogLikelihoodModel(X0, initUti2, incomeNetOfCommuting, groupLivingSpMatrix, dataDwellingSize, selectedDwellingSize, xData, yData, dataRent, selectedRents, dataHouseholdDensity, selectedDensity, predictorsAmenitiesMatrix, tableRegression, variablesRegression, CalculateDwellingSize, ComputeLogLikelihood, optionRegression);

% Optimization w/ lower and upper bounds
lowerBounds = [0.1; 3; 0; 0];
upperBounds = [1; 18; 10^6; 10^7];
optionsOptim = optimset('Display', 'iter');
[parameters, scoreTot, exitFlag] = fmincon(minusLogLikelihoodModel, initialVector, [], [], [], [], lowerBounds, upperBounds, [], optionsOptim);

% Estimate the function to get the parameters for amenities
optionRegression = 1;
[~, ~, ~, ~, ~, parametersAmenities, modelAmenity, parametersHousing] = LogLikelihoodModel(parameters, initUti2, incomeNetOfCommuting, groupLivingSpMatrix, dataDwellingSize, selectedDwellingSize, xData, yData, dataRent, selectedRents, dataHouseholdDensity, selectedDensity, predictorsAmenitiesMatrix, tableRegression, variablesRegression, CalculateDwellingSize, ComputeLogLikelihood, optionRegression);

disp('*** Estimation of beta and q0 done ***')

end 

function [interval] = confidence_interval(indices_max,quoi_indices,compute_score)
%% Confidence interval

d_beta=1.05;

fprintf('\n');
beta_interval = zeros(size(indices_max));
for index = 1:size(indices_max,1),
    indices_ici = indices_max;
    score_tmp = compute_score(indices_ici);
    
    indices_ici = indices_max;
    score_tmp2 = compute_score(indices_ici);
    
    indices_ici = indices_max;
    indices_ici(index) = indices_ici(index) - indices_ici(index)*(d_beta-1);
    score_tmp3 = compute_score(indices_ici);
    
    indices_ici = indices_max;
    dd_l_beta = -(score_tmp2 + score_tmp3 - 2*score_tmp) / (indices_ici(index)*(d_beta-1))^2;
    beta_interval(index) = 1.96 / (sqrt( abs(dd_l_beta)));
    fprintf('%s\t\t%g (%g ; %g)\n',quoi_indices{index},indices_ici(index),indices_ici(index)-beta_interval(index),indices_ici(index)+beta_interval(index))
end
interval = [indices_ici-beta_interval,indices_ici+beta_interval];

end



function [ weighting ] = create_weighting( X_data,Y_data,pas )
%create_weighting creates weights matrix for spatial autocorrelation
%analysis
%   Detailed explanation goes here

if size(X_data)~=size(Y_data),
    return;
end

index_data = true(size(X_data));

%put data as a row vector
if size(X_data,1)>size(X_data,2),
    X_data=X_data';
    Y_data=Y_data';
end

X_data2=repmat(X_data,size(X_data,2),1);
Y_data2=repmat(Y_data,size(Y_data,2),1);

X_data3=repmat(X_data',1,size(X_data,2));
Y_data3=repmat(Y_data',1,size(X_data,2));

dist=(X_data2-X_data3).^2+(Y_data2-Y_data3).^2;
weighting=(dist<pas^2);
weighting=logical(weighting-eye(size(weighting)));
end


