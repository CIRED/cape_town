function [parameters, scoreTot, parametersAmenities, modelAmenity, parametersHousing, selectedRents] = ...
    EstimateParametersByScanning(incomeNetOfCommuting, dataRent, dataDwellingSize, dataIncomeGroup, ...
    dataHouseholdDensity, selectedDensity, xData, yData, selectedSP, tableAmenities, variablesRegression, ...
    initRho, listBeta, listBasicQ, initUti2, listUti3, listUti4)

% Automated estimation of the parameters of NEDUM by maximizing log likelihood

% Here we scan a set of values for each parameters and determine the value
% of the log-likelihood (to see how the model behaves). 
% In EstimateParameters By Optimization, we use the minimization algorithm from Matlab to converge towards solution 

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

% Function that will be minimized 
optionRegression = 0;

% Initial value of parameters
combinationInputs = combvec(listBeta, listBasicQ, listUti3, listUti4); % So far, no spatial autocorrelation

% Scanning of the list
scoreAmenities = - 10000 .* ones(1, size(combinationInputs,2));
scoreDwellingSize = - 10000 .* ones(1, size(combinationInputs,2));
scoreIncomeSorting = - 10000 .* ones(1, size(combinationInputs,2));
scoreHousing = - 10000 .* ones(1, size(combinationInputs,2));
iterPrint = floor(size(combinationInputs,2) ./ 20);
fprintf('\nDone: ')
for index = 1:size(combinationInputs,2)
    [~, scoreAmenities(index), scoreDwellingSize(index), scoreIncomeSorting(index), scoreHousing(index)] = LogLikelihoodModel(combinationInputs(:,index), initUti2, incomeNetOfCommuting, groupLivingSpMatrix, dataDwellingSize, selectedDwellingSize, xData, yData, dataRent, selectedRents, dataHouseholdDensity, selectedDensity, predictorsAmenitiesMatrix, tableRegression, variablesRegression, CalculateDwellingSize, ComputeLogLikelihood, optionRegression);
    if floor(index / iterPrint) == index/iterPrint
        fprintf('%0.f%%  ', round(index / size(combinationInputs,2) .* 100));
    end
end
fprintf('\nScanning complete');
fprintf('\n');

scoreVect = scoreAmenities + scoreDwellingSize + scoreIncomeSorting + scoreHousing;
[scoreTot, which] = max(scoreVect); 
parameters = combinationInputs(:, which);

% Estimate the function to get the parameters for amenities
optionRegression = 1;
[~, ~, ~, ~, ~, parametersAmenities, modelAmenity, parametersHousing] = LogLikelihoodModel(parameters, initUti2, incomeNetOfCommuting, groupLivingSpMatrix, dataDwellingSize, selectedDwellingSize, xData, yData, dataRent, selectedRents, dataHouseholdDensity, selectedDensity, predictorsAmenitiesMatrix, tableRegression, variablesRegression, CalculateDwellingSize, ComputeLogLikelihood, optionRegression);
    
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


function [ utility ] = utilityFromRents(Ro, income, basic_q, beta)
% Precalculation of the utility
%   [ utili ] = utilityFromRents(rent,income,basic_q,beta) estimates
%   utility for a rent "Ro" and an income (net of transport costs)
%   "revenu". 
% Without amenities 

utility = (1-beta)^(1-beta)*beta^beta...
          .*sign(income - basic_q.*Ro)...
          .*abs(income - basic_q.*Ro)...
          ./(Ro.^beta);
utility((income - basic_q.*Ro) < 0) = 0;

% If the income is 0 then utility is 0 (we need to specify it to avoid
% NaNs later)
utility(income==0) = 0;


end

% Precalculations for rents, as a function
% The output of the function is a griddedInterpolant object, that gives the
% log of rents as a function of the log utility and the log income
function griddedRents = InterpolateRents(beta, basic_q, netIncome)
    
        % Decomposition for the interpolation (the more points, the slower the code)
        decompositionRent = [10.^([-9,-4,-3,-2]), 0.02:0.02:0.08, 0.1:0.05:1.4, 1.5:0.1:2.5, 100, 10^9]; 
        decompositionIncome = [10.^([-9,-4,-3,-2]), 0.02:0.02:0.08, 0.1:0.05:1.4, 1.5:0.1:2.5, 100, 10^9];

        % Min and Max values for the decomposition
        choiceIncome = 100000 .* decompositionIncome ;
        incomeMatrix = repmat(choiceIncome,length(choiceIncome),1)';
        if basic_q > 0.5
            choiceRent = choiceIncome./basic_q; % the maximum rent is the rent for which u = 0
        else
            choiceRent = 1000 .* decompositionRent;
        end
        rentMatrix = choiceRent' * decompositionRent;

        utilityMatrix = utilityFromRents(incomeMatrix, rentMatrix, basic_q, beta);
        solusRentTemp = @(x,y) griddata(incomeMatrix,utilityMatrix,rentMatrix,x,y);
        
        % Redefine a grid (to use griddedInterpolant)
        utilityVectLog = -1:0.1:log(max(max(10.*netIncome)));
        incomeLog = (-1:0.2:log(max(max(10.*netIncome))))';
        rentLog = log(solusRentTemp(exp(incomeLog), exp(utilityVectLog)));
        griddedRents = griddedInterpolant({utilityVectLog, incomeLog}, rentLog, 'linear', 'none');
        
end
