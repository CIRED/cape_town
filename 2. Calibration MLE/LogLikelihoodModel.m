function [scoreTotal, scoreAmenities, scoreDwellingSize, scoreIncomeSorting, scoreHousing, parametersAmenities, modelAmenities, parametersHousing] = LogLikelihoodModel(X0, Uo2, incomeNetOfCommuting, groupLivingMatrix, dataDwellingSize, selectedDwellingSize, xData, yData, dataRent, selectedRents, dataHouseholdDensity, selectedDensity, predictorsAmenitiesMatrix, tableRegression, variablesRegression, CalculateDwellingSize, ComputeLogLikelihood, optionRegression)
% Function to estimate the total likelihood of the model given the
% parameters

beta = X0(1);
basicQ = X0(2);
Uo = [Uo2; X0(3:4)];

% Errors on the amenity

% Calculate amenities as a residual
residualAmenities = log(Uo) - log((1-beta).^(1-beta).*beta.^beta.*(incomeNetOfCommuting(:,selectedRents) - basicQ.*dataRent(selectedRents)) ./ (dataRent(selectedRents).^beta));
residualAmenities = nansum(residualAmenities .* groupLivingMatrix(:,selectedRents));
residualAmenities(abs(imag(residualAmenities)) > 0) = NaN;
residualAmenities(residualAmenities == 0) = NaN;

% First possibility: amenities follow a log-normal law
% meanAmenities = nanmean(exp(residualAmenities));
% residualAmenities = residualAmenities - log(meanAmenities); % residualAmenities is the log of amenities 
% errorAmenities = residualAmenities;
% Uo = Uo ./ meanAmenities;

% Second possibility: residual for the regression of amenities follow a log-normal law
if optionRegression == 0
                    
    % Here regression as a matrix division (much faster)
    parametersAmenities = predictorsAmenitiesMatrix(~isnan(residualAmenities'),:) \ real(residualAmenities(~isnan(residualAmenities))');
    errorAmenities = real(residualAmenities(~isnan(residualAmenities))') - (predictorsAmenitiesMatrix(~isnan(residualAmenities'),:)*parametersAmenities); 
    
elseif optionRegression == 1
            
    % Compute regression with fitglm (longer)
    % Can only work if length(lists) = 1
    tableRegression.residu = real(residualAmenities');
    parametersAmenities = predictorsAmenitiesMatrix(~isnan(residualAmenities'),:) \ real(residualAmenities(~isnan(residualAmenities))');
    modelSpecification = ['residu ~ ', sprintf('%s + ', variablesRegression{1:end-1}), variablesRegression{end}];
    modelAmenities = fitglm(tableRegression, modelSpecification);
    errorAmenities = modelAmenities.Residuals.Raw;
            
end      

% Error on allocation of income groups

% First possibility: error_income is just the number of errors
% sorting = (1-beta).^(1-beta).*beta.^beta.*(incomeNetOfCommutingSelectedRents - basicQ .* dataRentMatrix) ./ (dataRentMatrix.^beta) < Uo./exp(residualAmenities);
% errorIncomeSorting = sum(sum(sorting(:,~isnan(residualAmenities)).*(~logical(incomeGroupSelectedRents(:,~isnan(residualAmenities))))) > 0) ./ sum(~isnan(residualAmenities));
            
% Second possibility: log-likelihood of a logit model on the location of income groups
griddedRents = InterpolateRents(beta, basicQ, incomeNetOfCommuting);
bidRents = exp(griddedRents(log(Uo) - residualAmenities, log(incomeNetOfCommuting(:,selectedRents))));

% Estimation of the scale parameter by maximization of the log-likelihood
selectedBidRents = nansum(bidRents) > 0;
incomeGroupSelectedRents = groupLivingMatrix(:,selectedRents);
likelihoodIncomeSorting = @(scaleParam) - (nansum(nansum(bidRents(:,selectedBidRents)./scaleParam .* incomeGroupSelectedRents(:,selectedBidRents))) - nansum(log(nansum(exp(bidRents(:,selectedBidRents)./scaleParam),1))));
scoreIncomeSorting = - likelihoodIncomeSorting(10000);
% optionsOptim = optimset('Display', 'notify');
% lowerBounds = 0;
% upperBounds = 10^10;
% initScale = 10^9;
% [scaleParameter, errorIncomeSorting, exitFlag] = fmincon(likelihoodIncomeSorting, initScale, [], [], [], [], lowerBounds, upperBounds, [], optionsOptim);
% if exitFlag < 1
%     disp('Scale parameter not found')
% end
% if exitFlag ~= 1 && optionRegression == 1
%     disp('Scale parameter not found - final state')
% end
% scoreIncomeSorting = - errorIncomeSorting; % we take the opposite of what we used in the minimization algorithm
    
% Errors on the dwelling sizes

% Option 1: real data, real sorting
% dwellingSizeTemp = CalculateDwellingSize(beta, basicQ, incomeNetOfCommuting(:, selectedDwellingSize), dataRent(selectedDwellingSize));
% dwellingSize = sum(dwellingSizeTemp .* groupLivingMatrix(:, selectedDwellingSize));

% Option 2: simulated rent, simulated sorting
% [simulatedBidRents, incomeGroupSimulation] = max(bidRents);
% dwellingSizeTemp = CalculateDwellingSize(beta, basicQ, incomeNetOfCommuting(:, selectedDwellingSize), simulatedBidRents(selectedDwellingSize(selectedRents)));
% dwellingSize = zeros(1, length(incomeGroupSimulation));
% for i = 1:3
%     dwellingSize(incomeGroupSimulation == i) = dwellingSizeTemp(i,incomeGroupSimulation == i);
% end

% Option 3: simulated rent, real sorting
simulatedRents = sum(bidRents(:,selectedDwellingSize(selectedRents)) .* groupLivingMatrix(:, selectedDwellingSize));
dwellingSize = CalculateDwellingSize(beta, basicQ, sum(incomeNetOfCommuting(:, selectedDwellingSize) .* groupLivingMatrix(:, selectedDwellingSize)), simulatedRents);
% Define errors
errorDwellingSize = log(dwellingSize) - log(dataDwellingSize(selectedDwellingSize));
      
% Errors on household density
% dataDensitySelected = dataHouseholdDensity(selectedRents);
% selectedResized = ~isnan(simulatedRents) & selectedDensity(selectedRents);
% logHousingSupply = log(dataDensitySelected(selectedResized)) + log(dwellingSize(selectedResized));
% parametersHousing = [ones(1,sum(selectedResized)); log(simulatedRents(selectedResized))]' \ logHousingSupply';
% errorHousing = logHousingSupply - parametersHousing' * [ones(1,sum(selectedResized)); log(simulatedRents(selectedResized))];
% scoreHousing = ComputeLogLikelihood(sqrt(nansum(errorHousing.^2) ./ sum(~isnan(errorAmenities))), errorHousing);
scoreHousing = 0;
parametersHousing = 0;

% Scores

scoreDwellingSize = ComputeLogLikelihood(sqrt(nansum(errorDwellingSize.^2) ./ sum(~isnan(errorDwellingSize))), errorDwellingSize);
scoreAmenities = ComputeLogLikelihood(sqrt(nansum(errorAmenities.^2) ./ sum(~isnan(errorAmenities))), errorAmenities);

scoreTotal = scoreAmenities + scoreDwellingSize + scoreIncomeSorting; %+ scoreHousing;

end


function [ utility ] = utilityFromRents(Ro, income, basic_q, beta)
% Precalculation of the utility
%   [ utili ] = utilityFromRents(rent,income,basic_q,beta) estimates
%   utility for a rent "Ro" and an income (net of transport costs)
%   "revenu". 
% Without amenities 

utility = (1-beta)^(1-beta)*beta^beta...
          .*(income - basic_q.*Ro)...
          ./(Ro.^beta);
utility((income - basic_q.*Ro) < 0) = 0;

% If the income is 0 then utility is 0 (we need to specify it to avoid
% NaNs later)
utility(income == 0) = 0;


end

% Precalculations for rents, as a function
% The output of the function is a griddedInterpolant object, that gives the
% log of rents as a function of the log utility and the log income
function griddedRents = InterpolateRents(beta, basic_q, netIncome)
    
        % Decomposition for the interpolation (the more points, the slower the code)
        decompositionIncome = [10.^([-9,-4:0.5:-2]), 0.03, 0.06:0.02:1.4,1.5:0.1:2.5, 4:2:10, 20, 10^9];
        decompositionRent = [10.^([-9,-4,-3,-2]), 0.02:0.01:0.79, 0.8:0.02:0.96, 0.98]; 
        
        % Min and Max values for the decomposition
        choiceIncome = 100000 .* decompositionIncome ;
        incomeMatrix = repmat(choiceIncome,length(choiceIncome),1)';
        choiceRent = choiceIncome./basic_q; % the maximum rent is the rent for which u = 0
        rentMatrix = choiceRent' * decompositionRent;

        utilityMatrix = utilityFromRents(rentMatrix, incomeMatrix, basic_q, beta);
        solusRentTemp = @(x,y) griddata(incomeMatrix, utilityMatrix, rentMatrix.^beta,x,y);
        
        % Redefine a grid (to use griddedInterpolant)
        utilityVectLog = -1:0.1:log(max(max(10.*netIncome)));
        incomeLog = (-1:0.2:log(max(max(10.*netIncome))))';
        rentLog = 1/beta.*log(solusRentTemp(exp(incomeLog), exp(utilityVectLog)));
        griddedRents = griddedInterpolant({utilityVectLog, incomeLog}, rentLog, 'linear', 'none');
        
end