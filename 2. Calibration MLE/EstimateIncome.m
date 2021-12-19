 function [modalSharesTot, incomeCentersSave, timeDistribution, distanceDistribution] = ...
    EstimateIncome(param, trans, poly, data, listLambda)
% Solve for income per employment centers for different values of lambda

disp('Estimation of local incomes, and lambda parameter')

annualToHourly = 1./(20*12*8);
bracketsTime = [0, 15, 30, 60, 90, nanmax(nanmax(nanmax(trans.timeOutput)))];
bracketsDistance = [0, 5, 10, 15, 20, 25, 30, 35, 40, 200];

timeCost = trans.timeCost(:,data.sp2011CapeTown, :);
timeCost(isnan(timeCost)) = 10^2;
monetaryCost = trans.monetaryCost(:,data.sp2011CapeTown,:,1) .* annualToHourly;
monetaryCost(isnan(monetaryCost)) = 10^3 .* annualToHourly;
transportTimes = trans.timeOutput(:,data.sp2011CapeTown) ./ 2; % Round-trip to direct
transportDistances = trans.distanceOutput(:,data.sp2011CapeTown, 1);

modalSharesTot = zeros(5, length(listLambda));
incomeCentersSave = zeros(length(poly.jobsCenters(:,1,1)), param.numberIncomeGroup, length(listLambda));
timeDistribution = zeros(length(bracketsTime) - 1, length(listLambda));
distanceDistribution = zeros(length(bracketsDistance) - 1, length(listLambda));

for i = 1:length(listLambda)

    lambda = listLambda(i);
    
    fprintf('    -> Estimating for lambda = %g', lambda);
    disp(' ')
    
    incomeCentersAll = -Inf .*  ones(length(poly.jobsCenters(:,1,1)), param.numberIncomeGroup);
    modalSharesGroup = zeros(5, param.numberIncomeGroup);
    timeDistributionGroup = zeros(length(bracketsTime) - 1, param.numberIncomeGroup);
    distanceDistributionGroup = zeros(length(bracketsDistance) - 1, param.numberIncomeGroup);

    for j = 1:param.numberIncomeGroup
        
        % Household size varies with transport costs
        householdSize = param.householdSizeTransport(j);
    
        averageIncomeGroup = poly.averageIncomeGroup(1,j) .* annualToHourly;
        
        fprintf('        incomes for group %g', j);
        
        whichJobsCenters = poly.jobsCenters(:,j,1) > 600; 
        popCenters = poly.jobsCenters(whichJobsCenters,j,1);
        popResidence = data.sp2011IncomeDistributionNClass(data.sp2011CapeTown,j)' .* sum(poly.jobsCenters(whichJobsCenters,j,1)) ./ sum(data.sp2011IncomeDistributionNClass(data.sp2011CapeTown,j));
                
        % Function to solve
        funSolve = @(incomeCentersTemp) fun0(incomeCentersTemp, averageIncomeGroup, popCenters, popResidence, monetaryCost(whichJobsCenters,:,:) .* householdSize, timeCost(whichJobsCenters,:,:) .* householdSize, lambda);

        % Solver
        % lowerBounds = zeros(length(incomeCenters0), 1);
        % upperBounds = Inf.*ones(length(incomeCenters0), 1);
        % optionsOptim = optimset('Display', 'iter');
        % [incomeCenters, scoreTot, exitFlag] = fmincon(funSolve, incomeCenters0, [], [], [], [], lowerBounds, upperBounds, [], optionsOptim);        
    
        maxIter = 700;
        tolerance = 0.1;
        if j == 1
            factorConvergenge = 0.008; % 0.0003;
        elseif j == 2
            factorConvergenge = 0.005;
        else
            factorConvergenge = 0.0005;
        end 
        
        iter = 0;
        error = zeros(length(popCenters), maxIter);
        scoreIter = zeros(1, maxIter);
        errorMax = 1;
        
        % Initializing the solver
        incomeCenters = zeros(sum(whichJobsCenters), maxIter);
        incomeCenters(:,1) =  averageIncomeGroup .* (popCenters ./ mean(popCenters)).^(0.1);
        error(:,1) = funSolve(incomeCenters(:,1));

        
        while (iter <= maxIter) && (errorMax > tolerance)
            
            iter = iter + 1;
            incomeCenters(:,iter) = incomeCenters(:,max(iter - 1,1)) + factorConvergenge .* averageIncomeGroup .* error(:, max(iter - 1,1))./popCenters;
            % incomeCenters(incomeCenters(:,iter) < 0, iter) = 0;
            error(:,iter) = funSolve(incomeCenters(:,iter));
            errorMax = max(abs(error(:,iter) ./ popCenters));
            scoreIter(iter) = mean(abs(error(:,iter) ./ popCenters));
            
        end
        
        if iter > maxIter
            [scoreBest, bestSolution] = min(scoreIter);
            incomeCenters(:,iter) = incomeCenters(:, bestSolution);
            fprintf(' - max iteration reached - mean error %d', scoreBest)
        else
            fprintf(' - computed - max error %d', errorMax)
        end
        
        disp(' ')
        
        
        incomeCentersRescaled = incomeCenters(:,iter) .* averageIncomeGroup ./ ((sum(incomeCenters(:,iter) .* popCenters) ./ sum(popCenters)));
        modalSharesGroup(:,j) = modalShares(incomeCentersRescaled, popCenters, popResidence, monetaryCost(whichJobsCenters,:,:) .* householdSize, timeCost(whichJobsCenters,:,:) .* householdSize, lambda);
        incomeCentersAll(whichJobsCenters,j) = incomeCentersRescaled;
        
        timeDistributionGroup(:,j) = computeDistributionCommutingTimes(incomeCentersRescaled, popCenters, popResidence, monetaryCost(whichJobsCenters,:,:) .* householdSize, timeCost(whichJobsCenters,:,:) .* householdSize, transportTimes(whichJobsCenters,:), bracketsTime, lambda);
        distanceDistributionGroup(:,j) = computeDistributionCommutingDistances(incomeCentersRescaled, popCenters, popResidence, monetaryCost(whichJobsCenters,:,:) .* householdSize, timeCost(whichJobsCenters,:,:) .* householdSize, transportDistances(whichJobsCenters,:), bracketsDistance, lambda);

        
    end
    
    modalSharesTot(:,i) = sum(modalSharesGroup, 2) ./ sum(sum(modalSharesGroup));
    incomeCentersSave(:,:,i) = incomeCentersAll(:,:) ./ annualToHourly;
    timeDistribution(:,i) = sum(timeDistributionGroup,2) ./ sum(sum(timeDistributionGroup));
    distanceDistribution(:,i) = sum(distanceDistributionGroup,2) ./ sum(sum(distanceDistributionGroup));

end


end

function [score] = ...
    fun0(incomeCenters, meanIncome, popCenters, popResidence, monetaryCost, timeCost, lambda)
% Computes error in employment allocation

% In order to control the average income of each group
incomeCentersFull = incomeCenters .* meanIncome ./ ((sum(incomeCenters .* popCenters) ./ sum(popCenters)));

% Transport costs and employment allocation
transportCostModes = double(monetaryCost + timeCost .* repmat(incomeCentersFull, 1, size(timeCost,2), size(timeCost,3)));

% Value max is to prevent the exp to diverge to infinity (in matlab: exp(800) = Inf)
valueMax = min(lambda .* transportCostModes,[],3) - 500;

% Transport costs
transportCost = - 1./lambda .* (log(nansum(exp(-lambda.*transportCostModes + valueMax),3)) - valueMax);

% minIncome is also to prevent diverging exponentials
minIncome = max(nanmax(lambda .* (incomeCentersFull - transportCost))) - 500;

% Differences in the number of jobs
score = popCenters - nansum( exp(lambda.*(incomeCentersFull - transportCost) - minIncome) ./ nansum(exp(lambda .* (incomeCentersFull - transportCost) - minIncome)) .* popResidence, 2);

end


function [modalSharesTot] = ...
    modalShares(incomeCenters, popCenters, popResidence, monetaryCost, timeCost, lambda)
% Computes total modal shares

% Transport cost by modes
transportCostModes = double(monetaryCost + timeCost .* repmat(incomeCenters, 1, size(timeCost,2), size(timeCost,3)));

% Value max is to prevent the exp to diverge to infinity (in matlab: exp(800) = Inf)
valueMax = min(lambda .* transportCostModes,[],3) - 500;

% Compute modal shares
modalSharesTemp = exp(- lambda .* transportCostModes + valueMax) ./ nansum(exp(- lambda .* transportCostModes + valueMax),3);

% Multiply by OD flows
transportCost = - 1./lambda .* (log(nansum(exp(-lambda.*transportCostModes + valueMax),3)) - valueMax);

% minIncome is also to prevent diverging exponentials
minIncome = max(nanmax(lambda .* (incomeCenters - transportCost))) - 500;

% Total modal shares
modalSharesTot = nansum(nansum(modalSharesTemp .* repmat(exp(lambda .* (incomeCenters - transportCost) - minIncome) ./ nansum(exp(lambda .* (incomeCenters - transportCost) - minIncome)), 1, 1, 5) .* popResidence, 2));
modalSharesTot(:) = permute(modalSharesTot, [3,1,2]);

end


function [nbCommuters] = ...
    computeDistributionCommutingTimes(incomeCenters, popCenters, popResidence, monetaryCost, timeCost, transportTime, bracketsTime, lambda)
% Computes total modal shares

% Transport cost by modes
transportCostModes = double(monetaryCost + timeCost .* repmat(incomeCenters, 1, size(timeCost,2), size(timeCost,3)));

% Value max is to prevent the exp to diverge to infinity (in matlab: exp(800) = Inf)
valueMax = min(lambda .* transportCostModes,[],3) - 600;

% Compute modal shares
modalSharesTemp = exp(- lambda .* transportCostModes + valueMax) ./ nansum(exp(- lambda .* transportCostModes + valueMax),3);

% Multiply by OD flows
transportCost = - 1./lambda .* (log(nansum(exp(-lambda.*transportCostModes + valueMax),3)) - valueMax);

% minIncome is also to prevent diverging exponentials
minIncome = max(nanmax(lambda .* (incomeCenters - transportCost))) - 600;

% Total distribution of times
nbCommuters = zeros(length(bracketsTime) - 1, 1);
for k = 1:length(bracketsTime)-1
    which = transportTime > bracketsTime(k) & transportTime <= bracketsTime(k + 1);
    nbCommuters(k) = nansum(nansum(nansum(which .* modalSharesTemp .* repmat(exp(lambda .* (incomeCenters - transportCost) - minIncome) ./ nansum(exp(lambda .* (incomeCenters - transportCost) - minIncome)), 1, 1, 5) .* popResidence, 2)));
end

end


function [nbCommuters] = ...
    computeDistributionCommutingDistances(incomeCenters, popCenters, popResidence, monetaryCost, timeCost, transportDistance, bracketsDistance, lambda)
% Computes total modal shares

% Transport cost by modes
transportCostModes = double(monetaryCost + timeCost .* repmat(incomeCenters, 1, size(timeCost,2), size(timeCost,3)));

% Value max is to prevent the exp to diverge to infinity (in matlab: exp(800) = Inf)
valueMax = min(lambda .* transportCostModes,[],3) - 500;

% Compute modal shares
modalSharesTemp = exp(- lambda .* transportCostModes + valueMax) ./ nansum(exp(- lambda .* transportCostModes + valueMax),3);

% Multiply by OD flows
transportCost = - 1./lambda .* (log(nansum(exp(-lambda.*transportCostModes + valueMax),3)) - valueMax);

% minIncome is also to prevent diverging exponentials
minIncome = max(nanmax(lambda .* (incomeCenters - transportCost))) - 500;

% Total distribution of times
nbCommuters = zeros(length(bracketsDistance) - 1, 1);
for k = 1:length(bracketsDistance)-1
    which = transportDistance > bracketsDistance(k) & transportDistance <= bracketsDistance(k + 1);
    nbCommuters(k) = nansum(nansum(nansum(which .* modalSharesTemp .* repmat(exp(lambda .* (incomeCenters - transportCost) - minIncome) ./ nansum(exp(lambda .* (incomeCenters - transportCost) - minIncome)), 1, 1, 5) .* popResidence, 2)));
end

end