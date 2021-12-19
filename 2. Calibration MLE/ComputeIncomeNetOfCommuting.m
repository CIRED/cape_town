function [incomeNetOfCommuting, modalShares, ODflows, averageIncome] = ...
    ComputeIncomeNetOfCommuting(param, trans, grid, poly, data, lambda, incomeCenters, areaInterp, year)
% Compute income net of commuting


annualToHourly = 1./(8*20*12);

timeCost = trans.timeCost;
timeCost(isnan(timeCost)) = 10^2;
monetaryCost = trans.monetaryCost(:,:,:,year) .* annualToHourly;
monetaryCost(isnan(monetaryCost)) = 10^3 .* annualToHourly;
incomeCenters = incomeCenters .* annualToHourly;

if areaInterp == 'GR' % grid
    xInterp = grid.xCoord;
    yInterp = grid.yCoord;
elseif areaInterp == 'SP'
    xInterp = data.spX;
    yInterp = data.spY;
end

modalShares = zeros(length(incomeCenters), size(trans.timeCost,2), trans.numberModes, param.numberIncomeGroup);
ODflows = zeros(length(incomeCenters), size(trans.timeCost,2), param.numberIncomeGroup);
incomeNetOfCommuting = zeros(param.numberIncomeGroup, size(trans.timeCost,2));
averageIncome = zeros(param.numberIncomeGroup, size(trans.timeCost,2));

for j = 1:param.numberIncomeGroup
        
        %% Household size varies with transport costs
        householdSize = param.householdSizeTransport(j);
        whichCenters = incomeCenters(:,j) > -100000;
        incomeCentersGroup = incomeCenters(whichCenters,j);
           
    %% Transport costs and employment allocation
        transportCostModes = double(householdSize.*monetaryCost(whichCenters,:,:) + timeCost(whichCenters,:,:) .* repmat(incomeCentersGroup, 1, size(timeCost,2), size(timeCost,3)));
        
        % Value max is to prevent the exp to diverge to infinity (in matlab: exp(800) = Inf)
        valueMax = min(lambda .* transportCostModes,[],3) - 500;
        
        % Modal shares
        modalShares(whichCenters,:,:,j) = exp(- lambda.*transportCostModes + valueMax) ./ nansum(exp(-lambda.*transportCostModes + valueMax), 3);
        
        % Transport costs
        transportCost = - 1./lambda .* (log(nansum(exp(-lambda.*transportCostModes + valueMax),3)) - valueMax);

        % minIncome is also to prevent diverging exponentials
        minIncome = nanmax(lambda .* (incomeCentersGroup - transportCost)) - 700;

        % OD flows
        ODflows(whichCenters,:,j) = exp(lambda .* (incomeCentersGroup - transportCost) - minIncome) ./ nansum(exp(lambda .* (incomeCentersGroup - transportCost) - minIncome), 1);

        % Income net of Commuting /!\ THIS FORMULA IS WRONG, see below
        % incomeNetOfCommuting(j,:) = nansum(ODflows(whichCenters,:,j) .* (incomeCentersGroup - transportCost),1);
        
        % Income net of commuting (correct formula)
        incomeNetOfCommuting(j,:) = 1./lambda .* (log(nansum(exp(lambda.*(incomeCentersGroup - transportCost) - minIncome), 1)) + minIncome);
        
        % Average income earned by a worker
        averageIncome(j,:) = nansum(ODflows(whichCenters,:,j) .* (incomeCentersGroup));
             
end

incomeNetOfCommuting = incomeNetOfCommuting ./ annualToHourly;
averageIncome = averageIncome ./ annualToHourly;

end
