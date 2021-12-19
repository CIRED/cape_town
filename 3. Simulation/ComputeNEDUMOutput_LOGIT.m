function [jobSimul,R,peopleInit,peopleCenter,housingSupply,dwellingSize,R_mat] = ComputeNEDUMOutput_LOGIT(Uo,param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParameter,interestRate,income,multiProbaGroup,transportCost,transportCostRDP,coeffLand,poly,amenities,solus_R, solus_Q, functionName, typeHousing)
% Works both for formal or informal housing


% carto2 = @(x) scatter(grid.xCoord, grid.yCoord, 400, x, '.');
basic_q_formal = param.basic_q;


%% Bid rents and dwelling sizes
% Calculation depends on whether we have basic_need / mini_lotsize

run(functionName)


% Income group in each location
proba = (R_mat == nanmax(R_mat)); 
proba(~isnan(multiProbaGroup)) = multiProbaGroup(~isnan(multiProbaGroup));
limit = (transTemp.incomeNetOfCommuting > 0) & (proba>0) & (~isnan(transTemp.incomeNetOfCommuting)) & (R_mat>0);
proba = proba.*limit;

[~, whichGroup] = max(R_mat);
whichGroup(~isnan(multiProbaGroup(1,:))) = sum(repmat([1:param.numberIncomeGroup]', 1, sum(~isnan(multiProbaGroup(1,:)))).*proba(:,~isnan(multiProbaGroup(1,:))));
temp = [0:size(transTemp.incomeNetOfCommuting,2)-1] * size(transTemp.incomeNetOfCommuting,1);
whichGroupTemp = whichGroup + temp; 

R = R_mat(whichGroupTemp);
dwellingSize = dwellingSize(whichGroupTemp);

%% Housing Construction 

switch typeHousing
    
case 'formal'
    housingSupply = CalculateHousingSupplyFormal(R, option, limitHousing, constructionParameter, param, agriculturalRent, referenceRent, interestRate);
case 'backyard'
    housingSupply = CalculateHousingSupplyBackyard(R, grid, param, basic_q_formal, income, transportCostRDP);
case 'informal'
    if option.doubleStoreyShacks == 0
        housingSupply = 1000000 .* ones(size(whichGroupTemp));
        housingSupply(R == 0) = 0;
    elseif option.double_storey_shacks == 1
        housingSupply = CalculateHousingSupplySettlement(R, grid, param, poly, income, transportCost, proba);
    end
end


peopleInit = housingSupply./dwellingSize .* (sum(limit,1)>0);
peopleInit(isnan(peopleInit)) = 0;
peopleInitLand = peopleInit .* coeffLand .* grid.sizeSquare^2;

peopleCenter = (ones(size(R_mat,1),1)*peopleInitLand) .* proba;
peopleCenter(isnan(peopleCenter)) = 0;
jobSimul = sum(peopleCenter,2)';
    
switch typeHousing 
case'formal'
    R = max(R,agriculturalRent);
end

end

