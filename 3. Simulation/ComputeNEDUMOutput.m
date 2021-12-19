function [jobSimul,R,peopleInit,peopleCenter,housingSupply,dwellingSize,R_mat] = ComputeNEDUMOutput(Uo,param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParameter,interestRate,income,multi_proba,transportCost,transportCostRDP,coeffLand,poly,amenities,solus_R, solus_Q, functionName, typeHousing)
% Works both for formal or informal housing


% carto2 = @(x) scatter(grid.xCoord, grid.yCoord, 400, x, '.');
basic_q_formal = param.basic_q;


%% Bid rents and dwelling sizes
% Calculation depends on whether we have basic_need / mini_lotsize

run(functionName)


%% Estimation of the share of households living in x

% These lines mean that the dwelling size is constant for each employment
% center at a location x
dwellingSizeMatrix = ones(size(Uo')) * dwellingSize;
spendingsHousingMatrix = ones(size(Uo')) * (dwellingSize .* R);

Z = (transTemp.incomeNetOfCommuting - spendingsHousingMatrix);
Z(Z<=0) = 0;
utility = ComputeUtilityWithAmenities(Z, dwellingSizeMatrix, param, amenities,income, 0);
maxUtility = Uo'; % Maximum utility is constant for each job center
maxUtilityMat = maxUtility * ones(size(R));
utility = (abs(utility)).^0.01;
maxUtilityMat = (abs(maxUtilityMat)).^0.01; %1000 si 0.01%2000 si 0.005;11000 si 0.001;110000 si 0.0001
lambda = param.lambda * ones(size(utility));%lambda(1:20,:)=lambda(1:20,:)*0.96;%for 0.01
probaLog = -(maxUtilityMat ./ utility - 1).*lambda;


pixelsNoProba = isnan(probaLog);
probaLog(pixelsNoProba) = -100000;

maxProbaLog = max(probaLog,[],1);
maxProbaLog(isnan(maxProbaLog)) = 0;
maxProbaLog(isinf(maxProbaLog)) = 0;
maxProbaLog = ones(size(R_mat,1),1)*maxProbaLog;
probaLog = probaLog - maxProbaLog;

% Number of jobs
probaLog = probaLog + log(double(multi_proba));

if sum(sum(~isreal(probaLog)))>=1
    disp('Error: complex numbers in proba - click pause')
    %pause
end

% Exponential form
proba = exp(probaLog);
proba(Z<=0) = 0;
proba(pixelsNoProba) = 0;
proba(R_mat<=0) = 0;

% Normalization of the proba
probaNormalized = sum(proba,1);
proba = proba./(ones(size(R_mat,1),1)*probaNormalized);
proba(((ones(size(R_mat,1),1)*probaNormalized))==0) = 0;
proba = single(proba);

%% Housing Construction 

switch typeHousing
    
case 'formal'
    housingSupply = CalculateHousingSupplyFormal(R, option, limitHousing, constructionParameter, param, agriculturalRent, referenceRent, interestRate);
case 'backyard'
    housingSupply = CalculateHousingSupplyBackyard(R, grid, param, basic_q_formal, income, transportCostRDP);
case 'informal'
    if option.doubleStoreyShacks == 0
        housingSupply = 1000000 .* ones(size(whichMaxTemp));
        housingSupply(R == 0) = 0;
    elseif option.double_storey_shacks == 1
        housingSupply = CalculateHousingSupplySettlement(R, grid, param, poly, income, transportCost, proba);
    end
end

limite = (income>transportCost) & (proba>0) & (~isnan(transportCost)) & (R_mat>0);
proba = proba.*limite;

peopleInit = housingSupply./dwellingSize .* (sum(limite,1)>0);
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


function [ utility ] = ComputeUtilityWithAmenities(Z, dwellingSize, param, amenities, income, Ro)
% Utility depending of rent "Ro" and income
% With amenities

if Ro == 0
    utility = Z.^(param.alpha).*(dwellingSize - param.basic_q).^param.beta;
else
    Ro = (ones(length(income(1,:)),1)*Ro)';
      utility=param.alpha^param.alpha*param.beta^param.beta...
            .*sign(income - param.basic_q.*Ro)...
            .*abs(income - param.basic_q.*Ro)...
            ./(Ro.^param.beta);
end

% Amenities are a factor of utility
utility = utility.*amenities;

% If the income is 0 then utility is 0 (we need to specify it to avoid
% NaNs later)
utility(income==0) = 0;

end