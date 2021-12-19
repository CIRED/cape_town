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




function initialState = RunEquilibriumSolverNEDUM(yearEquilibrium, trans, option, land, grid, macro, param, poly, Uo_perso, housingFixed)
% Solver with a Stone-Geary utility function, n income classes and informal
% housing (settlement and backyard shacks)

% Income for the year (varies in time)
incomeMatrix = InterpolateIncomeEvolution(macro,param,option,grid, poly,yearEquilibrium);
% Average income for the year of the simulation
incomeAverage = ppval(macro.splineIncome,yearEquilibrium);

% Final price for transport (not only pt)
transportCost = InterpolateGeneralizedTransportCostEvolution(trans,param,yearEquilibrium);
transportCost(transportCost==0) = NaN;
% Cost associated with the time spent commuting
transportCostTime = InterpolateCostEvolution(trans.yearTransport,trans.timeCost,param,yearEquilibrium);
transportCostTime(transportCostTime==0) = NaN;
interestRate = InterpolateInterestRateEvolution(macro, yearEquilibrium); 

%Population
population = InterpolatePopulationEvolution(macro,yearEquilibrium); %*coeff_temp;
totalRDP = ppval(macro.splineRDP, yearEquilibrium);

%Construction coefficient
constructionParam = InterpolateCoefficientConstruction(param,macro,incomeAverage);

% Evolution of coeffLand
land.coeffLand = InterpolateLandCoefficientEvolution(land, option, param, yearEquilibrium);
land.numberPropertiesRDP = ppval(land.splineEstimateRDP,yearEquilibrium)';

%Limit of housing construction
housingLimit = InterpolateHousingLimitEvolution(land,option,param,yearEquilibrium);
if nargin<10
    housingFixed = zeros(1,length(grid.distanceCBD));
end

% Transaction cost is the rent at the city limit (per year)
agriculturalRent = InterpolateAgriculturalRentEvolution(param,macro,incomeAverage);
rentReference = agriculturalRent;

% Tax outside the urban edge
param.taxUrbanEdgeMat = zeros(1, length(grid.distanceCBD));
if option.taxOutUrbanEdge == 1
    param.taxUrbanEdgeMat(land.urbanEdge == 0) = param.taxUrbanEdge .* interestRate;
end


% Computation of the initial state
if option.ownInitializationSolver
    initialState = ComputeEquilibrium(option,land,grid,macro,param,yearEquilibrium,rentReference,housingLimit,incomeMatrix,incomeAverage,transportCost,interestRate,population,agriculturalRent,constructionParam,poly,Uo_perso,transportCostTime,housingFixed, totalRDP);
else
    initialState = ComputeEquilibrium(option,land,grid,macro,param,yearEquilibrium,rentReference,housingLimit,incomeMatrix,incomeAverage,transportCost,interestRate,population,agriculturalRent,constructionParam,poly,1,transportCostTime,housingFixed, totalRDP);
end

end



%% Solver for the initial state
function initialState = ComputeEquilibrium(option,land,grid,macro,param,yearEquilibrium,referenceRent,limitHousing,incomeMatrix,incomeAverage,transportCost,interestRate,population,agriculturalRent,constructionParam,poly,Uo_init,transportCostTime,housingFixed, totalRDP)
nargin_bak = nargin;

if nargin_bak >= 18
    param.housing_in = housingFixed;
end

% Some options
printAdvanceBar = 1;

% Maximum number of iterations in the solver
maxIteration = param.maxIteration;
% Resolution precision
precision = param.precision;

fprintf('Static simulation %g\n',yearEquilibrium + param.yearBegin);

%% Preparation of the variables


% We add the depreciation rate and the interest rate
interestRate = interestRate + param.depreciationRate;

param.incomeYearReference = macro.incomeYearReference;

% Employment centers
numberJobsCenter = interp1(poly.year,poly.jobsCenters(1:length(poly.year),:),yearEquilibrium + param.yearBegin);

% Income of each class
averageIncome = interp1(poly.year, poly.averageIncome(1:length(poly.year), :), yearEquilibrium + param.yearBegin);
poly.incomeMult = averageIncome ./ ppval(macro.splineIncome, yearEquilibrium);

% Class of each center and housing types
formalDummy = poly.formal;
backyardDummy = poly.backyard;
settlementDummy = poly.settlement;
poly.formal = zeros(1, length(averageIncome));
poly.backyard = zeros(1, length(averageIncome));
poly.settlement = zeros(1, length(averageIncome));
poly.incomeGroup = poly.incomeGroup(1,:);

for i = 1:param.numberIncomeGroup
    if formalDummy(i) == 1
        poly.formal(poly.incomeGroup == i) = 1;   
    end
    if backyardDummy(i) == 1
        poly.backyard(poly.incomeGroup == i) = 1;
    end
    if settlementDummy(i) == 1
        poly.settlement(poly.incomeGroup == i) = 1;
    end
end

% Ajust the population to remove the population in RDP
ratio = population ./ sum(numberJobsCenter,2);
numberJobsCenter = numberJobsCenter * ratio;
totalRDP = totalRDP * ratio;
numberJobsCenter(poly.incomeGroup == 1) = numberJobsCenter(poly.incomeGroup == 1) - totalRDP .* numberJobsCenter(poly.incomeGroup == 1)./sum(numberJobsCenter(poly.incomeGroup == 1));
xJobCenter = poly.xCenter;
yJobCenter = poly.yCenter;
employmentCenters = [numberJobsCenter;xJobCenter;yJobCenter];
employmentCenters = single(employmentCenters);

% We define multi_proba (useful in the solver) here because it is faster
multi_proba = (employmentCenters(1,:)'*ones(1,size(grid.distanceCBD,2)));

% Commuting price for RDP households (useful for Backyarding)
% Note: Households in RDP are allocated randomly to job centers
transportCostRDP = sum(numberJobsCenter(poly.incomeGroup == 1)'.* transportCost(poly.incomeGroup == 1,:),1) ./ sum(numberJobsCenter(poly.incomeGroup == 1)) ;


%% Amenities

% Loading amenities
amenities = land.amenities;

% We transform amenities in a matrix with as many lines as employment
% centers
amenities = ones(size(incomeMatrix,1),1)*amenities;

%% What "bid_rents" function should we use

if param.basic_q > 0 && param.miniLotSize > 0
    functionName = 'ComputeBidRentsBasicNeedMinimumLotSize.m';
elseif param.basic_q > 0
    functionName = 'ComputeBidRentsBasicNeed.m';
elseif param.miniLotSize > 0
    functionName = 'ComputeBidRentsMinimumLotSize.m';
end


%% Pre-calculation of the utility / rent relationship

% Precalculations for rents

% Function that estimates the utility
utilityRent = @(Ro,income) ComputeUtilityFromRent(Ro, income, param.basic_q, param);

% Decomposition for the interpolation (the more points, the slower the code)
decompositionRent = [10.^([-9,-4,-3,-2]),0.02:0.015:0.08,0.1:0.02:1]; 
decompositionIncome = [10.^([-9,-4:0.5:-2]),0.03,0.06:0.02:1];

% Min and Max values for the decomposition
incomeVector = max(max(incomeMatrix)) .* decompositionIncome ;
incomeMat = repmat(incomeVector,length(incomeVector),1)';

if param.basic_q == 0
     rentVector = 800*12.*decompositionRent;
     rentMatrix = repmat(rentVector,length(incomeVector),1);
else
     rentVector = incomeVector./param.basic_q; % the maximum rent is the rent for which u = 0
     rentMatrix = rentVector' * decompositionRent;
end
XX = incomeMat;
YY_R = utilityRent(rentMatrix,incomeMat);
ZZ_R = rentMatrix;

% The function that will be used afterwards
solus_R = @(x,y) (griddata(XX,YY_R,ZZ_R.^param.beta,x,y)).^(1./param.beta);


%%  Precalculations for dwelling sizes
    
% Function that estimates the utility
utilitySize = @(q,income) ComputeUtilityFromDwellingSize(q, income, param.basic_q, param);

% Decomposition for the interpolation
decompositionQ = [10.^([-8:1:-1]), 0.11:0.01:0.14, 0.15:0.05:1.1, 1.2:0.1:3, 3.5:0.25:13, 15:0.5:59.5, 60:2.5:100, 110:10:200, 250, 300, 500, 1000, 2000, 200000, 1000000, 10^20]; 
decompositionIncome = [10.^([-9,-4:0.5:-2]), 0.03, 0.06:0.01:2, 2.2:0.2:2.6, 3:10, 100, 10^20];

incomeVector = max(max(incomeMatrix)) .* decompositionIncome ;
incomeMat = repmat(incomeVector,length(incomeVector),1)';

% Min and Max values for the decomposition
dwellingSizeVector = param.basic_q + decompositionQ*10;
dwellingSizeMatrix = repmat(dwellingSizeVector, length(incomeVector),1);

% Precalculation
XX = incomeMat;
YY_Q = utilitySize(dwellingSizeMatrix,incomeMat);
ZZ_Q = dwellingSizeMatrix;

param.max_U = max(max(YY_Q));
param.max_q = max(dwellingSizeVector);

% The function that will be used afterwards
solus_Q_temp = @(x,y) griddata(XX, YY_Q, ZZ_Q, x, min(y,param.max_U));

% Redefine a grid (to use griddedInterpolant)
logUtilityVect = -1:0.05:log(max(max(0.4.*incomeMatrix)));
logIncome = (-1:0.1:log(max(max(incomeMatrix.*1.02))))';
logDwellingSize = log(solus_Q_temp(exp(logIncome), exp(logUtilityVect)));
Q_interp = griddedInterpolant({logUtilityVect, logIncome}, logDwellingSize, 'linear', 'none');
solus_Q = @(income,utility) exp(Q_interp(log(utility),log(income)));

%% New dimensions to the grid (we remove the locations with coeffLand = 0)
% For speed

selectedPixels = (sum(land.coeffLand) > 0.01) & (max(incomeMatrix - transportCost) > 0);
land.coeffLand = land.coeffLand(:,selectedPixels);
gridTemp = grid;
grid.distanceCBD = grid.distanceCBD(selectedPixels);
limitHousing = limitHousing(selectedPixels);
multi_proba = multi_proba(:,selectedPixels);
transportCost = transportCost(:,selectedPixels);
transportCostRDP = transportCostRDP(:,selectedPixels);
transportCostTime = transportCostTime(:,selectedPixels);
param.minimumHousingSupply = param.minimumHousingSupply(selectedPixels);
param.housing_in = param.housing_in(selectedPixels);
param.taxUrbanEdgeMat = param.taxUrbanEdgeMat(selectedPixels);
incomeMatrix = incomeMatrix(:,selectedPixels);
amenities = amenities(:,selectedPixels);

carto2 = @(x) scatter(grid.xCoord(selectedPixels),grid.yCoord(selectedPixels),250.*0.5,x,'.');


% Transport
transTemp.generalizedCost = transportCost;
transTemp.timeCost = transportCostTime;


%% Solving the model

% Four types of housing:
% 1. formal
% 2. backyard
% 3. informal
% 4. RDP (not in the solver - exogenous, only for export)

% Useful variables for the solver
diffUtility = zeros(maxIteration,size(employmentCenters,2));
simulatedPeopleHousingTypes = zeros(maxIteration,3,size(grid.distanceCBD,2)); %3 is because we have 3 types of housing in the solver
simulatedJobs = zeros(maxIteration,3,size(employmentCenters,2));
totalSimulatedJobs = zeros(maxIteration,size(employmentCenters,2));
rentMatrix = zeros(maxIteration,3,size(grid.distanceCBD,2));
errorMaxAbs = zeros(1,maxIteration);
errorMax = zeros(1,maxIteration);
errorMean = zeros(1,maxIteration);
numberError = zeros(1,maxIteration);
error = zeros(maxIteration, size(employmentCenters,2));
canceling = 0; %canceling = 1 if we press "Cancel"

% Utility for each center: variable we will adjust in the solver
Uo = zeros(maxIteration,size(employmentCenters,2)); 

% impossible_population = 1 if we cannot reach the objective population
impossiblePopulation = false(1,size(employmentCenters,2)); 
numberImpossiblePopulation = 0;
% dummy that exits the solver if we cannot reach objective for the
% remaining centers
conditionPossible = true; 

% Definition of Uo
if option.ownInitializationSolver == 0
    Uo(1,:) = averageIncome .*0.2;  % Initially, utility is set above the expected level
else
    Uo(1,:) = Uo_init;
end
    
indexIteration = 1;
if option.polycentric == 0
    convergenceFactorInitial = 0.04; % Sets the speed of the adjustment of the utility
else 
    convergenceFactorInitial = 0.003 .* (mean(averageIncome) ./ param.incomeYearReference);
end
param.convergenceFactor = convergenceFactorInitial;

% Formal housing
[simulatedJobs(indexIteration,1,:),rentMatrix(indexIteration,1,:),simulatedPeopleHousingTypes(indexIteration,1,:),simulatedPeople(1,:,:),housingSupply(1,:),dwellingSize(1,:), R_mat(1,:,:)]...
    = ComputeNEDUMOutput(Uo(indexIteration,:),param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParam,interestRate,incomeMatrix,multi_proba,transportCost,transportCostRDP,...
    land.coeffLand(1,:), poly, amenities, solus_R, solus_Q, functionName, 'formal');

% Backyard housing
[simulatedJobs(indexIteration,2,:),rentMatrix(indexIteration,2,:),simulatedPeopleHousingTypes(indexIteration,2,:),simulatedPeople(2,:,:),housingSupply(2,:),dwellingSize(2,:), R_mat(2,:,:)]...
    = ComputeNEDUMOutput(Uo(indexIteration,:),param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParam,interestRate,incomeMatrix,multi_proba,transportCost,transportCostRDP, ...
    land.coeffLand(2,:), poly, amenities, solus_R, solus_Q, functionName, 'backyard');

% Informal settlements
[simulatedJobs(indexIteration,3,:),rentMatrix(indexIteration,3,:),simulatedPeopleHousingTypes(indexIteration,3,:),simulatedPeople(3,:,:),housingSupply(3,:),dwellingSize(3,:), R_mat(3,:,:)]...
    = ComputeNEDUMOutput(Uo(indexIteration,:),param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParam,interestRate,incomeMatrix,multi_proba,transportCost,transportCostRDP,...
    land.coeffLand(3,:), poly, amenities, solus_R, solus_Q, functionName, 'informal');

% Total simulated population
totalSimulatedJobs(indexIteration,:) = sum(simulatedJobs(indexIteration,:,:),2);

% deriv_U will be used to adjust the utility levels
diffUtility(indexIteration,:) = log((totalSimulatedJobs(indexIteration,:)+10)./(employmentCenters(1,:)+10));
diffUtility(indexIteration,:) = diffUtility(indexIteration,:) .* param.convergenceFactor;
diffUtility(indexIteration,diffUtility(indexIteration,:)>0) = diffUtility(indexIteration,diffUtility(indexIteration,:)>0).*1.1;

% Difference with reality
error(indexIteration,:) = (totalSimulatedJobs(indexIteration,:)./numberJobsCenter-1)*100;
errorMaxAbs(indexIteration) = max(abs(totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./numberJobsCenter(employmentCenters(1,:)~=0)-1));
errorMax(indexIteration) = -1;
errorMean(indexIteration) = mean(abs(totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./(numberJobsCenter(employmentCenters(1,:)~=0)+0.001)-1));
numberError(indexIteration) = sum(abs(totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./numberJobsCenter(employmentCenters(1,:)~=0)-1)>precision);

% Memory
indexMemory = indexIteration;
simulatedPeopleMemory = simulatedPeople;
housingStockMemory = housingSupply;
dwellingSizeMemory = dwellingSize;
errorMeanMemory = numberError(indexMemory);

% Display of the progression bar
switch printAdvanceBar
    case 1
        F = findall(0,'type','figure','tag','TMWWaitbar');
        delete(F);
        h = waitbar(0,'simulation',...
            'CreateCancelBtn',...
            'setappdata(gcbf,''canceling'',1)');
        setappdata(h,'canceling',0)
        frequencyPrint=5;
    case 0
        frequencyPrint=20;
end


while (indexIteration < maxIteration) && (errorMaxAbs(indexIteration) > precision) && conditionPossible 
    
    switch printAdvanceBar
        case 1
            if getappdata(h,'canceling')
                canceling=1;
                break
            end
    end
    
    % Iteration
    indexIteration = indexIteration + 1;
    
    % Adjusting the level of utility
    Uo(indexIteration,:) = exp(log(Uo(indexIteration-1,:)) + diffUtility(indexIteration-1,:)); 
    
    % Minimum and maximum levels of utility
    Uo(indexIteration, Uo(indexIteration,:)<0) = 10;
    Uo(indexIteration, impossiblePopulation) = 10; % For the centers for which the objective cannot be attained (impossible_population = 1), utility level is set at an arbitrary low level
    Uo(indexIteration, employmentCenters(1,:)==0) = 10000000;
    
    % Adjusting param.factor_convergence
    param.convergenceFactor = convergenceFactorInitial ./ (1 + 1.*abs((totalSimulatedJobs(indexIteration,:)+100) ./ (numberJobsCenter+100)-1)); %.*(Jval./mean(Jval)).^0.3 %We adjust the parameter to how close we are from objective 
    param.convergenceFactor = param.convergenceFactor .* (1 - 0.6.*indexIteration./maxIteration);
        
    
    % Formal housing
    [simulatedJobs(indexIteration,1,:),rentMatrix(indexIteration,1,:),simulatedPeopleHousingTypes(indexIteration,1,:),simulatedPeople(1,:,:),housingSupply(1,:),dwellingSize(1,:), R_mat(1,:,:)]...
        = ComputeNEDUMOutput(Uo(indexIteration,:),param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParam,interestRate,incomeMatrix,multi_proba,transportCost,transportCostRDP,...
        land.coeffLand(1,:), poly, amenities, solus_R, solus_Q, functionName, 'formal');

    % Backyard housing
    [simulatedJobs(indexIteration,2,:),rentMatrix(indexIteration,2,:),simulatedPeopleHousingTypes(indexIteration,2,:),simulatedPeople(2,:,:),housingSupply(2,:),dwellingSize(2,:), R_mat(2,:,:)]...
        = ComputeNEDUMOutput(Uo(indexIteration,:),param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParam,interestRate,incomeMatrix,multi_proba,transportCost,transportCostRDP,...
        land.coeffLand(2,:), poly, amenities, solus_R, solus_Q, functionName, 'backyard');

    % Informal settlements
    [simulatedJobs(indexIteration,3,:),rentMatrix(indexIteration,3,:),simulatedPeopleHousingTypes(indexIteration,3,:),simulatedPeople(3,:,:),housingSupply(3,:),dwellingSize(3,:), R_mat(3,:,:)]...
        = ComputeNEDUMOutput(Uo(indexIteration,:),param,option,transTemp,grid,agriculturalRent,limitHousing,referenceRent,constructionParam,interestRate,incomeMatrix,multi_proba,transportCost,transportCostRDP,...
        land.coeffLand(3,:), poly, amenities, solus_R, solus_Q, functionName, 'informal');

    % Total simulated population
    totalSimulatedJobs(indexIteration,:) = sum(simulatedJobs(indexIteration,:,:),2);
    
        
    % deriv_U will be used to adjust the utility levels
    diffUtility(indexIteration,:) = log((totalSimulatedJobs(indexIteration,:)+10)./(employmentCenters(1,:)+10));
    diffUtility(indexIteration,:) = diffUtility(indexIteration,:) .* param.convergenceFactor;
    diffUtility(indexIteration,diffUtility(indexIteration,:)>0) = diffUtility(indexIteration,diffUtility(indexIteration,:)>0).*1.1;

    % Variables to display
    error(indexIteration,:) = (totalSimulatedJobs(indexIteration,:)./numberJobsCenter-1)*100;
    [errorMaxAbs(indexIteration),m] = max(abs(totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./numberJobsCenter(employmentCenters(1,:)~=0) - 1));
    erreur_temp = (totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./numberJobsCenter(employmentCenters(1,:)~=0) - 1);
    errorMax(indexIteration) = erreur_temp(m);
    errorMean(indexIteration) = mean(abs(totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./(numberJobsCenter(employmentCenters(1,:)~=0)+0.001) - 1));
    numberError(indexIteration) = sum(abs(totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./numberJobsCenter(employmentCenters(1,:)~=0) - 1)>precision);
    
    
    % Advance of the progression bar
    if floor(indexIteration/frequencyPrint) == indexIteration/frequencyPrint
        textPrinted = sprintf('Erreur max: %1.2f %%, moyenne: %1.2f %%, %g au dessus',...
            errorMax(indexIteration)*100, errorMean(indexIteration)*100, numberError(indexIteration)); 
        switch printAdvanceBar
            case 1
                waitbar(indexIteration/maxIteration,h, textPrinted);
            case 0
                disp(textPrinted);
        end
    end   
    
    %In case, for one type of households, it is impossible to attain the
    %objective population (basic need effect)
    if ((sum(Uo(indexIteration,:)<1)>0) && (max((totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./numberJobsCenter(employmentCenters(1,:)~=0)-1))<precision))
        impossiblePopulation(Uo(indexIteration,:)<1) = true;
    end
    if (sum(impossiblePopulation) + sum(abs(totalSimulatedJobs(indexIteration,employmentCenters(1,:)~=0)./numberJobsCenter(employmentCenters(1,:)~=0)-1)<precision))>=length(poly.incomeMult) %If we have to stop the solver
        if sum(impossiblePopulation)==numberImpossiblePopulation
            conditionPossible = false; %We exit the solver
        else 
           numberImpossiblePopulation = sum(impossiblePopulation); %Gives the centers for which the model could not solve
        end
    end
    impossiblePopulation(totalSimulatedJobs(indexIteration,:)>(1+precision).*numberJobsCenter) = 0; %In case there are problems with initialization
    
    
    % The best solution attained is stored in memory
    if numberError(indexIteration)<=errorMeanMemory
        indexMemory = indexIteration;
        simulatedPeopleMemory = simulatedPeople;
        errorMeanMemory = numberError(indexMemory); 
        housingStockMemory = housingSupply;
        dwellingSizeMemory = dwellingSize;
    end
end

switch printAdvanceBar
    case 1
        F = findall(0,'type','figure','tag','TMWWaitbar');
        delete(F);
end

% Different exit options for the solver
if indexIteration==maxIteration
    disp('nombre max d''iterations atteint');
    indexIteration = indexMemory;
    simulatedPeople = simulatedPeopleMemory;
    housingSupply = housingStockMemory;
    dwellingSize = dwellingSizeMemory;
    canceling = 1; % So that it displays graphs
elseif canceling
    disp('simulation arretee avant la fin');
    indexIteration = indexMemory;
    simulatedPeople = simulatedPeopleMemory;
    housingSupply = housingStockMemory;
    dwellingSize = dwellingSizeMemory;
    toc
end

textPrinted=sprintf('Max error: %1.2f %%, mean: %1.2f %%, %g above tolerance',...
    errorMaxAbs(indexIteration)*100,errorMean(indexIteration)*100,numberError(indexIteration));
disp(textPrinted);

toc

%%%%%%%%%%%%%%%%%%%%%%%
% RDP houses 
%%%%%%%%%%%%%%%%%%%%%%%

householdsRDP = land.numberPropertiesRDP .* totalRDP ./ sum(land.numberPropertiesRDP);
constructionRDP = repmat(param.sizeRDP ./ (param.sizeRDP + param.sizeBackyard), 1, length(gridTemp.distanceCBD)) .* 1000000;
dwellingSizeRDP = repmat(param.sizeRDP, 1, length(gridTemp.distanceCBD));

simulatedPeopleWithRDP = zeros(4, length(poly.incomeMult), length(gridTemp.distanceCBD));
simulatedPeopleWithRDP(1:3,:,selectedPixels) = simulatedPeople;
simulatedPeopleWithRDP(4, poly.incomeGroup == 1, :) = repmat(householdsRDP, sum(poly.incomeGroup == 1), 1).*numberJobsCenter(poly.incomeGroup == 1)'./sum(numberJobsCenter(poly.incomeGroup == 1));

%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs of the solver 
%%%%%%%%%%%%%%%%%%%%%%%%

% Employment centers
initialState.error = error(indexIteration,:);
initialState.simulatedJobs(:,:) = simulatedJobs(indexIteration,:,:);

% Number of people
initialState.householdsHousingType(:,:) = sum(simulatedPeopleWithRDP,2);
initialState.householdsCenter(:,:) = sum(simulatedPeopleWithRDP,1);
initialState.households = simulatedPeopleWithRDP;

% Housing stock and dwelling size
housingSupplyExport = zeros(3, length(gridTemp.distanceCBD));
dwellingSizeExport = zeros(3, length(gridTemp.distanceCBD));
housingSupplyExport(:,selectedPixels) = housingSupply;
dwellingSizeExport(:,selectedPixels) = dwellingSize;
dwellingSizeExport(dwellingSize<=0) = NaN;
initialState.dwellingSize = [dwellingSizeExport; dwellingSizeRDP];
initialState.housingSupply = [housingSupplyExport; constructionRDP]; 

% Rents (hh in RDP pay a rent of 0)
rentTemp(:,:) = rentMatrix(indexIteration, :, :);
rentExport = zeros(3,length(gridTemp.distanceCBD));
rentExport(:,selectedPixels) = rentTemp;
rentExport(:,selectedPixels==0) = NaN;
initialState.rent(:,:) = [rentExport; zeros(1,length(gridTemp.distanceCBD))];
rentMatrixExport = zeros(3,length(poly.xCenter),length(gridTemp.distanceCBD));
rentMatrixExport(:,:,selectedPixels) = R_mat;
rentMatrixExport(:,:,selectedPixels == 0) = NaN;
initialState.rentMatrix = rentMatrixExport;

% Other outputs
initialState.capitalLand = (housingSupply./(param.coeff_grandA)).^(1/param.coeff_b);
initialState.incomeMatrix = incomeMatrix;
initialState.limitCity = (initialState.households > 1);
initialState.utility = Uo(indexIteration,:);
initialState.impossiblePopulation = impossiblePopulation;

%% Drawing control graphs

if canceling
    DrawGraphSolver(incomeMatrix, indexIteration, error)
end

end


function [ utility ] = ComputeUtilityFromRent(Ro, income, basic_q, param)
% Precalculation of the utility
%   [ utili ] = utilite( Ro,revenu,basic_q,param ) estimates
%   utility for a rent "Ro" and an income (net of transport costs)
%   "revenu". 
% Without amenities 

if (basic_q~=0)
      utility = param.alpha^param.alpha*param.beta^param.beta...
            .*sign(income - basic_q.*Ro)...
            .*abs(income - basic_q.*Ro)...
            ./(Ro.^param.beta);
    utility((income - basic_q.*Ro) < 0) = 0;
else
    utility = param.alpha^param.alpha*param.beta^param.beta...
        .*income./(Ro.^param.beta);
end

% If the income is 0 then utility is 0 (we need to specify it to avoid
% NaNs later)
utility(income==0) = 0;


end

function [ utility ] = ComputeUtilityFromDwellingSize(q, income, basic_q, param)
% Precalculation of the utility
%   [ utili ] = utilite( Ro,revenu,basic_q,param ) estimates
%   utility for a dwelling size "q" and an income (net of transport costs)
%   "revenu". 
% Without amenities 

if (basic_q~=0)
    utility = (param.alpha .* income).^param.alpha...
            .*(q - basic_q)...
            ./ (q - param.alpha.*basic_q).^param.alpha;
    utility(q < basic_q) = 0;
else
    utility = (param.alpha .* income).^param.alpha .* q .^ param.beta;
end

% If the income is 0 then utility is 0 (we need to specify it to avoid
% NaNs later)
utility(income == 0) = 0;


end

function DrawGraphSolver(income, indexIteration, error)
% Displays the evolution of errors during the solver computation

close all

cd = colormap('jet'); 
cd = interp1(linspace(min(log(income(:,1))),max(log(income(:,1))),length(cd)),cd,log(income(:,1)')); % map color to y values 
for i = 1:length(income(:,1))
     plot(error(1:indexIteration,i), 'Color', cd(i,:))
     hold on 
end


toc()
end


