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




function trans = LoadTransportCostCapeTown(option, grid, macro, param, poly, yearTrafic, extrapolate)
% Computes travel times and costs


global path_nedum


% nonLinearCostTime = 1 if we want non linear costs of time
nonLinearCostTime = 0;

tic
disp ('importing commuting data...');

referencement = poly.increment;

% Addendum to time spent for each mode (0 for now, can be estimated if needed)
supplementaryTimeCar = 0; 
supplementaryTimeWalking = 0;
supplementaryTimePubTransit = 0;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading transport time between TZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%**** EXPLANATION ****
%                       pixel of departure (grid)
%                           <----> [index1]
%                   ^   xxxxxxxxxxxxxxx
%                   |   xxxxxxxxxxxxxxx
%TZ of arrival      |   xxxxxxxxxxxxxxx
%                   |   xxxxxxxxxxxxxxx
%           [index2]|   xxxxxxxxxxxxxxx
%                   v   xxxxxxxxxxxxxxx
%                       xxxxxxxxxxxxxxx
%
%   Same format trans.generalized_cost...s


%lance l'algo de calcul des dur?e en metro
param.trainWaitingTime = 10;
param.walkingSpeed = 4; %5 km/h average walking speed
trans.numberModes = 5; 

% Number of distinct employment center
listCenter = unique(poly.corresp, 'stable');
coordinatesCenter = [unique(poly.xCenter, 'stable')', unique(poly.yCenter, 'stable')'];

% Old function for train and car (data from google maps and OSM). 
[distanceTrain, ~] = ImportSimulatedTransportTimeByTrain(poly,grid,param);
%[distance_car, duration_car]=import_donnees_voiture(poly, grille);

% Import transport data
load transportMatrices

if option.loadTransportTime == 0
    load timeCenters
else
    distanceCar = distance_vol_oiseau;
    durationCar = cars;
    durationTrain = train;
    distanceTrain = distance_vol_oiseau;
    durationMinibus = taxi;
    durationBus = bus;
    
    distanceCar_bak = distanceCar;
    durationCar_bak = durationCar;
    durationTrain_bak = durationTrain;
    durationTrain_bak(:,X > -8000) = NaN; %some weird data point that I remove
    distanceTrain_bak = distanceTrain;
    durationMinibus_bak = durationMinibus;
    durationBus_bak = durationBus;
    
    % Interpolation to the grid
    distanceCar_temp = zeros(length(listCenter),size(grid.distanceCBD,2));
    durationCar_temp = zeros(length(listCenter),size(grid.distanceCBD,2));
    distanceTrain_temp = zeros(length(listCenter),size(grid.distanceCBD,2));
    durationTrain_temp = zeros(length(listCenter),size(grid.distanceCBD,2));
    durationMinibus_temp = zeros(length(listCenter),size(grid.distanceCBD,2));
    durationBus_temp = zeros(length(listCenter),size(grid.distanceCBD,2));
    distanceCar = zeros(length(poly.corresp),size(grid.distanceCBD,2));
    durationCar = zeros(length(poly.corresp),size(grid.distanceCBD,2));
    distanceTrain = zeros(length(poly.corresp),size(grid.distanceCBD,2));
    durationTrain = zeros(length(poly.corresp),size(grid.distanceCBD,2));
    durationMinibus = zeros(length(poly.corresp),size(grid.distanceCBD,2));
    durationBus = zeros(length(poly.corresp),size(grid.distanceCBD,2));


    for i = 1:length(listCenter)
        whichModes = unique(poly.corresp(poly.corresp == listCenter(i)));
        distanceCar_temp(i,:) = griddata_extra(X/1000,Y/1000,distanceCar_bak(whichModes,:)',grid.xCoord,grid.yCoord, extrapolate, coordinatesCenter(i,:))';
        durationCar_temp(i,:) = griddata_extra(X/1000,Y/1000,durationCar_bak(whichModes,:)',grid.xCoord,grid.yCoord, extrapolate, coordinatesCenter(i,:))';
        distanceTrain_temp(i,:) = griddata_extra(X/1000,Y/1000,distanceTrain_bak(whichModes,:)',grid.xCoord,grid.yCoord, extrapolate, coordinatesCenter(i,:))';
        durationTrain_temp(i,:) = griddata_extra(X/1000,Y/1000,durationTrain_bak(whichModes,:)',grid.xCoord,grid.yCoord, extrapolate, coordinatesCenter(i,:))';
        durationMinibus_temp(i,:) = griddata_extra(X/1000,Y/1000,durationMinibus_bak(whichModes,:)',grid.xCoord,grid.yCoord, extrapolate, coordinatesCenter(i,:))';
        durationBus_temp(i,:) = griddata_extra(X/1000,Y/1000,durationBus_bak(whichModes,:)',grid.xCoord,grid.yCoord, extrapolate, coordinatesCenter(i,:))';
        % Correspondance
        distanceCar(poly.corresp == listCenter(i),:) = repmat(distanceCar_temp(i,:),sum(poly.corresp == listCenter(i)),1);
        durationCar(poly.corresp == listCenter(i),:) = repmat(durationCar_temp(i,:),sum(poly.corresp == listCenter(i)),1);
        distanceTrain(poly.corresp == listCenter(i),:) = repmat(distanceTrain_temp(i,:),sum(poly.corresp == listCenter(i)),1);
        durationTrain(poly.corresp == listCenter(i),:) = repmat(durationTrain_temp(i,:),sum(poly.corresp == listCenter(i)),1);
        durationMinibus(poly.corresp == listCenter(i),:) = repmat(durationMinibus_temp(i,:),sum(poly.corresp == listCenter(i)),1);
        durationBus(poly.corresp == listCenter(i),:) = repmat(durationBus_temp(i,:),sum(poly.corresp == listCenter(i)),1);
    end
    
    % trans.reliable gives us the data points where we have all the data
    trans.reliable = ones(1,length(grid.distanceCBD));
    for i=1:length(poly.xCenter)
        trans.reliable(isnan(durationCar(i,:))) = 0;
        trans.reliable(isnan(durationTrain(i,:))) = 0;
        trans.reliable(isnan(durationMinibus(i,:))) = 0;
        trans.reliable(isnan(durationBus(i,:))) = 0;
    end
    
    
    % Saving the transport times (saves time for future computation)
    save('./0. Precalculated inputs/timeCenters','distanceCar','durationCar','distanceTrain','durationTrain','durationMinibus', 'durationBus')
    
end

% Sinuosity parameter
LengthPrivateCar = distanceCar;

% For transport cost in metro, we use the computed distance
LengthInVehiculePubTransit = distanceTrain;


%%  Price for public transportation

% Inputs from own analysis (from Roux 2013)
priceTrainPerKMMonth = 1.5 / 40 .* ppval(macro.splineInflation, 2012 - param.yearBegin) / ppval(macro.splineInflation, 2015-param.yearBegin);  % 0.164
priceTrainFixedMonth = 121.98 .* ppval(macro.splineInflation, 2012 - param.yearBegin) / ppval(macro.splineInflation, 2015-param.yearBegin); %4.48*40; 
priceTaxiPerKMMonth = 0.785;
priceTaxiFixedMonth = 4.32*40;
priceBusPerKMMonth = 0.522; 
priceBusFixedMonth = 6.24*40;

% Inputs from Claus
% priceTrainPerKMMonth = 0.94;
% priceTrainFixedMonth = 0;
% priceTaxiPerKMMonth = 0.95;
% priceTaxiFixedMonth = 0;
% priceBusPerKMMonth = prix_taxi_2012_km; 
% priceBusFixedMonth = prix_taxi_2012_fixe_mois;


% Correct for inflation
inflation = ppval(macro.splineInflation,yearTrafic);
infla_2012 = ppval(macro.splineInflation, 2012 - param.yearBegin);
priceTrainPerKMMonth = priceTrainPerKMMonth.*inflation ./ infla_2012;
priceTrainFixedMonth = priceTrainFixedMonth .* inflation ./ infla_2012;
priceTaxiPerKMMonth = priceTaxiPerKMMonth .* inflation ./ infla_2012;
priceTaxiFixedMonth = priceTaxiFixedMonth .* inflation ./ infla_2012;
priceBusPerKMMonth = priceBusPerKMMonth .* inflation ./ infla_2012;
priceBusFixedMonth = priceBusFixedMonth .* inflation ./ infla_2012;


% Price from private car, from own analysis
priceFixedVehiculeMonth = 300;
priceFixedVehiculeMonth = priceFixedVehiculeMonth .* inflation./infla_2012;
priceFuel = ppval(macro.splineFuelCost, yearTrafic);

% From Claus
% priceFixedVehiculeMonth = 0 .* inflation./infla_2012;
% priceFuel = 4.85 .* inflation./infla_2012;


% Transport times
timePV = durationCar;
timeTrain = durationTrain; % duration_metro_2;
timeTaxi = durationMinibus;
timeBus = durationBus;

% For each year, we esimate the price per km for cars
priceFuelPerKMMonth = zeros(size(priceFuel));
for index=1:length(yearTrafic)
        
    % Determining taxes accross time
    if ~(param.Tax==0)
        %param.taxe est la taxe par km
        taxAccrossTime = ((param.year_begin+yearTrafic(index))>param.annee_taxe).*param.taxe;%.*ppval(macro.spline_revenu,t_trafic(index))/ppval(macro.spline_revenu,annee_taxe-param.year_begin);
    else
        taxAccrossTime = 0;
    end
    
    % Adding taxes, and multiplying by 2 (round-trip) and 20 (days per month) 
    priceFuelPerKMMonth(index) = (priceFuel(index) + taxAccrossTime)*2*20;
    
end

% Time by each mode
timeWalkingTemp = LengthPrivateCar/param.walkingSpeed*60 + supplementaryTimeWalking;
timeWalkingTemp(isnan(timePV)) = NaN;%si on ne fait pas ?a, on a des 0 au lieu d'avoir des nan
timeOutput(:,:,1) = timeWalkingTemp;
timeOutput(:,:,2) = timeTrain + supplementaryTimePubTransit;
timeOutput(:,:,3) = timePV + supplementaryTimeCar;
timeOutput(:,:,4) = timeTaxi + supplementaryTimePubTransit;
timeOutput(:,:,5) = timeBus + supplementaryTimePubTransit;
timeOutput = single(timeOutput);

% Monetary cost 
multiplierPrice(:,:,1) = zeros(size(timeOutput(:,:,1)));
multiplierPrice(:,:,2) = LengthInVehiculePubTransit;
multiplierPrice(:,:,3) = LengthPrivateCar;
multiplierPrice(:,:,4) = LengthPrivateCar;
multiplierPrice(:,:,5) = LengthPrivateCar;
pricePerKM(:,1) = ones(size(priceFuelPerKMMonth));
pricePerKM(:,2) = priceTrainPerKMMonth*20*2*12;
pricePerKM(:,3) = priceFuelPerKMMonth*12;
pricePerKM(:,4) = priceTaxiPerKMMonth*2*20*12;
pricePerKM(:,5) = priceBusPerKMMonth*2*20*12;


% Distances (not useful to calculate price but useful output)
distanceOutput(:,:,1) = LengthPrivateCar;
distanceOutput(:,:,2) = LengthInVehiculePubTransit;
distanceOutput(:,:,3) = LengthPrivateCar;
distanceOutput(:,:,4) = LengthPrivateCar;
distanceOutput(:,:,5) = LengthPrivateCar;

trans.distanceOutput = single(distanceOutput);

monetaryCost = zeros(length(poly.codeCentersPolycentric(referencement)),size(timeOutput,2),trans.numberModes);
generalizedCost = single(zeros(length(poly.codeCentersPolycentric(referencement)),size(timeOutput,2),length(yearTrafic)));
whichModes = uint8(zeros(size(poly.codeCentersPolycentric(referencement),1),size(timeOutput,2)));


% Household size (actually the number of commuters per household)
householdSizeMat = repmat(param.householdSizeTransport, 1, length(poly.codeCentersInitial) ./ param.numberIncomeGroup);
householdSizeMat = householdSizeMat(poly.selectedCenters)' * ones(1, length(grid.distanceCBD)); 

% For each year calculate the generalized cost
for index = 1:length(yearTrafic)
    for index2 = 1:trans.numberModes
        monetaryCost(:,:,index2) = pricePerKM(index,index2).*multiplierPrice(:,:,index2);
        monetaryCost(:,:,index2) = monetaryCost(:,:,index2).*householdSizeMat;
    end
    
    incomeCommuters = InterpolateIncomeEvolution(macro, param, option, grid, poly, yearTrafic(index));
    incomeCommuters = repmat(incomeCommuters,[1 1 size(monetaryCost,3)]);
    clear('transi.generalized_cost'); 
    %ajout des couts fixes
    monetaryCost(:,:,2) = monetaryCost(:,:,2) + priceTrainFixedMonth(index)*12.*householdSizeMat; % train (monthly fare)
    monetaryCost(:,:,3) = monetaryCost(:,:,3) + priceFixedVehiculeMonth(index)*12*householdSizeMat; % private car
    monetaryCost(:,:,4) = monetaryCost(:,:,4) + priceTaxiFixedMonth(index)*12.*householdSizeMat; % minibus-taxi
    monetaryCost(:,:,5) = monetaryCost(:,:,5) + priceBusFixedMonth(index)*12.*householdSizeMat; % bus

    numberHourWorkedPerWeek = 40;
    numberWeeksPerYear = 52;
    incomePerHour = incomeCommuters/numberWeeksPerYear/numberHourWorkedPerWeek;
    costTime = timeOutput .* param.timeCost .* incomePerHour/60*2*20*12;
    if nonLinearCostTime == 1
        costTime(timeOutput > param.limite_temps) = (param.timeLimit.*param.timeCost + (timeOutput(timeOutput > param.timeLimit) - param.timeLimit).*param.timeCost2)...
            .*incomePerHour(timeOutput > param.timeLimit)/60*2*20*12;
    end
    finalCost = monetaryCost + costTime;
        
    if index==1
        trans.monetaryCostInitial = monetaryCost;
        trans.timeCostInitial = costTime;
    end
    
    if option.logit == 1
        minimumCost(:,:,1) = min(finalCost,[],3);
        minimumCost(:,:,2) = minimumCost(:,:,1);
        minimumCost(:,:,3) = minimumCost(:,:,1);
        minimumCost(:,:,4) = minimumCost(:,:,1);
        minimumCost(:,:,5) = minimumCost(:,:,1);

        coeffLogit = param.logitFactor./minimumCost;
        modeLogit = logit(coeffLogit, finalCost, trans.numberModes);
        
        generalizedCost(:,:,index) = single(ComputeLogitMean(coeffLogit, finalCost)./param.logitFactor);
        generalizedCost(:,:,index) = single(ComputeLogitMean(coeffLogit, finalCost)./coeffLogit(:,:,1));
                
        whichModes = single(modeLogit);
        trans.timeCost(:,:,index) = sum(whichModes.*costTime,3);
        
    else
        
        [generalizedCost(:,:,index),whichModes(:,:,index)] = min(single(finalCost),[],trans.numberModes);
        
    end
    fprintf('%d %%, year %d\n',round((index)/length(yearTrafic)*100),yearTrafic(index) + param.yearBegin);
end

% Output
trans.generalizedCost = generalizedCost;
disp('Estimation of accessibility done');


% If we want a general accessibility index
% (useful for calibration of amenities)    
disp('Also computing a general accessibility index');
destinationWeightGridTemp = zeros(length(listCenter), length(grid.distanceCBD)); 
destinationWeightGrid = zeros(length(poly.corresp), length(grid.distanceCBD)); 
for i = 1:length(listCenter)
    destinationWeightGridTemp(i,:) = griddata_extra(X/1000,Y/1000, poly.destinationWeight(i,:)', grid.xCoord, grid.yCoord, 0);
    destinationWeightGrid(poly.corresp == listCenter(i),:) = repmat(destinationWeightGridTemp(i,:), sum(poly.corresp == listCenter(i)), 1);
end
% One value per income
generalizedCostIndex = zeros(param.numberIncomeGroup, length(grid.distanceCBD));
for i = 1:param.numberIncomeGroup
    generalizedCostIndex(i,:) = nansum(generalizedCost(poly.incomeGroup(1,:) == i,:,1) .* destinationWeightGrid(poly.incomeGroup(1,:) == i,:)) ./ ...
                                  nansum(destinationWeightGrid(poly.incomeGroup(1,:) == i,:));
end
generalizedCostIndex(:,trans.reliable == 0) = NaN;
trans.generalizedCostIndex = generalizedCostIndex;

% Other outputs
trans.yearTransport = yearTrafic + param.yearBegin;
trans.timeOutput = timeOutput;
trans.whichModes = whichModes;

toc()

end



function sortie = ComputeLogitMean(coeff_logit,price)
%coeff c'est le facteur dans le calcul et prix est le vecteur des prix

AAA = exp(-coeff_logit.*price);
sortie = -log(sum(AAA,3)) + 0.5772;

end


function sortie = logit(coeff_logit,cost, numberMode)%coeff c'est le facteur dans le calcul et prix est le vecteur des prix

AAA = exp(- coeff_logit.*cost);
BBB = sum(AAA,3);

CCC = zeros(size(AAA));
for i = 1:numberMode
    CCC(:,:,i)=BBB;
end

sortie = AAA./CCC;

end



function output = griddata_extra(a,b,c,x,y, extrapolate, coordinate_center)
% Like a griddata but with extrapolation outside the area of data
% availability

test = ~isnan(c);

% scatteredInterpolant does not extrapolate if 'none' is the second method
% parameter

if extrapolate == 1

    x_center = coordinate_center(1);
    y_center = coordinate_center(2);

    % 'linear' is a linear extrapolation based on the gradient at the border
    % (can cause problem locally - decreasing transport times)
    surface_linear = scatteredInterpolant(a(test), b(test), c(test), 'linear', 'none');
    
    % Extrapolation by linear trend on angle slices
    interpolated_values = surface_linear(x',y');
    extrapolated_values = zeros(length(interpolated_values), 1);
    interv_angle = pi/20;
    angles = [interv_angle/2:interv_angle:(2*pi - interv_angle/2)];
    grid_angles = atan2(y - y_center, x - x_center) + pi;
    grid_distance = sqrt((x - x_center).^2 + (y - y_center).^2);
    
    for theta = 1:length(angles)
         which_angles = (grid_angles >= angles(theta) - interv_angle / 2) & (grid_angles <= angles(theta) + interv_angle / 2);
         coefficient = nanmean(interpolated_values(which_angles) ./ grid_distance(which_angles)');
         [max_value, which_max] = max(interpolated_values(which_angles' & ~isnan(interpolated_values)));
         max_distance = grid_distance(which_angles & ~isnan(interpolated_values'));
         max_distance = max_distance(which_max);
         extrapolated_values(which_angles) = max_value + coefficient .* (grid_distance(which_angles) - max_distance)';
    end
    
    output = interpolated_values;
    output(isnan(output)) = extrapolated_values(isnan(output));
    
else
    
    surface_linear = scatteredInterpolant(a(test), b(test), c(test), 'linear', 'none');
    output = surface_linear(x', y');

end

end