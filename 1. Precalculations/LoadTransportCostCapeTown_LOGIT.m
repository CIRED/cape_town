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




function trans = LoadTransportCostCapeTown_LOGIT(option, grid, macro, param, poly, data, yearTraffic, areaInterp, extrapolate)
% Computes travel times and costs
% FOR A LOGIT MODEL OF ACCESSIBILITY


tic
disp ('importing commuting data...');

% Addendum to time spent for each mode (0 for now, can be estimated if needed)
supplementaryTimeCar = 0; 
supplementaryTimeWalking = 0;
supplementaryTimePubTransit = 0;

% Parameters
param.walkingSpeed = 4; % in km/h average walking speed
trans.numberModes = 5; 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loading transport time between TZ

%**** EXPLANATION ****
%                       Area of departure
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

% Number of distinct employment center
listCenter = unique(poly.corresp, 'stable');
coordinatesCenter = [unique(poly.xCenter, 'stable'), unique(poly.yCenter, 'stable')];

% Old function for train and car (data from google maps and OSM). 
% param.trainWaitingTime = 0;
% [distanceTrain, ~] = ImportSimulatedTransportTimeByTrain(poly,grid,param);
% [distance_car, duration_car]=import_donnees_voiture(poly, grille);

% Import transport data
load transportMatrices

distanceCar = distance_vol_oiseau;
durationCar = cars;
durationTrain = train;
durationTrain(:,X > -8000) = NaN; %some weird data point that I remove
durationMinibus = taxi;
durationBus = bus;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Keeping only interesting employment centers

% Selected employment centers
distanceCar_bak = distanceCar(poly.selectedCenters, :);
durationCar_bak = durationCar(poly.selectedCenters, :);
durationTrain_bak = durationTrain(poly.selectedCenters, :);
durationMinibus_bak = durationMinibus(poly.selectedCenters, :);
durationBus_bak = durationBus(poly.selectedCenters, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Interpolating at the desired level

if areaInterp == 'GR' % grid
    xInterp = grid.xCoord;
    yInterp = grid.yCoord;
    loadAdress = './0. Precalculated Inputs/Transport_times_GRID.mat';
elseif areaInterp == 'SP'
    xInterp = data.spX;
    yInterp = data.spY;
    loadAdress = './0. Precalculated Inputs/Transport_times_SP.mat';
end

distanceCar = zeros(length(poly.corresp), length(xInterp));
durationCar = zeros(length(poly.corresp), length(xInterp));
durationTrain = zeros(length(poly.corresp), length(xInterp));
durationMinibus = zeros(length(poly.corresp), length(xInterp));
durationBus = zeros(length(poly.corresp), length(xInterp));
%%
if option.loadTransportTime == 0    
    for i = 1:length(listCenter)
        distanceCar(i,:) = griddata_extra(X/1000,Y/1000, distanceCar_bak(i,:)', xInterp, yInterp, extrapolate, coordinatesCenter(i,:))';
        durationCar(i,:) = griddata_extra(X/1000,Y/1000, durationCar_bak(i,:)', xInterp, yInterp, extrapolate, coordinatesCenter(i,:))';
        durationTrain(i,:) = griddata_extra(X/1000,Y/1000, durationTrain_bak(i,:)', xInterp, yInterp, extrapolate, coordinatesCenter(i,:))';
        durationMinibus(i,:) = griddata_extra(X/1000,Y/1000, durationMinibus_bak(i,:)', xInterp, yInterp, extrapolate, coordinatesCenter(i,:))';
        durationBus(i,:) = griddata_extra(X/1000,Y/1000, durationBus_bak(i,:)', xInterp, yInterp, extrapolate, coordinatesCenter(i,:))';
    end
    save(loadAdress, 'distanceCar', 'durationCar', 'durationTrain', 'durationMinibus', 'durationBus');
else 
    load(loadAdress)
    disp('Transport data loaded directly')
end

%% Distance in train
distanceTrain = distanceCar;

% Sinuosity parameter
LengthPrivateCar = distanceCar .* 1; % Already in the data

% For transport cost in metro, we use the computed distance
LengthInVehiculePubTransit = distanceTrain;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Price per km and fixed costs

% Inputs from own analysis (from Roux 2013)
priceTrainPerKMMonth = 0.164 .* ppval(macro.splineInflation, 2011 - param.yearBegin) / ppval(macro.splineInflation, 2013 - param.yearBegin); % 1.5 / 40 .* ppval(macro.splineInflation, 2012 - param.yearBegin) / ppval(macro.splineInflation, 2015-param.yearBegin);  
priceTrainFixedMonth = 4.48*40 .* ppval(macro.splineInflation, 2011 - param.yearBegin) / ppval(macro.splineInflation, 2013-param.yearBegin); % 121.98 .* ppval(macro.splineInflation, 2012 - param.yearBegin) / ppval(macro.splineInflation, 2015-param.yearBegin); 
priceTaxiPerKMMonth = 0.785 .* ppval(macro.splineInflation, 2011 - param.yearBegin) / ppval(macro.splineInflation, 2013-param.yearBegin);
priceTaxiFixedMonth = 4.32*40 .* ppval(macro.splineInflation, 2011 - param.yearBegin) / ppval(macro.splineInflation, 2013-param.yearBegin);
priceBusPerKMMonth = 0.522 .* ppval(macro.splineInflation, 2011 - param.yearBegin) / ppval(macro.splineInflation, 2013-param.yearBegin); 
priceBusFixedMonth = 6.24*40 .* ppval(macro.splineInflation, 2011 - param.yearBegin) / ppval(macro.splineInflation, 2013-param.yearBegin);

% Inputs from Claus
% priceTrainPerKMMonth = 0.94;
% priceTrainFixedMonth = 0;
% priceTaxiPerKMMonth = 0.95;
% priceTaxiFixedMonth = 0;
% priceBusPerKMMonth = prix_taxi_2012_km; 
% priceBusFixedMonth = prix_taxi_2012_fixe_mois;

% Correct for inflation
inflation = ppval(macro.splineInflation,yearTraffic);
infla_2012 = ppval(macro.splineInflation, 2012 - param.yearBegin);
priceTrainPerKMMonth = priceTrainPerKMMonth.*inflation ./ infla_2012;
priceTrainFixedMonth = priceTrainFixedMonth .* inflation ./ infla_2012;
priceTaxiPerKMMonth = priceTaxiPerKMMonth .* inflation ./ infla_2012;
priceTaxiFixedMonth = priceTaxiFixedMonth .* inflation ./ infla_2012;
priceBusPerKMMonth = priceBusPerKMMonth .* inflation ./ infla_2012;
priceBusFixedMonth = priceBusFixedMonth .* inflation ./ infla_2012;

priceFixedVehiculeMonth = 400; 
priceFixedVehiculeMonth = priceFixedVehiculeMonth .* inflation./infla_2012;
priceFuel = ppval(macro.splineFuelCost, yearTraffic);

% From Claus
% priceFixedVehiculeMonth = 0 .* inflation./infla_2012;
% priceFuel = 4.85 .* inflation./infla_2012;

%% Transport times
timePV = durationCar;
timeTrain = durationTrain; % duration_metro_2;
timeTaxi = durationMinibus;
timeBus = durationBus;

%% For each year, we esimate the price per km for cars
priceFuelPerKMMonth = zeros(size(priceFuel));
for index=1:length(yearTraffic)
        
    % Determining taxes accross time
    if ~(param.Tax==0)
        %param.taxe is the tax per km
        taxAccrossTime = ((param.year_begin+yearTraffic(index))>param.annee_taxe).*param.taxe;%.*ppval(macro.spline_revenu,t_trafic(index))/ppval(macro.spline_revenu,annee_taxe-param.year_begin);
    else
        taxAccrossTime = 0;
    end
    
    % Adding taxes 
    priceFuelPerKMMonth(index) = priceFuel(index) + taxAccrossTime;
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Transport times and costs as matrices

% Time by each mode
timeWalkingTemp = LengthPrivateCar/param.walkingSpeed*60.*1.2*2 + supplementaryTimeWalking; % including a sinuosity factor + round-trip
timeWalkingTemp(isnan(timePV)) = NaN; % Without this, we have 0 instead of NaN
timeOutput(:,:,1) = timeWalkingTemp;
timeOutput(:,:,2) = timeTrain + supplementaryTimePubTransit;
timeOutput(:,:,3) = timePV + supplementaryTimeCar;
timeOutput(:,:,4) = timeTaxi + supplementaryTimePubTransit;
timeOutput(:,:,5) = timeBus + supplementaryTimePubTransit;
timeOutput = single(timeOutput);

%% Monetary cost 
multiplierPrice(:,:,1) = zeros(size(timeOutput(:,:,1)));
multiplierPrice(:,:,2) = LengthInVehiculePubTransit;
multiplierPrice(:,:,3) = LengthPrivateCar;
multiplierPrice(:,:,4) = LengthPrivateCar;
multiplierPrice(:,:,5) = LengthPrivateCar;

%% Number of worked days per year
numberDaysPerYear = 235;

% Multiplying by 20 (days per month)
pricePerKM(:,1) = zeros(size(priceFuelPerKMMonth));
pricePerKM(:,2) = priceTrainPerKMMonth*numberDaysPerYear;
pricePerKM(:,3) = priceFuelPerKMMonth*numberDaysPerYear;
pricePerKM(:,4) = priceTaxiPerKMMonth*numberDaysPerYear;
pricePerKM(:,5) = priceBusPerKMMonth*numberDaysPerYear;

% Distances (not useful to calculate price but useful output)
distanceOutput(:,:,1) = LengthPrivateCar;
distanceOutput(:,:,2) = LengthInVehiculePubTransit;
distanceOutput(:,:,3) = LengthPrivateCar;
distanceOutput(:,:,4) = LengthPrivateCar;
distanceOutput(:,:,5) = LengthPrivateCar;

trans.distanceOutput = single(distanceOutput);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Household sizes

% (actually the number of commuters per household)
% householdSizeMat = repmat(param.householdSizeTransport, 1, length(poly.codeCentersInitial) ./ param.numberIncomeGroup);
% householdSizeMat = householdSizeMat(poly.selectedCenters)' * ones(1, length(grid.distanceCBD)); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Monetary cost

% For each year calculate the monetary cost
monetaryCost = zeros(length(poly.codeCentersPolycentric),size(timeOutput,2),trans.numberModes);
trans.monetaryCost = zeros(length(poly.codeCentersPolycentric),size(timeOutput,2),trans.numberModes, length(yearTraffic));
for index = 1:length(yearTraffic)
    
    for index2 = 1:trans.numberModes
        monetaryCost(:,:,index2) = pricePerKM(index,index2) .* multiplierPrice(:,:,index2);
        monetaryCost(:,:,index2) = monetaryCost(:,:,index2);  %.*householdSizeMat;
    end
    
    % Adding fixed costs
    monetaryCost(:,:,2) = monetaryCost(:,:,2) + priceTrainFixedMonth(index)*12; % train (monthly fare)
    monetaryCost(:,:,3) = monetaryCost(:,:,3) + priceFixedVehiculeMonth(index)*12; % private car
    monetaryCost(:,:,4) = monetaryCost(:,:,4) + priceTaxiFixedMonth(index)*12; % minibus-taxi
    monetaryCost(:,:,5) = monetaryCost(:,:,5) + priceBusFixedMonth(index)*12; % bus

    trans.monetaryCost(:,:,:,index) = monetaryCost;
    
end


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Cost associated with time

numberHourWorkedPerDay= 8;
incomePerHour = 1/numberDaysPerYear/numberHourWorkedPerDay;
costTime = timeOutput .* param.timeCost .* incomePerHour /60 * numberDaysPerYear;
        
trans.timeCost = costTime;
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  If areaInterp == "Grid", we compute the income net of Commuting Cost,
%  ODflows, and modalShares, used in the model

if areaInterp == 'GR'
    
    lambda = param.lambda;
    incomeNetOfCommuting = zeros(param.numberIncomeGroup, length(grid.distanceCBD), length(yearTraffic));
    averageIncome = zeros(param.numberIncomeGroup, length(grid.distanceCBD), length(yearTraffic));
    modalShares = zeros(length(poly.incomeCentersInit), length(grid.distanceCBD), trans.numberModes, param.numberIncomeGroup, length(yearTraffic));
    ODflows = zeros(length(poly.incomeCentersInit), length(grid.distanceCBD), param.numberIncomeGroup, length(yearTraffic));
    
    for index = 1:length(yearTraffic)
        incomeGroup = InterpolateIncomeEvolution(macro,param,option,grid, poly, yearTraffic(index));
        incomeGroupRef = InterpolateIncomeEvolution(macro, param, option, grid, poly, 0);
        incomeCenters = poly.incomeCentersInit .* incomeGroup(:,1)' ./ incomeGroupRef(:,1)';
        [incomeNetOfCommuting(:,:,index), modalShares(:,:,:,:,index), ODflows(:,:,:,index), averageIncome(:,:,index)] = ComputeIncomeNetOfCommuting(param, trans, grid, poly, data, lambda, incomeCenters, areaInterp, index);
        fprintf('%d %%, year %d\n',round((index)/length(yearTraffic)*100),yearTraffic(index) + param.yearBegin);
    end
    
    %trans.incomeNetOfCommuting = incomeNetOfCommuting;
    %trans.modalShares = single(modalShares); 
    %trans.ODflows = single(ODflows);
    %trans.averageIncome = single(averageIncome);
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  End of code

% Other outputs
trans.yearTransport = yearTraffic + param.yearBegin;
trans.timeOutput = single(timeOutput);

disp('Estimation of accessibility done');
toc()

end




function output = griddata_extra(a,b,c,x,y, extrapolate, coordinate_center)
% Like a griddata but with extrapolation outside the area of data
% availability

test = ~isnan(c);
if sum(test) < 10
   
    output = NaN(length(x), 1);
    disp(sprintf('problematic employment: %g', coordinate_center(1), coordinate_center(2)))

else

% scatteredInterpolant does not extrapolate if 'none' is the second method
% parameter

if extrapolate == 1

    x_center = coordinate_center(1);
    y_center = coordinate_center(2);

    % 'linear' is a linear extrapolation based on the gradient at the border
    % (can cause problem locally - decreasing transport times)
    surface_linear = scatteredInterpolant(a(test), b(test), c(test), 'linear', 'none');
    
    % Extrapolation by linear trend on angle slices
    interpolated_values = surface_linear(x, y);
    extrapolated_values = NaN(length(interpolated_values), 1);
    interv_angle = pi/10;
    angles = [interv_angle/2:interv_angle:(2*pi - interv_angle/2)];
    grid_angles = atan2(y - y_center, x - x_center) + pi;
    grid_distance = sqrt((x - x_center).^2 + (y - y_center).^2);
    
    for theta = 1:length(angles)
         which_angles = (grid_angles >= angles(theta) - interv_angle / 2) & (grid_angles <= angles(theta) + interv_angle / 2);
         if sum(which_angles .* ~isnan(interpolated_values)) ~= 0
            coefficient = nanmean(interpolated_values(which_angles) ./ grid_distance(which_angles));
            [max_value, which_max] = max(interpolated_values(which_angles & ~isnan(interpolated_values)));
            max_distance = grid_distance(which_angles & ~isnan(interpolated_values));
            max_distance = max_distance(which_max);
            extrapolated_values(which_angles) = max_value + coefficient .* (grid_distance(which_angles) - max_distance);
         end
    end
    
    output = interpolated_values;
    output(isnan(output)) = extrapolated_values(isnan(output));
    
else
    
    surface_linear = scatteredInterpolant(a(test), b(test), c(test), 'linear', 'none');
    output = surface_linear(x', y');

end

% end of if
end
end