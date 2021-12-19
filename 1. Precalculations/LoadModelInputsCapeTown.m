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




function output = LoadModelInputsCapeTown(param,option)
% Imports simulation inputs 
% Claus' inputs for future scenarios, except for price of transit by cars (R/km)

disp('import parameters...');

% Location of function
global path_nedum
global slash

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% importation du sc?nario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Methods to interpolate
method = 'linear';

% Path for the scenarios
pathScenarios = strcat(path_nedum, 'Scenarios', slash);

%%%  1. Import distribution of incomes %%% 

% Path for the scenario file
pathScenarioIncomeDistribution = strcat(pathScenarios, 'Scenario_inc_distrib_', option.scenarioIncomeDistribution, '.csv');
importfile(pathScenarioIncomeDistribution);


% Import data for 2001 and 2011
importfile([path_nedum,'Income_distribution_2011.csv']);

% Create export
macro.splinePopulationIncomeDistribution = interp1([2001, 2011, 2040]' - param.yearBegin, [Households_nb_2001, Households_nb, Households_nb_2040]', 'linear', 'pp');
averageIncome2001 = sum(Households_nb_2001 .* INC_med) ./ sum(Households_nb_2001);
averageIncome2011 = sum(Households_nb .* INC_med) ./ sum(Households_nb);
averageIncome2040 = sum(Households_nb_2040 .* INC_med_2040) ./ sum(Households_nb_2040);

%%% 2. Import other data from the scenario %%% 

% Import of data and scenarios

% Ex-files
% importfile([path_nedum,'Scenario_CAPE_TOWN_Claus_inputs.csv']);
% importfile([path_nedum,'Scenario_CoCT_2019.csv']);

% New files
pathScenarioPop = strcat(pathScenarios, 'Scenario_pop_', option.scenarioPop, '.csv');
pathScenarioInterestRate = strcat(pathScenarios, 'Scenario_interest_rate_', option.scenarioInterestRate, '.csv');
pathScenarioPriceFuel = strcat(pathScenarios, 'Scenario_price_fuel_', option.scenarioPriceFuel, '.csv');
pathScenarioInflation = strcat(pathScenarios, 'Scenario_inflation_', option.scenarioInflation, '.csv');
importfile(pathScenarioPop);
importfile(pathScenarioInterestRate);
importfile(pathScenarioPriceFuel);
importfile(pathScenarioInflation);

% Inflation - for transport costs 
macro.splineInflation = interp1(Year_infla(~isnan(inflation_base_2010)) - param.yearBegin, inflation_base_2010(~isnan(inflation_base_2010)), method, 'pp');

% Income - after 2011, income evolves as inflation (and it is the same for
% income for each income group)
% Inc_avg_infla = Inc_avg(~isnan(Inc_avg));
% Year_inc = Year(~isnan(Inc_avg));
% infla_inc = ppval(macro.splineInflation, Year_inc(Year_inc > param.yearBegin) - param.yearBegin) ./ ppval(macro.splineInflation, 0);
% Inc_avg_infla(Year_inc > param.yearBegin) = Inc_avg_infla(Year_inc == param.yearBegin) .* infla_inc;
% macro.splineIncome = interp1(Year_inc - param.yearBegin, Inc_avg_infla,  method, 'pp');
% macro.incomeYearReference = ppval(macro.splineIncome, param.yearReference - param.yearBegin);

% Income - after 2011, income evolves as inflation (and it is the same for
% income for each income group)
yearInc = Year_infla(~isnan(inflation_base_2010));
Inc_year_infla = ppval(interp1([2001;2011;2040], [averageIncome2001; averageIncome2011; averageIncome2040], method, 'pp'), yearInc);
inflaRef = ppval(macro.splineInflation, yearInc(yearInc > param.yearBegin) - param.yearBegin) ./ ppval(macro.splineInflation, 0);
Inc_year_infla(yearInc > param.yearBegin) = Inc_year_infla(yearInc == param.yearBegin) .* inflaRef;
macro.splineIncome = interp1(yearInc - param.yearBegin, Inc_year_infla,  method, 'pp');
macro.incomeYearReference = ppval(macro.splineIncome, param.yearReference - param.yearBegin);

% Some bug from matlab (??)
% Year_pop = x0xFFFD0xFFFD0xFFFDYear_pop;

% Rescale income by inflation (after 2011)
incomeDistribution = [INC_med'; INC_med'; INC_med'; INC_med_2040'];
incomeDistribution(Year_pop > 2011, :) = incomeDistribution(Year_pop > 2011, :) .* repmat(ppval(macro.splineInflation, Year_pop(Year_pop > 2011) - param.yearBegin) / ppval(macro.splineInflation, 2011 - param.yearBegin), 1, size(incomeDistribution,2));
macro.splineIncomeDistribution = interp1(Year_pop(~isnan(Year_pop)) - param.yearBegin, incomeDistribution(~isnan(Year_pop),:), method, 'pp');


% Interest rate
macro.splineInterestRate = interp1(Year_interest_rate(~isnan(real_interest_rate)) - param.yearBegin, real_interest_rate(~isnan(real_interest_rate)), method, 'pp');

% Price of fuel
macro.splineFuelCost = interp1(Year_fuel(~isnan(price_fuel)) - param.yearBegin, price_fuel(~isnan(price_fuel))/100, method, 'pp');

% Total population
macro.splinePopulation = interp1(Year_pop(~isnan(HH_total)) - param.yearBegin, HH_total(~isnan(HH_total)), method, 'pp');



%%% 3. Import the scenario for RDP/BNG  houses %%%

RDP_2011 = min(2.2666e+05, sum(formal(param.incomeDistribution == 1))); %(estimated as sum(data.gridFormal(data.countRDPfromGV > 0)))    % RDP_2011 = 320969; %227409; % Where from?
RDP_2001 = min(1.1718e+05, sum(Households_nb_2001(param.incomeDistribution == 1))); %(estimated as sum(data.gridFormal(data.countRDPfromGV > 0)))  % 262452; % Estimated by nb inc_1 - BY - settlement in 2001
if option.futureConstructionRDP == 1
    macro.splineRDP = interp1([2001 - param.yearBegin; 2011 - param.yearBegin; 2018 -  param.yearBegin; 2041 - param.yearBegin], [RDP_2001; RDP_2011; RDP_2011 + 7*5000; RDP_2011 + 7*5000 + 23*param.futureRatePublicHousing], method, 'pp');
else 
    macro.splineRDP = interp1([2001 - param.yearBegin; 2011 - param.yearBegin; 2018 - param.yearBegin; 2041 - param.yearBegin], [RDP_2001; RDP_2011; RDP_2011 + 7*5000; RDP_2011 + 7*5000], method, 'pp');
end


%%% 3. Import evolution of agricultural land %%%

agriculturalRent2100 = param.agriculturalRent2011 .* ppval(macro.splineInflation, 2100 - param.yearBegin) ./ ppval(macro.splineInflation, 2011 - param.yearBegin);
macro.splineAgriculturalRent = interp1([2001 - param.yearBegin; 2011 - param.yearBegin; 2100 - param.yearBegin], [param.agriculturalRent2001; param.agriculturalRent2011; agriculturalRent2100], method, 'pp');


%%% 4. Output %%%
output = macro;

end

