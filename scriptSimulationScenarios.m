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

clear;
clc;
close all;
tic();
disp ('*** NEDUM-Cape-Town - Polycentric Version - Formal and Informal housing ***');
disp ('*************** Scenarios 2040 with and without Urban Edge ****************');


global slash
slash = '/'; % To be changed when on Mac or Windows
% Load function addresses
pathProgram;
% Load data address
pathData;

%% Parameters and options

%%% Import Parameters %%%
param = ChoiceParameters;

%%% Set options for the simulation %%% 
option = ChoiceOptions; 

%%% Parameters for the scenario %%%
param.yearUrbanEdge = 2015; % in case option.urban_edge = 0, the year the constraint is removed
param.taxOutUrbanEdge = 10000;
param.coeffDoubleStorey = 0.02; % Random for now

%%% Years for the simulation %%%
param.yearBegin = 2011;
t = 0:1:29; % Up to 2040

%%% Import inputs for the simulation %%%
ImportInputsData_LOGIT;

% Transport data
yearTrafic = [0, 10, 20, 29]; % for speed
trans = LoadTransportCostCapeTown_LOGIT(option, grid, macro, param, poly, data, yearTrafic, 'GR', 1);
disp ('*** Data and parameters imported succesfully ***')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Memory for income distribution
scenarioIncomeMemory = option.scenarioIncomeDistribution;

% Combine scenarios
listScenarioPop = [1, 2, 3];
listScenarioIncomeDistribution = [1, 2, 3];
scenarios.definition = combvec(listScenarioPop, listScenarioIncomeDistribution);
numberScenario = size(scenarios.definition, 2);

% Outputs of the simulations
scenarios.numberBackyard_UE = zeros(1, numberScenario);
scenarios.numberSettlement_UE = zeros(1, numberScenario);
scenarios.numberFormal_UE = zeros(1, numberScenario);
scenarios.numberRDP_UE = zeros(1, numberScenario);
scenarios.priceCBD_UE = zeros(1, numberScenario);            
scenarios.urbanFootprint_UE = zeros(1, numberScenario);

% Outputs of the simulations
scenarios.numberBackyard_noUE = zeros(1, numberScenario);
scenarios.numberSettlement_noUE = zeros(1, numberScenario);
scenarios.numberFormal_noUE = zeros(1, numberScenario);
scenarios.numberRDP_noUE = zeros(1, numberScenario);
scenarios.priceCBD_noUE = zeros(1, numberScenario);            
scenarios.urbanFootprint_noUE = zeros(1, numberScenario);

% Land use w/ UE
option.urbanEdge = 1; %1 means we keep the urban edge
land_UE = ImportLandUseInputsCapeTown(grid, option, param, macro, data);
land_UE.amenities = land.amenities;

% Land w/ UE
option.urbanEdge = 0; %1 means we keep the urban edge
param.yearUrbanEdge = 2015; %in case option.urban_edge = 0, the year the constraint is removed
land_noUE = ImportLandUseInputsCapeTown(grid, option, param, macro, data);
land_noUE.amenities = land.amenities;

% Initial state
Uo_init = 100000;
initialState = RunEquilibriumSolverNEDUM_LOGIT(t(1),trans,option,land_UE,grid,macro,param,poly,Uo_init,param.housing_in);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numberScenario
   
    % Scenarios
    option.scenarioPop = num2str(scenarios.definition(1,s));
    option.scenarioIncomeDistribution = num2str(scenarios.definition(2,s));
    
    % Define macro scenarios
    macro = LoadModelInputsCapeTown(param, option);
    
    if (option.scenarioIncomeDistribution ~= scenarioIncomeMemory)
        
        % Income groups and employment centers
        poly = ImportEmploymentCentersCapeTown_LOGIT(grid, param, option, macro, data, t);
        poly.incomeCentersInit = incomeCentersKeep;

        % Transport times
        trans = LoadTransportCostCapeTown_LOGIT(option, grid, macro, param, poly, data, yearTrafic, 'GR', 1);
        
        % Save new memory
        scenarioIncomeMemory = option.scenarioIncomeDistribution;
        
    end
    
    % Run simulation w/ UE
    disp('*** beginning of evolution - urban edge ***');
    simulation_UE = RunDynamicEvolutionNEDUM_LOGIT(t,initialState,trans,grid,land_UE,poly,param,macro,option);
    statDynamics_UE = ComputeFinalStatisticsCapeTown(macro,option,trans,land_UE,grid,param,poly,simulation_UE);
    disp('*** end of evolution ***');

    % Save outputs w/ UE
    scenarios.numberBackyard_UE(s) = statDynamics_UE.totalBackyard(end);
    scenarios.numberSettlement_UE(s) = statDynamics_UE.totalSettlement(end);
    scenarios.numberFormal_UE(s) = statDynamics_UE.totalPrivate(end);
    scenarios.numberRDP_UE(s) = statDynamics_UE.totalRDP(end);
    scenarios.priceCBD_UE(s) = statDynamics_UE.avgRentFormalCBD(end);            
    scenarios.urbanFootprint_UE(s) = statDynamics_UE.urbanFootprint(end);

        
    % Run simulation w/out UE
    disp('*** beginning of evolution - no urban edge ***');
    simulation_noUE = RunDynamicEvolutionNEDUM_LOGIT(t,initialState,trans,grid,land_noUE,poly,param,macro,option);
    statDynamics_noUE = ComputeFinalStatisticsCapeTown(macro,option,trans,land_noUE,grid,param,poly,simulation_noUE);
    disp('*** end of evolution ***');

    % Save outputs w/out UE
    scenarios.numberBackyard_noUE(s) = statDynamics_noUE.totalBackyard(end);
    scenarios.numberSettlement_noUE(s) = statDynamics_noUE.totalSettlement(end);
    scenarios.numberFormal_noUE(s) = statDynamics_noUE.totalPrivate(end);
    scenarios.numberRDP_noUE(s) = statDynamics_noUE.totalRDP(end);
    scenarios.priceCBD_noUE(s) = statDynamics_noUE.avgRentFormalCBD(end);            
    scenarios.urbanFootprint_noUE(s) = statDynamics_noUE.urbanFootprint(end);
end


%% Saving

save 'simulations scenarios - 201908'

