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
disp ('**************** Scenario BAU with and without Urban Edge *****************');


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
option.polycentric = 1;
option.futureConstructionRDP = 1;
option.logit = 1;
option.taxOutUrbanEdge = 0;
option.doubleStoreyShacks = 0;
option.urbanEdge = 0; %1 means we keep the urban edge
option.ownInitializationSolver = 1;
option.adjustHousingSupply = 1;
option.loadTransportTime = 1;
option.sortEmploymentCenters = 1;
option.constructionFunction = 'C-D'; % 'C-D' or 'CES'

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
yearTrafic = [t(1):sign(t(length(t)))*2:(t(length(t)))];
trans = LoadTransportCostCapeTown_LOGIT(option, grid, macro, param, poly, data, yearTrafic, 'GR', 1);
disp ('*** Data and parameters imported succesfully ***')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation with urban edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option.urbanEdge = 1; %1 means we keep the urban edge
param.yearUrbanEdge = 2015; %in case option.urban_edge = 0, the year the constraint is removed

land_UE = ImportLandUseInputsCapeTown(grid, option, param, macro, data);
land_UE.amenities = land.amenities;

option.ownInitializationSolver = 0;
Uo_init = 0;
option.ajust_bati = 1;

initialState = RunEquilibriumSolverNEDUM_LOGIT(t(1),trans,option,land_UE,grid,macro,param,poly,Uo_init,param.housing_in);
statInitialState = ComputeStatInitialState_LOGIT(trans, land_UE, grid, macro, param, option, initialState, poly, t(1));

disp('*** beginning of evolution - urban edge ***');
simulation_UE = RunDynamicEvolutionNEDUM_LOGIT(t,initialState,trans,grid,land_UE,poly,param,macro,option);
statDynamics_UE = ComputeFinalStatisticsCapeTown(macro,option,trans,land_UE,grid,param,poly,simulation_UE);
disp('*** end of evolution ***');


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation with NO urban edge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option.urbanEdge = 0; %1 means we keep the urban edge
param.yearUrbanEdge = 2015; %in case option.urban_edge = 0, the year the constraint is removed

land_noUE = ImportLandUseInputsCapeTown(grid, option, param, macro, data);
land_noUE.amenities = land.amenities;

option.ownInitializationSolver = 0;
Uo_init = 0;
option.ajust_bati = 1;

disp('*** beginning of evolution - no urban edge ***');
simulation_noUE = RunDynamicEvolutionNEDUM_LOGIT(t,initialState,trans,grid,land_noUE,poly,param,macro,option);
statDynamics_noUE = ComputeFinalStatisticsCapeTown(macro,option,trans,land_noUE,grid,param,poly,simulation_noUE);
disp('*** end of evolution ***');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation with Low public housing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option.urbanEdge = 0; %1 means we keep the urban edge
param.yearUrbanEdge = 2015; %in case option.urban_edge = 0, the year the constraint is removed

param.futureRatePublicHousing = 2500; % Number per year
macro_lowPH = LoadModelInputsCapeTown(param, option);
land_lowPH = ImportLandUseInputsCapeTown(grid, option, param, macro_lowPH, data);
land_lowPH.amenities = land.amenities;

disp('*** beginning of evolution - low public housing ***');
simulation_lowPH = RunDynamicEvolutionNEDUM_LOGIT(t,initialState,trans,grid,land_lowPH,poly,param,macro_lowPH,option);
statDynamics_lowPH = ComputeFinalStatisticsCapeTown(macro_lowPH,option,trans,land_lowPH,grid,param,poly,simulation_lowPH);
disp('*** end of evolution ***');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation with High public housing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option.urbanEdge = 0; %1 means we keep the urban edge
param.yearUrbanEdge = 2015; %in case option.urban_edge = 0, the year the constraint is removed

param.futureRatePublicHousing = 10000; % Number per year
macro_highPH = LoadModelInputsCapeTown(param, option);
land_highPH = ImportLandUseInputsCapeTown(grid, option, param, macro_highPH, data);
land_highPH.amenities = land.amenities;

disp('*** beginning of evolution - high public housing ***');
simulation_highPH = RunDynamicEvolutionNEDUM_LOGIT(t,initialState, trans, grid, land_highPH, poly, param, macro_highPH, option);
statDynamics_highPH = ComputeFinalStatisticsCapeTown(macro_highPH, option, trans, land_highPH, grid, param, poly, simulation_highPH);
disp('*** end of evolution ***');

%% Saving

save 'simulations - 201908'

