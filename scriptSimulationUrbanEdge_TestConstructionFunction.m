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

%%% Parameters for the scenario %%%
param.yearUrbanEdge = 2015; % in case option.urban_edge = 0, the year the constraint is removed
param.taxOutUrbanEdge = 10000;
param.coeffDoubleStorey = 0.02; % Random for now

%%% Years for the simulation %%%
param.yearBegin = 2011;
t = 0:1:28; % Up to 2039

%%% Import inputs for the simulation %%%
ImportInputsData_LOGIT;

% Transport data
yearTrafic = [t(1):sign(t(length(t)))*2:(t(length(t)))];
trans = LoadTransportCostCapeTown_LOGIT(option, grid, macro, param, poly, data, yearTrafic, 'GR', 1);
disp ('*** Data and parameters imported succesfully ***')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation with NO urban edge and Cobb-Douglas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option.urbanEdge = 0; %1 means we keep the urban edge
param.yearUrbanEdge = 2015; %in case option.urban_edge = 0, the year the constraint is removed

land_noUE = ImportLandUseInputsCapeTown(grid, option, param, macro, data);
land_noUE.amenities = land.amenities;

option.ownInitializationSolver = 0;
Uo_init = 0;
option.ajust_bati = 1;

option.constructionFunction = 'C-D'; % 'C-D' or 'CES'

initialState = RunEquilibriumSolverNEDUM_LOGIT(t(1),trans,option,land_noUE,grid,macro,param,poly,Uo_init,param.housing_in);
statInitialState = ComputeStatInitialState_LOGIT(trans, land_noUE, grid, macro, param, option, initialState, poly, t(1));


disp('*** beginning of evolution - no urban edge ***');
simulation_CD = RunDynamicEvolutionNEDUM_LOGIT(t,initialState,trans,grid,land_noUE,poly,param,macro,option);
statDynamics_CD = ComputeFinalStatisticsCapeTown(macro,option,trans,land_noUE,grid,param,poly,simulation_CD);
disp('*** end of evolution ***');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation with NO urban edge and CES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

option.constructionFunction = 'CES'; % 'C-D' or 'CES'

initialState = RunEquilibriumSolverNEDUM_LOGIT(t(1),trans,option,land_noUE,grid,macro,param,poly,Uo_init,param.housing_in);
statInitialState = ComputeStatInitialState_LOGIT(trans, land_noUE, grid, macro, param, option, initialState, poly, t(1));

disp('*** beginning of evolution - no urban edge ***');
simulation_CES = RunDynamicEvolutionNEDUM_LOGIT(t,initialState,trans,grid,land_noUE,poly,param,macro,option);
statDynamics_CES = ComputeFinalStatisticsCapeTown(macro,option,trans,land_noUE,grid,param,poly,simulation_CES);
disp('*** end of evolution ***');
