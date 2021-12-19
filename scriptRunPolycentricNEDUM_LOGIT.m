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
disp ('**************** NEDUM-Cape-Town - Polycentric Version - Formal and Informal housing ****************');

global slash
slash = '/'; % To be changed when on Mac or Windows
% Load function addresses
pathProgram;
% Load data address
pathData;

%% Parameters and options

%%% Import Parameters %%%
param = ChoiceParameters;
option = ChoiceOptions;
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
option.constructionFunction = 'CES'; % 'C-D' or 'CES'

%%% Parameters for the scenario %%%
param.yearUrbanEdge = 2016; % in case option.urban_edge = 0, the year the constraint is removed
param.taxOutUrbanEdge = 10000;
param.coeffDoubleStorey = 0.02; % Random for now

%%% Years for the simulation %%%
param.yearBegin = 2011;
t = 0:1:2; 

%%% Import inputs for the simulation %%%
ImportInputsData_LOGIT;

% Transport data
yearTrafic = [t(1):sign(t(length(t)))*2:t(length(t))];
%%
trans = LoadTransportCostCapeTown_LOGIT(option, grid, macro, param, poly, data, yearTrafic, 'GR', 1);
disp ('*** Data and parameters imported succesfully ***')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initial state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if option.ownInitializationSolver
    load('./0. Precalculated inputs/calibratedUtilities.mat')
    Uo_init  = [1000;10000;utilitiesCorrected] .* 4;
else
    Uo_init = 0;
end

%%% Solver %%%
initialState = RunEquilibriumSolverNEDUM_LOGIT(t(1),trans,option,land,grid,macro,param,poly,Uo_init,param.housing_in);

disp('*** End of static resolution ***');

%%% Compute initial statistics %%%
statInitialState = ComputeStatInitialState_LOGIT(trans, land, grid, macro, param, option, initialState, poly, t(1));


%% Plot Graphs for the Initial State and Comparison with the Data

nameFolder = strcat('..', slash, 'results_Cape_Town', slash, 'calibration_072019');

%%% Graphs for main outputs of the model %%%
PlotResultsInitialState(grid, initialState, statInitialState, macro, data, land, param, macro, option, t(1), nameFolder)
PlotResultsHousingInitial(grid, initialState, statInitialState, data, land, param, macro, poly, option, t(1), nameFolder)

%%% Export data for maps in R, save parameters %%%
ExportDataForMapsInitialState(grid, poly, land, initialState, statInitialState, data, nameFolder)
save(strcat(nameFolder, slash, 'parameters.mat'), 'param')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Test for the uniqueness of the equilibrium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% /!\ long to run /!\

% New sets of precision
param.maxIteration = 5000; 
param.precision = 0.01; %0.02;
option.ownInitializationSolver = 1;

Uo_1 = 10.^(2:3:8);
Uo_init_vect = combvec(Uo_1, Uo_1, Uo_1, Uo_1);

util_equi = zeros(size(Uo_init_vect));
impossible_equi = zeros(size(Uo_init_vect));
error_equi = zeros(size(Uo_init_vect));

for i = 1:size(Uo_init_vect, 2)
    initialState_temp = RunEquilibriumSolverNEDUM_LOGIT(t(1),trans,option,land,grid,macro,param,poly,Uo_init_vect(:,i),param.housing_in);
    util_equi(:,i) = initialState_temp.utility';
    impossible_equi(:,i) = initialState_temp.impossiblePopulation';
    error_equi(:,i) = initialState_temp.error';
end

% maximum error in percentage
for i = 1:param.numberIncomeGroup
    max_error(i) = max(abs((util_equi(i,:) ./ mean(util_equi(i,:)) - 1) .* 100));
end
