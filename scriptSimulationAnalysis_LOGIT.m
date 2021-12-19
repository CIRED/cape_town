%% Script to analyse the results of scenarios

nameFileGeneral = strcat('..', slash, 'results_Cape_Town', slash, '201907 - simulations', slash);

%% Comparison between Urban Edge and no Urban Edge

nameFile = strcat(nameFileGeneral, 'UE_no_UE', slash);


% Longitudinal charts for scenario analysis

statDynamics_UE = ComputeFinalStatisticsCapeTown(macro,option,trans,land_UE,grid,param,poly,simulation_UE);
statDynamics_noUE = ComputeFinalStatisticsCapeTown(macro,option,trans,land_noUE,grid,param,poly,simulation_noUE);

PlotResultsComparativeDynamicScenario(param, macro, simulation_UE, statDynamics_UE, land_UE, simulation_noUE, statDynamics_noUE, land_noUE, 'Urban edge', 'No urban edge', nameFile);

% Comparison of final states

PlotResultsComparativeFinalStates(param, option, grid, macro, poly, trans, simulation_UE, statDynamics_UE, land_UE, simulation_noUE, statDynamics_noUE, land_noUE, 'Urban edge', 'No urban edge', nameFile)
close all

% Export for maps
yearUrbanized = 2011;
ExportDataForMaps(grid, param, land_UE, data, poly, t, simulation_UE, yearUrbanized, nameFile, 'Urban Edge');
ExportDataForMaps(grid, param, land_noUE, data, poly, t, simulation_noUE, yearUrbanized, nameFile, 'No Urban Edge');


%% Comparison between Low public housing and high public housing (no Urban Edge)

nameFile = strcat(nameFileGeneral, 'Public housing', slash);


% Longitudinal charts for scenario analysis

statDynamics_lowPH = ComputeFinalStatisticsCapeTown(macro_lowPH,option,trans,land_lowPH,grid,param,poly,simulation_lowPH);
statDynamics_highPH = ComputeFinalStatisticsCapeTown(macro_highPH,option,trans,land_highPH,grid,param,poly,simulation_highPH);
PlotResultsComparativeDynamicScenario(param, macro_lowPH, simulation_lowPH, statDynamics_lowPH, land_lowPH,...
    simulation_highPH, statDynamics_highPH, land_highPH, 'RDP/BNG - Low scenario', 'RDP/BNG - High scenario', nameFile);

% Comparison of final states

PlotResultsComparativeFinalStates(param, option, grid, macro, poly, trans, simulation_lowPH, statDynamics_lowPH, land_lowPH,...
    simulation_highPH, statDynamics_highPH, land_highPH, 'RDP/BNG - Low scenario', 'RDP/BNG - BAU', nameFile)
close all

% Export for maps
yearUrbanized = 2011;
ExportDataForMaps(grid, param, land_lowPH, data, poly, t, simulation_lowPH, yearUrbanized, nameFile, 'Low public housing');
ExportDataForMaps(grid, param, land_highPH, data, poly, t, simulation_highPH, yearUrbanized, nameFile, 'High public housing');

