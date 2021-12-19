
% Import grid  
grid = createGridCapeTown(option);

% Macro trajectories
macro = LoadModelInputsCapeTown(param, option);

% Import of the main datasets
data = LoadDataCapeTown(grid, param);
macro.splinePopulationIncomeDistribution = interp1([2001; 2011; 2100] - param.yearBegin, [sum(data.sp2001Distribution12Class); sum(data.sp2011IncomeDistribution12Class); sum(data.sp2011IncomeDistribution12Class)], 'linear', 'pp');

% Income groups and employment centers
poly = ImportEmploymentCentersCapeTown(grid, param, option, macro, data, t);

% Housing types: 
% 1 - Formal
% 2 - Informal in Backyard
% 3 - Informal in Settlements
% 4 - RDP/BNG housing (exogenously given)

% Import coefficient for land use
land = ImportLandUseInputsCapeTown(grid, option, param, macro, data);

% Construction 
param.housing_in = data.gridFormalDensityHFA./land.coeffLand(1,:).*1.1;
param.housing_in(~isfinite(param.housing_in)) = 0;
param.housing_in(param.housing_in>2*10^6) = 2*10^6;
param.housing_in(param.housing_in<0) = 0;

% In Mitchells Plain, housing supply is given exogenously (planning). 
param.minimumHousingSupply = zeros(1,length(grid.distanceCBD));
param.minimumHousingSupply(data.MitchellsPlain) = data.gridFormalDensityHFA(data.MitchellsPlain)./land.coeffLand(1,data.MitchellsPlain);
param.minimumHousingSupply(land.coeffLand(1,:) < 0.1 | isnan(param.minimumHousingSupply)) = 0;

% Transport data
yearTrafic = [t(1):sign(t(length(t))):t(length(t))];
trans = LoadTransportCostCapeTown(option, grid, macro, param, poly, data, yearTrafic, 'SP', 1);

carto2 = @(x) MakeMapCape(grid, poly, x);

disp ('*** Data and parameters imported succesfully ***')