
% Import grid  
grid = createGridCapeTown(option);

% Macro trajectories
macro = LoadModelInputsCapeTown(param, option);

% Import of the main datasets
data = LoadDataCapeTown(grid, param);

% Income groups and employment centers
poly = ImportEmploymentCentersCapeTown_LOGIT(grid, param, option, macro, data, t);

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

% In Mitchells Plain, housing supply is given exogenously (planning), and
% household of group 2 live there (Coloured neighborhood). 
param.minimumHousingSupply = zeros(1,length(grid.distanceCBD));
param.minimumHousingSupply(data.MitchellsPlain) = data.gridFormalDensityHFA(data.MitchellsPlain)./land.coeffLand(1,data.MitchellsPlain);
param.minimumHousingSupply(land.coeffLand(1,:) < 0.1 | isnan(param.minimumHousingSupply)) = 0;
param.multiProbaGroup = NaN .* ones(param.numberIncomeGroup, length(grid.distanceCBD));
% param.multiProbaGroup(:, data.MitchellsPlain) = [0;1;0;0] * ones(1, sum(data.MitchellsPlain));

%%%  Import Calibrated parameters %%% 

% Load amenities
load('./0. Precalculated inputs/calibratedAmenities.mat')
land.amenities = amenities ./ mean(amenities);

% Load utility function coefficients
load('./0. Precalculated inputs/calibratedUtility_beta.mat')
load('./0. Precalculated inputs/calibratedUtility_q0.mat')
param.beta = calibratedUtility_beta;
param.alpha = 1 - param.beta; 
param.basic_q = calibratedUtility_q0;

% Load coefficients for construction function
load('./0. Precalculated inputs/calibratedHousing_b.mat')
load('./0. Precalculated inputs/calibratedHousing_kappa.mat')
param.coeff_b = coeff_b;
param.coeff_a = 1 - param.coeff_b;
param.coeff_grandA = coeffKappa;
load('./0. Precalculated inputs/calibratedHousing_b_CES')
load('./0. Precalculated inputs/calibratedHousing_kappa_CES')
load('./0. Precalculated inputs/calibratedHousing_sigma_CES')
param.coeff_b_CES = coeff_b_CES;
param.coeff_a_CES = 1 - param.coeff_b_CES;
param.coeff_grandA_CES = coeffKappa_CES;
param.coeff_sigma_CES = coeffSigma_CES;

% Load amenities for informal housing
load('./0. Precalculated inputs/calibratedParamAmenities.mat')
param.amenityBackyard = calibratedParamAmenities(1);
param.amenitySettlement = calibratedParamAmenities(2);

% Define minimum lot-size 
param.miniLotSize = min(data.spDwellingSize((data.spInformalSettlement + data.spInformalBackyard)./data.spTotalDwellings < 0.1));
% We transform the rent per unit of land to a rent per unit of floor area:
% param.agriculturalRent2011 = param.agriculturalRent2011.^(param.coeff_a) .* (param.depreciationRate + InterpolateInterestRateEvolution(macro, 0)) ./ (param.coeff_grandA .* param.coeff_b .^ param.coeff_b);

% Import local incomes and lambda
load('./0. Precalculated inputs/lambda.mat')
param.lambda = lambdaKeep;
load('./0. Precalculated inputs/incomeCentersKeep.mat')
poly.incomeCentersInit = incomeCentersKeep;

carto2 = @(x) MakeMapCape(grid, poly, x);
