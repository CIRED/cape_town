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




% Choice parameters for Cape Town

function [param] = ChoiceParameters
% Parameters for the simulation
% For Cape Town

param.yearEquilibrium = 2011;
param.yearReference = 2011;

% Useless
param.yearTax = 2018;
param.Tax = 0;

%%  Utility function

% We use calibrated values below
param.beta = 0.25; 
param.alpha = 1 - param.beta;
param.basic_q = 15; % only for formal dwellings
param.miniLotSize = 30;  


% Heterogeneity of households
param.numberIncomeGroup = 4; 
% Income distribution
param.incomeDistribution = [0 1 1 1 1 2 3 3 4 4 4 4]; % allocation of the 12 groups from the Census to income groups in the simulation
% param.incomeDistribution = [0 1 1 1 1 2 3 3 4 4 4 4]; % allocation of the 12 groups from the Census to income groups in the simulation

% Household size for transport costs
% param.householdSizeTransport = [0.6 1.2 2 2]; 
% Data from Claus: 
param.householdSizeTransport = [1.14 1.94 1.94 1.94]; 



%%  Parameters for transport costs

param.logitFactorMin = 6;
param.logitFactor = 3;

param.timeLimit = 20.2000;
param.meanSlopeSmooth = 3;

param.timeCost = 1;
param.timeCost2 = 1; %param.prix_temps*0.6;


%%  Supply of housing

% Construction function 
% Coeff b is the elasticity of housing supply to capital

% Old value of parameters
% param.coeff_grandA = 0.69; %1.23; %0.4; %0.03;%0.002;0.005;
% param.coeff_b = 0.55; %0.6; %0.78;%0.9;%0.88;

% Coeff estimated from R script
param.coeff_b = 0.35927; %0.32958; %0.37429;
param.coeff_a = 1 - param.coeff_b;
coeffRegression = 10.96326;
param.coeff_grandA = (1./param.coeff_b.^param.coeff_b) .* exp(coeffRegression * param.coeff_b);  % 1.3 .* (1./param.coeff_b.^param.coeff_b) .* exp(10.39611 * param.coeff_b) ; %10.84211

% Discount rate
param.depreciationRate = 0.025;  %0.05; %0.0600;
param.interestRate = 0.0250;

% Agricultural rent
param.agriculturalRent2011 = 807.2; %471.8284; % 533.7; % This is the price per unit of land
param.agriculturalRent2001 = 70.7; % From data

% Limit of the city (density)
param.minDensityUrban = 30000;
param.limitPeopleCityEdge = 40;

% Radius of the city "center" for housing maximum height
param.historicRadius = 100; 
% Maximum housing density in the center
param.limitHeightCenter = 10; % very high => as if there were no limit
% Maximum housing density outside the center
param.limitHeightOut = 10;

% Land occupation max
param.maxCoeffLand = 0.7;

% Parameters for informality
param.sizeShack = 20; %in m2, the size of a backyard shack
param.maxCoeffLandBackyard = 0.45; %0.5408;
param.maxCoeffLandSettlement = 0.4; %0.2;
param.amenityBackyard = 0.38; %0.25;
param.amenitySettlement = 0.37;% 0.25; %0.25;

% Public housing 
param.sizeRDP = 40; %in m2 ; the land area occupied by a RDP house
param.sizeBackyard = 70; %in m2 ; size of the backyard of a RDP house

% Future scenario for public housing
param.sizeBackyardFuture = param.sizeBackyard;
param.futureRatePublicHousing = 5000;


%%  Parameters for dynamic evolution

% Inertia on vertical supply of housing
param.timeInvestHousing = 3;
param.timeDepreciationBuildings = 100;

% If we want to add inertia on horizontal expansion of the city
%param.timeInfraKm = 1;


%% Parameter for the computation

% Maximum of iterations
param.maxIteration = 2500; 

% Precision of the resolution (0.02 = max 2% error)
param.precision = 0.01; %0.02;

param.yearBegin = 2011;

% Threshold (n of jobs) above which we keep the employment center
param.thresholdJobs = 20000;
param.step = 2;



end


