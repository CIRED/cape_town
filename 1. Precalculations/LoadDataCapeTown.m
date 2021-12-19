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

function data = LoadDataCapeTown(grid, param)
% Function to import main datasets useful to compare outputs

global path_nedum
global slash

% Importing data takes time: if option.LOAD_DATA == 0, we just load an
% saved version of data_courbe

%%%%%
loadData = 0;
%%%%%

if loadData == 0
    % If the data are already imported
    load(strcat('.', slash, '0. Precalculated inputs', slash, 'data'))
    disp('Data loaded directly');
    
else    

%% donn?es issues du recensement 2011

fprintf('importing Census data...');

% Subplace data from the Census 2011
importfile([path_nedum,'sub-places-dwelling-statistics.csv']);
data.spCode = SP_CODE;
data.spX = CoordX./1000;
data.spY = CoordY./1000;
data.sp2011Distance = sqrt((data.spX - grid.xCenter).^2+(data.spY - grid.yCenter).^2);
data.sp2011Area = ALBERS_ARE; % in km2
data.sp2011CapeTown = (MN_CODE == 199);
data.sp2011MitchellsPlain = (MP_CODE == 199039);

% Income distribution
data.gridAverageIncome = InterpolateDataSPtoGrid(Data_census_dwelling_INC_average, data, grid);
data.sp2011AverageIncome = Data_census_dwelling_INC_average;
data.sp2011IncomeDistribution12Class = [Data_census_dwelling_INC_0 Data_census_dwelling_INC_1_4800	Data_census_dwelling_INC_4801_9600 Data_census_dwelling_INC_9601_19600 Data_census_dwelling_INC_19601_38200 Data_census_dwelling_INC_38201_76400 Data_census_dwelling_INC_76401_153800 Data_census_dwelling_INC_153801_307600 Data_census_dwelling_INC_307601_614400 Data_census_dwelling_INC_614001_1228800 Data_census_dwelling_INC_1228801_2457600 Data_census_dwelling_INC_2457601_more];
data.sp2011IncomeDistributionNClass = zeros(length(SP_CODE), param.numberIncomeGroup);
for i = 1:param.numberIncomeGroup
    data.sp2011IncomeDistributionNClass(:,i) = sum(data.sp2011IncomeDistribution12Class(:, param.incomeDistribution == i), 2);
end
middleClass = floor(param.numberIncomeGroup./2);
poor = param.incomeDistribution <= middleClass;
rich = param.incomeDistribution > middleClass;

% Interpolation to the grid
data.gridNumberPoor = InterpolateDataSPtoGridAlternative(sum(data.sp2011IncomeDistribution12Class(:,poor),2), data, grid);
data.gridNumberRich = InterpolateDataSPtoGridAlternative(sum(data.sp2011IncomeDistribution12Class(:,rich),2), data, grid);
data.MitchellsPlain = InterpolateDataSPtoGridAlternative(data.sp2011MitchellsPlain, data, grid) > 0;

% Income distribution at the City level
importfile([path_nedum,'Income_distribution_2011.csv']);
data.medianIncome = INC_med;
for j = 1:param.numberIncomeGroup
    data.thresholdIncomeDistribution(j) = max(INC_max(param.incomeDistribution==j)); 
end


% Type de logement
data_type = [Data_census_dwelling_House_concrete_block_structure Data_census_dwelling_Traditional_dwelling Data_census_dwelling_Flat_apartment Data_census_dwelling_Cluster_house	Data_census_dwelling_Townhouse Data_census_dwelling_Semi_detached_house Data_census_dwelling_House_flat_room_in_backyard Data_census_dwelling_Informal_dwelling_in_backyard	Data_census_dwelling_Informal_dwelling_settlement Data_census_dwelling_Room_flatlet_on_property_or_larger_dw Data_census_dwelling_Caravan_tent Data_census_dwelling_Other Data_census_dwelling_Unspecified Data_census_dwelling_Not_applicable];
data.spTotalDwellings = sum(data_type,2);
data.spInformalBackyard = data_type(:,8);
data.spInformalSettlement = data_type(:,9);
data.spInformalBackyard(isnan(data.spInformalBackyard)) = 0;
data.spInformalSettlement(isnan(data.spInformalSettlement)) = 0;
data.gridInformalBackyard = InterpolateDataSPtoGridAlternative(data.spInformalBackyard, data, grid);
data.gridInformalSettlement = InterpolateDataSPtoGridAlternative(data.spInformalSettlement, data, grid);
data.gridFormal = InterpolateDataSPtoGridAlternative(data.spTotalDwellings - data.spInformalSettlement - data.spInformalBackyard, data, grid);




%% Donn?es densit? pop 2001

importfile([path_nedum,'Census_2001_income.csv'])
data.sp2001X = X_2001./1000;
data.sp2001Y = Y_2001./1000;
data.sp2001Code = SP_CODE;
data.sp2001DistanceCenter = sqrt((data.sp2001X - grid.xCenter).^2+(data.sp2001Y - grid.yCenter).^2);
data.sp2001Area = Area_sqm./1000000; %in km2
data.sp2001CapeTown = (data.sp2001Code > 17000000) & (data.sp2001Code < 18000000);

% Income Classes
data.sp2001Distribution12Class = [Census_2001_inc_No_income Census_2001_inc_R1_4800 Census_2001_inc_R4801_9600 Census_2001_inc_R9601_19200 Census_2001_inc_R19201_38400 Census_2001_inc_R38401_76800 Census_2001_inc_R76801_153600 Census_2001_inc_R153601_307200 Census_2001_inc_R307201_614400 Census_2001_inc_R614401_1228800 Census_2001_inc_R1228801_2457600 Census_2001_inc_R2457601_more];
for i=1:param.numberIncomeGroup
    data.sp2001DistributionNClass(:,i) = sum(data.sp2001Distribution12Class(:,param.incomeDistribution==i),2);
end

% Density of people
data.sp2001NumberPoor = sum(data.sp2001Distribution12Class(:,poor),2);
data.sp2001NumberRich = sum(data.sp2001Distribution12Class(:,rich),2); %We keep the same categories than for 2011 (is it relevant?)
% data.grid2001Households = InterpolateDataSP2001toGridAlternative(data.sp2001NumberPoor + data.sp2001NumberRich, data, grid);


importfile([path_nedum,'Census_2001_dwelling_type.csv'])
formal = House_brick_structure_separate_stand + Flat_in_block + semi_detached_house + House_flat_in_backyard + Room_flatlet_shared_property + Caravan_tent + Ship_boat;
backyard = 	Informal_dwelling_in_backyard;
informal = Traditional_dwelling_traditional_materials + Informal_dwelling_NOT_backyard;

data.sp2001Formal = zeros(1,length(data.sp2001Code));
data.sp2001Backyard = zeros(1,length(data.sp2001Code));
data.sp2001Settlement = zeros(1,length(data.sp2001Code));
for i = 1:length(data.sp2001Code)
    match = SP_Code == data.sp2001Code(i);
    data.sp2001Formal(i) = formal(match);
    data.sp2001Backyard(i) = backyard(match);
    data.sp2001Settlement(i) = informal(match);
end

data.grid2001Formal = InterpolateDataSP2001toGridAlternative(data.sp2001Formal, data, grid);
data.grid2001Backyard = InterpolateDataSP2001toGridAlternative(data.sp2001Backyard, data, grid);
data.grid2001Settlement = InterpolateDataSP2001toGridAlternative(data.sp2001Settlement, data, grid);


%% Total number of people per income class

% For 2001 and 2011
data.totalNumberPerIncomeGroup = [sum(data.sp2001DistributionNClass(data.sp2001CapeTown,:)); sum(data.sp2011IncomeDistributionNClass(data.sp2011CapeTown, :))];
data.totalNumberPerIncomeBracket = [sum(data.sp2001Distribution12Class(data.sp2001CapeTown,:)); sum(data.sp2011IncomeDistribution12Class(data.sp2011CapeTown, :))];


%% Data on dwelling prices in 2011
fprintf('importing real estate data...');

% Data from the Sales data of the City of Cape Town
% Sales data were previously treaty and aggregated at the SP level on R
% (see other script)
importfile([path_nedum,'SalePriceStat_SP.csv'])
data.spPrice=zeros(3,length(data.spCode));
for i = 1:length(data.spCode)
   if sum(SP_CODE==data.spCode(i)) == 1
     data.spPrice(3,i) = Median_2011(SP_CODE == data.spCode(i));
     data.spPrice(2,i) = Median_2006(SP_CODE == data.spCode(i));
     data.spPrice(1,i) = Median_2001(SP_CODE == data.spCode(i));
   end
end
data.yearPrice = [2001 2006 2011];
data.spPrice(data.spPrice==0) = NaN;
data.xPrice = data.spX';
data.yPrice = data.spY';
data.distancePrice = sqrt((data.xPrice - grid.xCenter).^2+(data.yPrice - grid.yCenter).^2);


%% Data dwelling sizes and available area per Subplace (aggregated with R)

importfile([path_nedum, 'SP_dwelling_size_area_urb_constrained.csv']);
data.spDwellingSize = 10000 .* ones(1,length(data.spCode));
data.spFormalDensityHFA = zeros(1,length(data.spCode));
data.SP_is_CT = zeros(1,length(data.spCode));
data.spUnconstrainedArea = zeros(1,length(data.spCode));
for i=1:length(data.spCode)
   if sum(SP_CODE == data.spCode(i)) == 1
     which = SP_CODE == data.spCode(i);
     data.spDwellingSize(i) = MEAN_RES_D(which);
     data.spFormalDensityHFA(i) = dens_HFA(which);
     data.SP_is_CT(i) = is_Cape_Town(which);
     data.spUnconstrainedArea(i) = unconstrained(which);
   end
end
data.spDwellingSize(data.spDwellingSize == 10000) = NaN;

% data.gridDwellingSize = InterpolateDataSPtoGrid(data.spDwellingSize, data, grid);
data.gridFormalDensityHFA = InterpolateDataSPtoGrid(data.spFormalDensityHFA, data, grid);
data.limitCapeTown = InterpolateDataSPtoGrid(data.SP_is_CT, data, grid) > 0;
                                 
%% RDP houses from GV2012

importfile([path_nedum, 'GV2012_grid_RDP_count2.csv']);
data.gridCountRDPfromGV = count_RDP';
data.gridAreaRDPfromGV = area_RDP';


%% Household size as a function of income groups 

% Data from 2011 Census (Claus' estimation)
data.householdSizeIncomeGroup = [6.556623149, 1.702518978, 0.810146856, 1.932265222];

save(strcat('.', slash, '0. Precalculated inputs', slash, 'data'), 'data')

end

end


function [dataGrid]=...
    InterpolateDataSPtoGrid(data_SP, data, grid)

global path_nedum

importfile([path_nedum,'grid_SP_intersect.csv'])

dataGrid = zeros(1,length(grid.distanceCBD));

for index=1:length(grid.distanceCBD)
    
        listCode = unique(SP_CODE(ID_grille==grid.ID(index)));
        areaExcluded = 0;
        for i=1:length(listCode)
            if isempty(data_SP(data.spCode==listCode(i)))
                areaExcluded = areaExcluded + sum(Area(ID_grille==grid.ID(index) & SP_CODE==listCode(i)));
            else
            dataGrid(index)=dataGrid(index)+sum(Area(ID_grille==grid.ID(index) & SP_CODE==listCode(i)))...
            *data_SP(data.spCode==listCode(i));
            end
        end
        if areaExcluded>0.9*sum(Area(ID_grille==grid.ID(index)))
            dataGrid(index)=NaN;
        else
            dataGrid(index)=dataGrid(index)/(sum(Area(ID_grille==grid.ID(index))) - areaExcluded);
        end 
            
end


end


function [dataGrid]=...
    InterpolateDataSPtoGridAlternative(data_SP, data_courbe, grid)
% For "extensive" variables

global path_nedum

importfile([path_nedum,'grid_SP_intersect.csv'])

dataGrid = zeros(1,length(grid.distanceCBD));

for index=1:length(grid.distanceCBD)
    
        ici = unique(SP_CODE(ID_grille == grid.ID(index)));
       
        for i=1:length(ici)
            if ~isempty(data_SP(data_courbe.spCode == ici(i)))
                dataGrid(index) = dataGrid(index) + sum(Area(ID_grille == grid.ID(index) & SP_CODE == ici(i)))...
                *data_SP(data_courbe.spCode == ici(i)) / sum(Area(SP_CODE == ici(i)));
            end
        end
            
end


end

function [dataGrid]=...
    InterpolateDataSP2001toGridAlternative(data_SP, data_courbe, grid)
% For "extensive" variables

global path_nedum

importfile([path_nedum,'grid_SP2001_intersect.csv'])

dataGrid = zeros(1,length(grid.distanceCBD));

for index=1:length(grid.distanceCBD)
    
        ici = unique(SP_CODE(ID_grille == grid.ID(index)));
       
        for i=1:length(ici)
            if ~isempty(data_SP(data_courbe.sp2001Code == ici(i)))
                dataGrid(index) = dataGrid(index) + sum(area_intersection(ID_grille == grid.ID(index) & SP_CODE==ici(i)))...
                *data_SP(data_courbe.sp2001Code == ici(i)) / sum(area_intersection(SP_CODE==ici(i)));
            end
        end
            
end


end
