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




function poly = ImportEmploymentCentersCapeTown(grid, param, option, macro, data, t)
% Using Transport Zones data

global path_nedum


% Import data
[X,Y,TZ2013,Ink1,Ink2,Ink3,Ink4] = ImportTZData(strcat(path_nedum, 'TAZ_amp_2013_proj_centro2.csv'));

TZ2013Centers = importfile_centers(strcat(path_nedum, 'TAZ_final.csv'));
listCenters = TZ2013Centers;
xCenters = zeros(1,length(listCenters));
yCenters = zeros(1,length(listCenters));
corresp = zeros(1,length(listCenters));
increment = 1:length(TZ2013);
for i = 1:length(listCenters)
    xCenters(i) = X(TZ2013 == listCenters(i));
    yCenters(i) = Y(TZ2013 == listCenters(i));
    corresp(i) = increment(TZ2013 == listCenters(i));
end

% Number of employees in each TZ for the 12 income classes
jobsCenters12Class = [zeros(length(Ink1),1) Ink1./3 Ink1./3 Ink1./3 Ink2./2 Ink2./2 Ink3./3 Ink3./3 Ink3./3 Ink4./3 Ink4./3 Ink4./3];

% To use TZ data
% poly.codeCentersInitial = TZ2013;
% poly.CodeCentersInitialSimple = TZ2013;
ID_centre = TZ2013;
xCoord = X;
yCoord = Y;


%% Total number of households per class

yearIncomeDistribution = param.yearBegin + t;
if length(t) == 1
    totalBracket = ppval(macro.splinePopulationIncomeDistribution, t)';
    avgIncomeBracket = ppval(macro.splineIncomeDistribution, t)';
else
    totalBracket = ppval(macro.splinePopulationIncomeDistribution, t);
    avgIncomeBracket = ppval(macro.splineIncomeDistribution, t);
end
avgIncomeGroup = zeros(length(yearIncomeDistribution), param.numberIncomeGroup);
totalGroup = zeros(length(yearIncomeDistribution), param.numberIncomeGroup);
 
%total income distribution in the city
for j=1:param.numberIncomeGroup
     totalGroup(:,j) = sum(totalBracket(:, param.incomeDistribution == j), 2);
     avgIncomeGroup(:,j) = sum(avgIncomeBracket(:,param.incomeDistribution==j).*totalBracket(:, param.incomeDistribution==j), 2)./totalGroup(:,j);
end


%% Duplication of employment centers for the several income groups

poly.codeCentersInitial = zeros(1, param.numberIncomeGroup*length(ID_centre));
xCenters = zeros(1, param.numberIncomeGroup*length(ID_centre));
yCenters = zeros(1,param.numberIncomeGroup*length(ID_centre));
jobsCenters = zeros(length(yearIncomeDistribution),param.numberIncomeGroup*length(ID_centre));
avgIncomeCenters = zeros(length(yearIncomeDistribution), param.numberIncomeGroup*length(ID_centre));
incomeGroupCenters = zeros(length(yearIncomeDistribution), param.numberIncomeGroup*length(ID_centre));

for i = 1:length(ID_centre)
    for j = 1:param.numberIncomeGroup
        poly.codeCentersInitial(1, param.numberIncomeGroup*(i-1)+j) = ID_centre(i);
        xCenters(param.numberIncomeGroup*(i-1)+j) = xCoord(i)/1000;
        yCenters(param.numberIncomeGroup*(i-1)+j) = yCoord(i)/1000;
        jobsCenters(:, param.numberIncomeGroup*(i-1)+j) = repmat(sum(jobsCenters12Class(i,param.incomeDistribution==j)), length(yearIncomeDistribution), 1);
        avgIncomeCenters(:,param.numberIncomeGroup*(i-1)+j) = avgIncomeGroup(:,j);
        incomeGroupCenters(:,param.numberIncomeGroup*(i-1)+j) = j;
    end
end

IDcenter = 1:length(poly.codeCentersInitial);

%% Selection of employment centers to keep

if option.polycentric == 1 
    
    poly.selectedCenters = false(1, length(poly.codeCentersInitial));
    
    % Manual sorting 
    poly.selectedCenters(poly.codeCentersInitial == 5101) = true; % CBD
    poly.selectedCenters(poly.codeCentersInitial == 2002) = true; % Bellville
    poly.selectedCenters(poly.codeCentersInitial == 1201) = true; % Epping
    poly.selectedCenters(poly.codeCentersInitial == 1553) = true; % Claremont
    poly.selectedCenters(poly.codeCentersInitial == 3509) = true; % Sommerset West
    %poly.selectedCenters(poly.codeCentersInitial == 5508) = true; % Century City
    %poly.selectedCenters(poly.codeCentersInitial == 5526) = true; % Table View
    poly.selectedCenters(poly.codeCentersInitial == 5523) = true; % Table View + Century City
    
    if option.sortEmploymentCenters == 1
        poly.selectedCenters(jobsCenters(1,:) <= 0) = false;
    end
    
else
    
    % We only keep the CBD
    poly.selectedCenters = false(1, length(poly.codeCentersInitial));
    poly.selectedCenters(poly.codeCentersInitial == 5101) = true; %CBD

end



%% Rescale of number of jobs after selection

% Rescale to include for each center all the jobs located in a defined buffer zone
distanceBuffer = 4;
jobsCentersTemp = jobsCenters(:,poly.selectedCenters);
sumIncomeGroup = zeros(length(yearIncomeDistribution), param.numberIncomeGroup);
distanceCenters = sqrt((xCenters - xCenters(poly.selectedCenters)').^2 + (yCenters - yCenters(poly.selectedCenters)').^2);

for i = 1:length(yearIncomeDistribution)
    for j = 1:param.numberIncomeGroup
        jobs_i = jobsCenters(i,:);
        group_i = incomeGroupCenters(i,:);
        group_selected_i = incomeGroupCenters(i,poly.selectedCenters);
        jobsCentersTemp(i, group_selected_i == j) = jobs_i(group_i == j) * (distanceCenters(group_selected_i == j,group_i == j)' < distanceBuffer);
        sumIncomeGroup(i,j) = sum(jobsCentersTemp(i, group_selected_i == j));
    end
end

% Remove the employment centers that are not significant enough
IDCentersSelected = IDcenter(poly.selectedCenters);
whichTemp = ones(1,sum(poly.selectedCenters));
if option.sortEmploymentCenters == 1
    for j = 1:param.numberIncomeGroup
        whichTemp(jobsCentersTemp(1, group_selected_i == j) ./ sumIncomeGroup(1,j) < 0.1) = 0;
    end 
end 

whichTemp = logical(whichTemp);
jobsCentersTemp = jobsCentersTemp(:,whichTemp);
IDCenterRemoved = IDCentersSelected(whichTemp == 0);
poly.selectedCenters(ismember(IDcenter, IDCenterRemoved)) = 0;
poly.xCenter = xCenters(poly.selectedCenters);
poly.yCenter = yCenters(poly.selectedCenters);
poly.incomeGroup = incomeGroupCenters(:,poly.selectedCenters);
poly.averageIncome = avgIncomeCenters(:,poly.selectedCenters);

% Rescale to keep the correct global income distribution

sumHouseholdsSelectedCenters = zeros(length(yearIncomeDistribution), param.numberIncomeGroup);
jobsCentersSelected = zeros(size(jobsCentersTemp));

for j = 1:param.numberIncomeGroup
    
    sumHouseholdsSelectedCenters(:,j) = sum(jobsCentersTemp(:, poly.incomeGroup(1,:) == j), 2);
    numberCentersSelected = sum(poly.incomeGroup(1,:) == j);
    
    jobsCentersSelected(:,poly.incomeGroup(1,:) == j)...
        = jobsCentersTemp(:,poly.incomeGroup(1,:) == j) .*...
            repmat(totalGroup(:, j)./sumHouseholdsSelectedCenters(:, j), 1, numberCentersSelected);
        
end

% Export 

yearCenters = yearIncomeDistribution; 
poly.totalHouseholdsGroup = totalGroup;

poly.year = yearCenters;

increment = 1:length(poly.averageIncome(1,:));
poly.increment = increment;

poly.codeCentersPolycentric = poly.codeCentersInitial(poly.selectedCenters);
poly.averageIncome = avgIncomeCenters(:,poly.selectedCenters);
poly.incomeGroup = incomeGroupCenters(:,poly.selectedCenters);

increment = 1:length(TZ2013);
corresp = zeros(1,length(poly.xCenter));
for i = 1:length(poly.xCenter)
    corresp(i) = increment(TZ2013==poly.codeCentersPolycentric(i));
end
poly.corresp = corresp';

poly.jobsCentersMemory = jobsCenters;
poly.jobsCenters = jobsCentersSelected;

%% Area of influence of each center (for calibration of amenities)

% List of employment centers and their 'catchment' area 
poly.listCenter = unique(poly.codeCentersPolycentric,'stable');
poly.listCatchment = sqrt((xCoord./1000 - unique(poly.xCenter, 'stable')).^2 + (yCoord./1000 - unique(poly.yCenter, 'stable')).^2) < distanceBuffer;

% Loading OD_matrix
load([path_nedum,'Emme_OD_matrices.mat'])
destinationWeight = zeros(length(poly.listCenter), length(poly.listCatchment(:,1)));
for i = 1:length(poly.listCenter)
     destinationWeight(i,:) = nansum(OD_matrix(:,poly.listCatchment(:,i)'),2)';
end
% poly.ODMatrixAll = OD_matrix; 
poly.destinationWeight = destinationWeight ./ nansum(destinationWeight,1);

%% Defining who can live in each type of housing

% Who can live in informal / formal settlements
poly.formal = [1 1 1 1];
poly.backyard = [1 1 0 0];
poly.settlement = [1 1 0 0];




end