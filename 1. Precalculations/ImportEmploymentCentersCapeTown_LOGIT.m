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




function poly = ImportEmploymentCentersCapeTown_LOGIT(grid, param, option, macro, data, t)
% Using Transport Zones data

global path_nedum


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import data

% Import data
[X,Y,TZ2013,Ink1,Ink2,Ink3,Ink4] = ImportTZData(strcat(path_nedum, 'TAZ_amp_2013_proj_centro2.csv'));

% Number of employees in each TZ for the 12 income classes
jobsCenters12Class = [zeros(length(Ink1),1) Ink1./3 Ink1./3 Ink1./3 Ink2./2 Ink2./2 Ink3./3 Ink3./3 Ink3./3 Ink4./3 Ink4./3 Ink4./3];

% To use TZ data
poly.codeCentersInitial = TZ2013;
ID_centre = TZ2013;
xCoord = X ./ 1000;
yCoord = Y ./ 1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Total number of households per class 

yearIncomeDistribution = param.yearBegin + t;
if length(t) == 1
    totalBracket = ppval(macro.splinePopulationIncomeDistribution, t)';
    avgIncomeBracket = ppval(macro.splineIncomeDistribution, t)';
else
    totalBracket = ppval(macro.splinePopulationIncomeDistribution, t);
    avgIncomeBracket = ppval(macro.splineIncomeDistribution, t);
end

% Total income distribution in the city
avgIncomeGroup = zeros(length(yearIncomeDistribution), param.numberIncomeGroup);
totalGroup = zeros(length(yearIncomeDistribution), param.numberIncomeGroup);
for j=1:param.numberIncomeGroup
     totalGroup(:,j) = sum(totalBracket(:, param.incomeDistribution == j), 2);
     avgIncomeGroup(:,j) = sum(avgIncomeBracket(:,param.incomeDistribution==j).*totalBracket(:, param.incomeDistribution==j), 2)./totalGroup(:,j);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Selection of employment centers to keep

poly.selectedCenters = sum(jobsCenters12Class,2) > 2500;

% Where we don't have reliable transport data
poly.selectedCenters(xCoord > -10) = false;
poly.selectedCenters(yCoord > -3719) = false;
poly.selectedCenters(xCoord > -20 & yCoord > -3765) = false;
poly.selectedCenters(poly.codeCentersInitial == 1010) = false;
poly.selectedCenters(poly.codeCentersInitial == 1012) = false;
poly.selectedCenters(poly.codeCentersInitial == 1394) = false;
poly.selectedCenters(poly.codeCentersInitial == 1499) = false;
poly.selectedCenters(poly.codeCentersInitial == 4703) = false;

poly.xCenter = xCoord(poly.selectedCenters);
poly.yCenter = yCoord(poly.selectedCenters);

% Number of workers per group for the selected 
jobsCentersNgroup = zeros(length(xCoord), param.numberIncomeGroup);
for j=1:param.numberIncomeGroup
     jobsCentersNgroup(:,j) = sum(jobsCenters12Class(:, param.incomeDistribution == j), 2);
end
jobsCentersNgroup = jobsCentersNgroup(poly.selectedCenters,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rescale of number of jobs after selection

% Rescale to keep the correct global income distribution
jobsCentersNGroupRescaled = zeros(size(jobsCentersNgroup,1), size(jobsCentersNgroup, 2), length(yearIncomeDistribution));
for i = 1:length(yearIncomeDistribution)
    jobsCentersNGroupRescaled(:,:,i) = jobsCentersNgroup .* totalGroup(i,:) ./ sum(jobsCentersNgroup, 1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Export 

yearCenters = yearIncomeDistribution; 
poly.totalHouseholdsGroup = totalGroup;

poly.year = yearCenters;

poly.codeCentersPolycentric = poly.codeCentersInitial(poly.selectedCenters);
poly.averageIncomeGroup = avgIncomeGroup;

increment = 1:length(TZ2013);
poly.increment = increment(poly.selectedCenters);

corresp = zeros(1,length(poly.xCenter));
for i = 1:length(poly.xCenter)
    corresp(i) = increment(TZ2013==poly.codeCentersPolycentric(i));
end
poly.corresp = corresp';

poly.jobsCentersMemory = jobsCenters12Class;
poly.jobsCenters = jobsCentersNGroupRescaled;

[~,whichYearInit] = min(abs(param.yearBegin - poly.year));
poly.jobsCenterInit(:,:) = poly.jobsCenters(:,:,whichYearInit);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining who can live in each type of housing

% Who can live in informal / formal settlements
poly.formal = [1 1 1 1];
poly.backyard = [1 1 0 0];
poly.settlement = [1 1 0 0];
poly.incomeGroup = repmat([1,2,3,4], length(poly.year), 1);



end