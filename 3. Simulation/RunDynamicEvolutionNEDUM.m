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




function simulation = RunDynamicEvolutionNEDUM(t,initialState,trans,grid,land,poly,param,macro,option)
% Dynamic evolution of NEDUM
% Computes equilibrium for each year, then adds inertia on the building
% stock

iterCalcLite = 1; % One iteration every X 
fprintf(['*** Nedum Cape Town lite: one iteration every %g year ***\n'],(t(2)-t(1))/iterCalcLite);


option.adjustHousingInit = option.adjustHousingSupply;
option.adjustHousingSupply = 0;

% New time step
yearsSimulations = t(1):(t(2)-t(1))/iterCalcLite:t(length(t));

option.ownInitializationSolver = 1;

%matrice de la solution
simulation.dwellingSize = zeros([length(t), size(initialState.dwellingSize)]);
simulation.rent = zeros([length(t), size(initialState.dwellingSize)]);
simulation.households = zeros([length(t), size(initialState.households)]);
simulation.housingSupply = zeros([length(t), size(initialState.dwellingSize)]);

for indexIter = [1:(length(yearsSimulations))]
    
    yearTemp = yearsSimulations(indexIter);
    statTemp = initialState;
    
    if indexIter > 1     
                
        if indexIter == length(t)
            disp('stop')
        end
        
        % Simulation with equilibrium housing stock
        disp('Simulation without constraint');
        option.adjustHousingSupply = 1;
        incomeTemp = InterpolateIncomeEvolution(macro,param,option,grid,poly,yearTemp);
        incomeTemp = incomeTemp(:,1)';
        Uo_unconstrained = (statTemp.utility) ./ statTemp.incomeMatrix(:,1)'.*incomeTemp;
        tmpi = RunEquilibriumSolverNEDUM(yearTemp,trans,option,land,grid,macro,param,poly,Uo_unconstrained);
       
        % Estimation of the derivation of housing supply between t and t+1
        derivHousingTemp = EvolutionHousingSupply(land,param,option, yearsSimulations(indexIter), yearsSimulations(indexIter-1), tmpi.housingSupply(1,:), statTemp.housingSupply(1,:));
        param.housing_in = statTemp.housingSupply(1,:) + derivHousingTemp;

        % Run a new simulation with fixed housing
        disp('Simulation with constraint');
        option.adjustHousingSupply = 0;   
        Uo_simulation = (tmpi.utility + Uo_unconstrained)./2;
        initialState = RunEquilibriumSolverNEDUM(yearTemp,trans,option,land,grid,macro,param,poly,Uo_simulation,param.housing_in);
       
        %Ro de la simulation libre
        statTemp.utility = tmpi.utility;
        statTemp.derivHousing = derivHousingTemp;

    else
        
        statTemp.derivHousing = zeros(size(statTemp.rent(1,:)));
        
    end
    
    
    if (indexIter-1)/iterCalcLite - floor((indexIter-1)/iterCalcLite) == 0 

        simulation.householdsCenter((indexIter-1)/iterCalcLite+1,:,:) = initialState.householdsCenter;
        simulation.householdsHousingType((indexIter-1)/iterCalcLite+1,:,:) = initialState.householdsHousingType;
        simulation.dwellingSize((indexIter-1)/iterCalcLite+1,:,:) = initialState.dwellingSize;
        simulation.rent((indexIter-1)/iterCalcLite+1,:,:) = initialState.rent;
        simulation.households((indexIter-1)/iterCalcLite+1,:,:,:) = initialState.households;
        simulation.error((indexIter-1)/iterCalcLite+1,:,:) = initialState.error;
        simulation.housingSupply((indexIter-1)/iterCalcLite+1,:,:) = initialState.housingSupply;
        simulation.utility((indexIter-1)/iterCalcLite+1,:,:) = initialState.utility;
        simulation.derivHousing((indexIter-1)/iterCalcLite+1,:,:) = statTemp.derivHousing;
        
    end
end

% End of the code

if length(t) < length(yearsSimulations)
    T = t';
else
    T = yearsSimulations';
end
simulation.T = T;


% I set adjustHousing to its initial value
option.adjustHousingSupply = option.adjustHousingInit;


toc()
disp('*** End of computation ***');


end

function output = EvolutionHousingSupply(land,param,option,t1,t0,housingSupply1, housingSupply0)

T = t1';

% Interpolate for the simulation year
housingLimitSimulation = InterpolateHousingLimitEvolution(land,option,param,T);

% New housing supply (accounting for inertia and depreciation w/ time)
if t1 - t0 > 0
    diffHousing = (housingSupply1 - housingSupply0) .* (housingSupply1 > housingSupply0) .* (t1 - t0) ./ param.timeInvestHousing - housingSupply0 .* (t1 - t0)  ./ param.timeDepreciationBuildings;
else
    diffHousing = (housingSupply1 - housingSupply0) .* (housingSupply1 < housingSupply0) .* (t1 - t0) ./ param.timeInvestHousing - housingSupply0 .* (t1 - t0)  ./ param.timeDepreciationBuildings;
end
housingSupplyTarget = housingSupply0 + diffHousing;

% Housing height is limited by potential regulations
housingSupplyTarget = min(housingSupplyTarget, housingLimitSimulation);
minimumHousingSupplyInterp = interp1([2001; 2011; 2100] - param.yearBegin, [zeros(1,length(param.minimumHousingSupply)); param.minimumHousingSupply; param.minimumHousingSupply], t1);
housingSupplyTarget = max(housingSupplyTarget, minimumHousingSupplyInterp);

% Derivation with inertia
outputHousing = housingSupplyTarget - housingSupply0;
output = real(outputHousing);

end

