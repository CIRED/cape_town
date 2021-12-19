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






function [option] = ChoiceOptions
% Ste the options for the simulations

%% Assumptions of the model

% Polycentric or monocentric
option.polycentric = 1;
option.sortEmploymentCenters = 1;

% Whether there is construction of RDP houses in the future
option.futureConstructionRDP = 1;

% Logit for the modal allocation and ``cross-commuting''
option.logit = 1;

% The function for housing supply (Cobb-Douglas or CES)
option.constructionFunction = 'C-D'; % 'C-D' or 'CES'

% Initialization of the solver
option.ownInitializationSolver = 1;
option.adjustHousingSupply = 1;

% For speed, need to be set to 0 if there is a change in the employment
% centers or the grid
option.loadTransportTime = 1;

%% Scenarios for prospective simulation

% Whether we keep the Urban Edge in the future simulations
option.urbanEdge = 0; %1 means we keep the urban edge

% Scenario numbers
option.scenarioPop = '2'; %(1, 2, 3)
option.scenarioIncomeDistribution = '2'; %(1, 2, 3)
option.scenarioInflation = '1'; %(1)
option.scenarioInterestRate = '1'; %(1)
option.scenarioPriceFuel = '1'; %(1)

% For now, cannot be used
option.taxOutUrbanEdge = 0;
option.doubleStoreyShacks = 0;



end




