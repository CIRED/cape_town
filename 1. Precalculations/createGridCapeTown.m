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




function output = createGridCapeTown(option)
% Create grid for Cape Town, with useful variables associated

global path_nedum


% Import Grid
importfile([path_nedum,'grid_NEDUM_Cape_Town_500.csv'])

% ID
grid.ID = ID;

% Coordinates of grid centroid (crs = CAPE_NO_19, unit = km).
grid.xCoord = X'/1000;
grid.yCoord = Y'/1000;

grid.sizeSquare = 0.5; % in km

% Coordinates of the CBD
grid.xCenter = -53267.944572790904/1000;
grid.yCenter = -3754855.1309322729/1000; 

% Distance from the CBD
grid.distanceCBD = ((grid.xCoord - grid.xCenter).^2 + (grid.yCoord - grid.yCenter).^2).^0.5;

output = grid;

end