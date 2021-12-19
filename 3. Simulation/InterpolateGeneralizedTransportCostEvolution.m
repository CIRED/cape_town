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




function output = InterpolateGeneralizedTransportCostEvolution(trans,param,t)
% computes transport generalized cost for a given year, by interpolation, using variable trans as input

for index=1:length(t)
    
    
    [index1,index2,ponder1,ponder2] = CreatePonderation(t(index)+param.yearBegin, trans.yearTransport);
    
    output(:,:,index) = ponder1*trans.generalizedCost(:,:,index1)...
        +ponder2*trans.generalizedCost(:,:,index2);
    
end
end



function [index1,index2,ponder1,ponder2] = CreatePonderation(value,vector)
vectorCenter = vector - value;

[valueMin,index] = min(abs(vectorCenter));

if valueMin == 0
    index1 = index;
    index2 = index;
    ponder1 = 1;
    ponder2 = 0;
else
    vecteurNeg = vectorCenter;
    vecteurNeg(vecteurNeg>0) = NaN;
    [close1,index1] = max(vecteurNeg);
    
    vecteurPos = vectorCenter;
    vecteurPos(vecteurPos<0) = NaN;
    [close2,index2] = min(vecteurPos);
    
    ponder1 = abs(close1)/(close2 - close1);
    ponder2 = 1-ponder1;
end
end