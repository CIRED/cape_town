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




function [z2,a2,weight_sum,L,U] = GenerateSmoothMean(X,Y,weight,precision)

% Removing NaNs
whichRemove = (isnan(X))|(isnan(Y))|(isnan(weight));
whichRemove = ~whichRemove;
X = X(whichRemove);
Y = Y(whichRemove);

if length(weight)>1
    weight=weight(whichRemove);
end
if length(precision)>1
    precision=precision(whichRemove);
end

% Weighting or not?
if length(weight)==1
    if weight==1
        weight = ones(size(X));
    end
end

% Approximation (precision is the scale
A = round(X/precision)*(precision);
% Sorting A
[a_sort,IX] = sort(A);
% Sorting Y
z_sort = Y(IX);
% We regroupe unique values for A
[a2, m, n] = unique(a_sort, 'first');

leng = size(m,1);
m(leng+1) = size(a_sort,1)+1;

% z2 is the output (moving average)
z2 = zeros(size(a2));
% weight_sum sum of weigths
weight_sum = zeros(size(a2));
L = zeros(size(a2));
U = zeros(size(a2));
for index = 1:leng
    pop_tot = sum(weight(m(index):m(index+1)-1));
    z2(index) = sum(z_sort(m(index):m(index+1)-1) .* weight(m(index):m(index+1)-1)) / pop_tot;
    L(index) = z2(index) - quantile(z_sort(m(index):m(index+1)-1),0.25);
    U(index) = -z2(index) + quantile(z_sort(m(index):m(index+1)-1),0.75);
    weight_sum(index) = pop_tot;
end

end