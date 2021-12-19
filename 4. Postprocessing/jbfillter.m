% Copyright or © or Copr. Ecole des Ponts ParisTech 2012
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




%sert pour tracer les courbes (faire des surfaces transparentes)

if nargin<7;transparency=.5;end %default is to have a transparency of .5
if nargin<6;add=1;end     %default is to add to current plot
if nargin<5;edge='k';end  %dfault edge color is black
if nargin<4;color='b';end %default color is blue

if length(upper)==length(lower) && length(lower)==length(xpoints)
    %msg='';
    filled=[upper,fliplr(lower)];
    xpoints=[xpoints,fliplr(xpoints)];
%     if add
%         hold on
%     end
    fillhandle=fill3(xpoints,filled,-1*ones(size(xpoints)),color*(transparency)+[1 1 1]*(1-transparency));%plot the data
    transparency_bak=transparency;
    transparency=1;
    if nargin==8,
    set(fillhandle,'EdgeColor','None','FaceAlpha',transparency_bak,'EdgeAlpha',transparency_bak,'DisplayName',name);%set edge color
    else
        set(fillhandle,'EdgeColor','None','FaceAlpha',transparency_bak,'EdgeAlpha',transparency_bak);%set edge color
    end
%     if add
%         hold off
%     end

% fillhandle=fill3(xpoints,filled,1*ones(size(xpoints)),color);%plot the data
%     transparency=transparency_bak;
%     if nargin==8,
%     set(fillhandle,'EdgeColor','None','FaceAlpha',transparency,'EdgeAlpha',transparency,'DisplayName',name);%set edge color
%     else
%         set(fillhandle,'EdgeColor','None','FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
%     end
else
    msg='Error: Must use the same number of points in each vector';
end

