function [distance_metro, duration_metro] = ImportSimulatedTransportTimeByTrain(poly, grid, param)
% import and estimate transport time by the metro
global path_nedum
importfile([path_nedum,'metro_station_poly.csv']);

station_line_time = [Bellvill1_B Bellvill2_M Bellvill3_S Bonteheuwel1_C Bonteheuwel2_B Bonteheuwel3_K Capeflats Malmesbury Simonstown Worcester];

duration = zeros(length(ID_station));
for i = 1:length(ID_station) %matrice des temps O-D entre les stations de m?tro
   for j = 1:i
       if i==j
           duration(i,j)=0;
       elseif sum(station_line_time(i,:).*station_line_time(j,:)) > 0 %pas besoin de faire un changement
           temps = abs(station_line_time(j,station_line_time(i,:).*station_line_time(j,:)>0) - station_line_time(i,station_line_time(i,:).*station_line_time(j,:)>0));
           duration(i,j) = min(temps) + param.trainWaitingTime;
           duration(j,i) = duration(i,j);
       else %il faut faire un changement
           line_i = station_line_time(i,:)>0;
           line_j = station_line_time(j,:)>0;
           noeud = false(length(ID_station),1);
           for k = 1:length(ID_station)
             if sum(station_line_time(k,:).*station_line_time(i,:))>0 && sum(station_line_time(k,:).*station_line_time(j,:))>0
                noeud(k) = true;
             end    
           end
           temps1 = (abs(repmat(station_line_time(j,line_j>0),sum(noeud),1) - station_line_time(noeud,line_j>0)));
           temps2 = (abs(repmat(station_line_time(i,line_i>0),sum(noeud),1) - station_line_time(noeud,line_i>0)));
           duration(i,j) = min(min(temps1,[],2)+min(temps2,[],2));
           duration(i,j) = duration(i,j) + 2.*param.trainWaitingTime;
           duration(j,i) = duration(i,j);
       end
   end
end

% For each grid point we take the closest station, compute distance
ID_station_grille = griddata(X_cape./1000,Y_cape./1000,ID_station,grid.xCoord, grid.yCoord,'nearest');
distance_grille = zeros(length(grid.distanceCBD),1);

% For each employment center we take the closest station, compute distance
ID_station_center = griddata(X_cape./1000, Y_cape./1000, ID_station, poly.xCenter, poly.yCenter,'nearest');
distance_center = zeros(length(poly.xCenter),1);
for i=1:length(poly.xCenter)
    distance_center(i) = sqrt((poly.xCenter(i) - X_cape(ID_station_center(i))./1000).^2+(poly.yCenter(i) - Y_cape(ID_station_center(i))./1000).^2);
end

% Matrix of distances and durations
duration_metro = zeros(length(grid.distanceCBD),length(poly.xCenter));
distance_metro = zeros(length(grid.distanceCBD),length(poly.xCenter));
for i=1:length(grid.distanceCBD)
   distance_grille(i) = sqrt((grid.xCoord(i)-X_cape(ID_station_grille(i))./1000).^2+(grid.yCoord(i)-Y_cape(ID_station_grille(i))./1000).^2); 
   for j=1:length(poly.xCenter)
       duration_metro(i,j) = (distance_grille(i)+distance_center(j)).*1.2./(param.walkingSpeed./60) + duration(ID_station_grille(i),ID_station_center(j));
       distance_metro(i,j) = max(sqrt((grid.xCenter-X_cape(ID_station_grille(i))./1000).^2+(grid.yCenter-Y_cape(ID_station_grille(i))./1000).^2), sqrt((grid.xCenter-X_cape(ID_station_center(j))./1000).^2+(grid.yCenter-Y_cape(ID_station_center(j))./1000).^2));
   end
end
duration_metro = duration_metro';
distance_metro = distance_metro';

end 