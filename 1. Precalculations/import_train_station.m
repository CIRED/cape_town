function [X_train_station,Y_train_station,ID_train_station]=import_train_station(param)
%import the location of train stations

global path_nedum
importfile([path_nedum,'metro_station_poly.csv']);
ID_train_station=ID_station;
X_train_station=X_cape./1000;
Y_train_station=Y_cape./1000;

end