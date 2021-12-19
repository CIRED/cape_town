function [ output_args ] = MakeMapCape(grille,poly, input)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

close all
hold on
scatter(grille.xCoord,grille.yCoord,510,input,'.');
scatter(poly.xCenter, poly.yCenter, sum(poly.jobsCenters(:,:,1),2)./400, 'o', 'w');
ylim([-3810 -3695])

end

