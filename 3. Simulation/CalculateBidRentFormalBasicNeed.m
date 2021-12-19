function sortie = CalculateBidRentFormalBasicNeed(Uo,param,transTemp,grid,income,amenities,solus)
% Determines the numerical relationship between bid-rents and utility
% Only for formal housing

% Stone Geary utility function
R_mat = solus(...
   income - double(transTemp.cout_generalise),...
   Uo'*ones(1,size(income,2))./amenities);


sortie = R_mat;

end