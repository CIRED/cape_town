%% Dwelling sizes

switch typeHousing
    
case 'formal'
        
    dwellingSize = CalculateDwellingSizeFormal(Uo, param, transTemp, income, amenities, solus_Q);
    
    % If we want to add the cost associated with time (as it is not paid in
    % practice):
    % dwellingSize = hous + param.beta*transTemp.timeCost./R;
    
    % Here we introduce the minimum lot-size 
    dwellingSize = max(dwellingSize, param.miniLotSize, 'includenan');
    dwellingSize(poly.formal == 0,:) = NaN;
    
case 'backyard'
    dwellingSize = param.sizeShack.*ones(length(poly.incomeMult), length(grid.distanceCBD));
    dwellingSize(poly.backyard == 0,:) = NaN;
    
case 'informal'
    dwellingSize = param.sizeShack.*ones(length(poly.incomeMult), length(grid.distanceCBD));
    dwellingSize(poly.settlement == 0,:) = NaN;
    
end

%% Bid rents
% We estimate bid rents from dwelling sizes

switch typeHousing
    
case 'formal'
    R_mat = param.beta.*(transTemp.incomeNetOfCommuting)./(dwellingSize - param.alpha.*param.basic_q); 
    R_mat(transTemp.incomeNetOfCommuting < 0) = 0;
    R_mat(poly.formal == 0,:) = 0;
    
case 'backyard'
    amenities = amenities.*param.amenityBackyard;
    [R_mat] = CalculateRentInformal(Uo,param,transTemp,income,amenities);
    R_mat(poly.backyard == 0,:) = 0;
    
case 'informal'
    amenities = amenities.*param.amenitySettlement;
    [R_mat] = CalculateRentInformal(Uo,param,transTemp,income,amenities);
    R_mat(poly.settlement == 0,:) = 0;

end

R_mat = single(R_mat);
R_mat(R_mat<0) = 0;

