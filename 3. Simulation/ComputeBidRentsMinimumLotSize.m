%% Dwelling sizes

switch typeHousing
    
case 'formal'
        
    dwellingSize = (repmat(Uo', 1, length(grid.dist))./amenity).^(1/param.beta).*(param.alpha.*abs(transTemp.incomeNetOfCommuting)).^(-param.alpha/param.beta);
    dwellingSize(transTemp.incomeNetOfCommuting < 0) = 0;
    
    % If we want to add the cost associated with time (as it is not paid in
    % practice):
    % hous = hous+param.beta*trans_tmp.prix_temps./R;
    
    % Here we introduce the minimum lot-size 
    dwellingSize = max(dwellingSize, param.miniLotSize, 'includenan');
    dwellingSize(poly.formal == 0,:) = NaN;
    
case 'backyard'
    dwellingSize = param.sizeShack.*ones(length(poly.income_mult), length(grid.dist));
    dwellingSize(poly.backyard == 0,:) = NaN;
    
case 'informal'
    dwellingSize = param.size_shack.*ones(length(poly.income_mult), length(grid.dist));
    dwellingSize(poly.settlement == 0,:) = NaN;
    
end

%% Bid rents
% We estimate bid rents using precalculated matrices
switch typeHousing
    
case 'formal'
    R_mat = param.beta.*(transTemp.incomeNetOfCommuting)./dwellingSize; 
    R_mat(income < transTemp.cout_generalise) = 0;
    R_mat(poly.formal == 0,:) = 0;
    
case 'backyard'
    amenity = amenity.*param.amenityBackyard;
    [R_mat] = CalculateRentInformal(Uo,param,transTemp,income,amenity);
    R_mat(poly.backyard == 0,:) = 0;
    
case 'informal'
    amenity = amenity.*param.amenitySettlement;
    [R_mat] = CalculateRentInformal(Uo,param,transTemp,income,amenity);
    R_mat(poly.settlement == 0,:) = 0;

end

R_mat = single(R_mat);
R_mat(R_mat<0) = 0;

% Rents
[R,whichMax] = max(R_mat(:,:),[],1);

% whichMax is the type of households with the highest bid-rent in each location
temp = [0:size(transTemp.incomeNetOfCommuting,2)-1]*size(transTemp.incomeNetOfCommuting,1);
whichMaxTemp = whichMax+temp; 

% Dwelling size of the highest bidder
dwellingSize = dwellingSize(whichMaxTemp);

