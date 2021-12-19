%% Bid rents
% We estimate bid rents using precalculated matrices
switch typeHousing
    
case 'formal'
    [R_mat] = CalculateBidRentFormalBasicNeed(Uo, param, transTemp, grid, income, amenities, solus_R);
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

% Rents
[R,whichMax] = max(R_mat(:,:),[],1);

%quel is the type of households with the highest bid-rent in each location

%% Dwelling sizes

temp = [0:size(transTemp.incomeNetOfCommuting,2)-1]*size(transTemp.incomeNetOfCommuting,1);
whichMaxTemp = whichMax+temp; 

switch typeHousing
    
case 'formal'
    dwellingSize = param.beta.*(transTemp.incomeNetOfCommuting(whichMaxTemp))./R + param.alpha.*param.basic_q;

    % If we want to add the cost associated with time (as it is not paid in
    % practice):
    % hous = hous+param.coeff_beta*trans_tmp.prix_temps(quel_mat)./R;
case 'backyard'
    dwellingSize = param.sizeShack.*ones(size(whichMaxTemp));
case 'informal'
    dwellingSize = param.sizeShack.*ones(size(whichMaxTemp));

end

