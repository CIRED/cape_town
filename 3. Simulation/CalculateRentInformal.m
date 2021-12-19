function R_mat = CalculateRentInformal(Uo,param,transTemp,income,amenity)
% Determines the numerical relationship between bid-rents and utility
% Only for informal housing, for which the size of the dwelling is fixed

R_mat = 1./param.sizeShack.*(transTemp.incomeNetOfCommuting - (repmat(Uo',1,size(income,2))./(amenity.*(param.sizeShack - param.basic_q).^param.beta)).^(1./param.alpha));

end