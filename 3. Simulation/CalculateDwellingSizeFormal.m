function output = CalculateDwellingSizeFormal(Uo,param,trans,income,amenities,solus)
% Determines the numerical relationship between dwelling sizes q and utility
% Only for formal housing

% Stone Geary utility function
income_temp = trans.incomeNetOfCommuting;
income_temp(income_temp < 0) = NaN;
hous = solus(...
   income_temp,...
   Uo'*ones(1,size(income,2))./amenities);

hous(Uo'*ones(1,size(income,2))./amenities > param.max_U) = param.max_q;

output = hous;

end