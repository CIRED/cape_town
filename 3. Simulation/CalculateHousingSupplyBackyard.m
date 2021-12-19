function housing = CalculateHousingSupplyBackyard(R, grid, param, basic_q_formal, income, transportCostRDP)
% Calculates the backyard available for construction as a function of rents

housing = param.alpha.*(param.sizeRDP + param.sizeBackyard - basic_q_formal)./(param.sizeBackyard) - param.beta.*(income(1,:) - transportCostRDP)./((param.sizeBackyard).*R);
housing(income(1,:) < transportCostRDP) = param.alpha.*(param.sizeRDP + param.sizeBackyard - basic_q_formal)./(param.sizeBackyard) - param.beta.*(income(1, income(1,:) < transportCostRDP))./((param.sizeBackyard).*R(income(1,:) < transportCostRDP));
housing(R == 0) = 0;
housing = min(housing, 1);
housing = max(housing, 0);
housing = 1000000 .* housing;

end