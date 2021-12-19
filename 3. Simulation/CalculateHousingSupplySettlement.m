function housingSupply = CalculateHousingSupplySettlement(R, grid, param, poly, income, transportCost, proba)
% Calculates the informal available for construction as a function of rents

netIncome = sum(proba(poly.class == 1, :) .* (income(poly.class == 1, :) - transportCost(poly.class == 1, :))) ./ sum(proba(poly.class == 1, :));
housingSupply = 1 + param.alpha ./ param.coeff_mu - netIncome ./ R;
housingSupply = max(housingSupply, 1);
housingSupply = min(housingSupply, 2);
housingSupply(R == 0) = 0;
housingSupply = 1000000 .* housingSupply;

end