function output = InterpolateAgriculturalRentEvolution(option, param, macro, t)
% Rent on the city limit

% There are different ways to make it vary with time

% 1. Agricultural rent proportional to average income (= inflation given our
% scenarios)
% output = ppval(macro.splineIncome, t)./macro.incomeYearReference .* param.agriculturalRent2011;

% 2. Direct spline for agricultural rent
output = ppval(macro.splineAgriculturalRent, t);
option.constructionFunction = 'C-D';
coeffKappaT = InterpolateCoefficientConstruction(option, param, macro, ppval(macro.splineIncome, t));
output = output.^(param.coeff_a) .* (param.depreciationRate + InterpolateInterestRateEvolution(macro, t)) ./ (coeffKappaT .* param.coeff_b .^ param.coeff_b);


end