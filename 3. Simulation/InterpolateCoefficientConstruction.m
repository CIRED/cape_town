function output = InterpolateCoefficientConstruction(option,param,macro,income)

% Set the scale coefficient of the construction function: varies over time as the
% average income, power -b 
if option.constructionFunction == 'C-D'
    coeff_A = param.coeff_grandA;
    coeff_b = param.coeff_b;
elseif option.constructionFunction == 'CES'
    coeff_A = param.coeff_grandA_CES;
    coeff_b = param.coeff_b;
end
output = (income./macro.incomeYearReference).^(-coeff_b).*coeff_A;

end