function housingSupply = CalculateHousingSupplyFormal(R, option, limitHousing, constructionParameter, param, agriculturalRent, referenceRent, interestRate, limitProfil)
% Calculates the housing construction as a function of rents

if option.adjustHousingSupply == 1
    
    if option.constructionFunction == 'C-D'
        housingSupply = 1000000 .* constructionParameter.^(1/param.coeff_a)...
            .*(param.coeff_b/interestRate).^(param.coeff_b/param.coeff_a)...
            .*(R).^(param.coeff_b/param.coeff_a);
    elseif option.constructionFunction == 'CES'
        housingSupply = 1000000 .* constructionParameter .* param.coeff_a_CES^(-param.coeff_sigma_CES/(1-param.coeff_sigma_CES)) .*...
            (1 - param.coeff_b_CES^param.coeff_sigma_CES .* (constructionParameter .* R/interestRate).^...
            (param.coeff_sigma_CES - 1)).^(param.coeff_sigma_CES/(1 - param.coeff_sigma_CES));
    end
    
    % Outside the agricultural rent, no housing (accounting for a tax)
    housingSupply(R < agriculturalRent + param.taxUrbanEdgeMat) = 0;
    
    housingSupply(isnan(housingSupply)) = 0;
    housingSupply(imag(housingSupply) ~= 0) = 0;
    housingSupply(housingSupply < 0) = 0;
    housingSupply = min(housingSupply, (ones(size(housingSupply,1),1)*limitHousing));
    
    % To add the construction on Mitchells_Plain
    housingSupply = max(housingSupply, param.minimumHousingSupply .* 1000000);
    
else    
    
    housingSupply = param.housing_in;
    
end

end