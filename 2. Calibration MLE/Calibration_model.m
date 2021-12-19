function [score_amenities, score_hous, indexes, parameters_amenities, residual_amenities_mat] = Calibration_model(...
    revenu1, data_rent, data_density,data_dwelling_size, data_income_class, ...
    x_data, y_data, which_general, table_amenities, which_regression)

% Automated calibration of the parameters of NEDUM
%   Using MLE for different values of the parameters


%% Neighbor matrix

[ weighting ] = create_weighting(x_data, y_data, 4);


%% Data as matrices, where should we regress (remove where we have no data)

% Where is which class
which_sample_mat = revenu1 > 0;
for i = 1:4
    which_sample_mat(i,data_income_class ~= i) = false;
end
which_sample_does_exist = sum(which_sample_mat) == 1;


which_rent = ~isnan(data_rent) & which_sample_does_exist & which_general;
which_rent_mat = which_sample_mat(:, which_rent);
revenu_rent = revenu1(:,which_rent);
which_hous = ~isnan(data_dwelling_size) & ~isnan(data_rent) & which_sample_does_exist & which_general;
which_hous_mat = which_sample_mat(:, which_hous);
revenu_hous = revenu1(:,which_hous);

% For the regression
% model_spec = 'residu ~ distance_distr_parks + distance_ocean + distance_urban_herit + airport_cone2 + slope + distance_protected_envir + RDP_proximity'; %+ distance_power_station';
table_reg_global = table_amenities(which_rent,:);
matrix_predict = table2array(table_reg_global(:,which_regression));
matrix_predict = [ones(size(matrix_predict,1),1), matrix_predict];

data_dwelling_size_which = data_dwelling_size(which_hous);
data_rent_which = data_rent(which_rent);
income_class_which = data_income_class(which_rent);

%% Definition of indices

% List of parameters to identify
which_indexes = {'rho', 'beta', 'basic_q', 'Uo1', 'Uo2', 'Uo3', 'Uo4'};
% Uo is a vector of indexes
% The order for indexes is the same as which_indexes

% Rho is the coefficient for spatial autocorrelation
list_rho = 0; 

% Coefficients of the model
list_beta = 0.2:0.02:0.5; %0.1:0.025:0.45; %necessarily between 0 and 1
list_basic_q = 10:1:18;

% Utilities
decomposition_U = 0.7:0.1:1.3;
Uo_cible = [300; 1000; 2800; 12000];
list_Uo1 = Uo_cible(1); % Useless to have variations of this one
list_Uo2 = decomposition_U .* Uo_cible(2); 
list_Uo3 = decomposition_U .* Uo_cible(3);
list_Uo4 = decomposition_U .* Uo_cible(4);

%% Functions of the model 

%%%% Dwelling sizes %%%%

% % Decomposition for the interpolation (the more points, the slower the code)
decomposition_income = [10.^([-8,-4:0.5:-2]),0.03,0.06:0.02:1];

% % Min and Max values for the decomposition
choice_income = max(max(revenu1)) .* decomposition_income ;
income = repmat(choice_income,length(choice_income),1)';
  
% Decomposition for the interpolation
% decomposition_q = [10.^([-9,-2,-1]), 0.2:0.1:1.1, 1.2:0.2:3, 3.5:0.5:13, 15, 20:5:60, 100, 200000]; 
% choice_q = min(list_basic_q) + decomposition_q*10;
% q = repmat(choice_q, length(choice_income),1);

% The function that will be used afterwards
% calcule_hous = @(Uo, beta, basic_q, revenu_temp) griddata(income, utilite_q(q, income, basic_q, beta), q, revenu_temp, Uo);
%function hous = calcule_hous(Uo, beta, basic_q, revenu_temp)
%   hous = griddata(income, utilite_q(q, income, basic_q, beta), q, revenu_temp, Uo);
%   hous(Uo./revenu_temp < (1-beta).^(1-beta).*(min(min(q)) - basic_q) ./ (min(min(q)) - (1-beta).*basic_q)) = basic_q;
%   hous(Uo'*ones(1,size(revenu_temp,2)) > max(max(utilite_q(q, income, basic_q, beta)))) = max(max(q));
%end

% If we estimate calcule_hous directly from data from rents (no extrapolation)
calcule_hous = @(beta, basic_q, revenu_temp, rent_model) beta .* revenu_temp ./ rent_model + (1-beta) .* basic_q;

%%%% Rents %%%%

% Decomposition for the interpolation (the more points, the slower the code)
decomposition_rent = [10.^([-5,-4,-3,-2]),0.02:0.015:0.08,0.1:0.02:1]; 
choice_rent = choice_income./min(list_basic_q); % the maximum rent is the rent for which u = 0
rent = choice_rent' * decomposition_rent;
% 
calcule_rent = @(Uo, beta, basic_q, revenu_temp) griddata(income, utilite_R(rent, income, basic_q, beta), rent, revenu_temp, Uo);

% With indexes
calcule_hous2 = @(indices, revenu_temp, rent_model) calcule_hous(indices(2), indices(3), revenu_temp, rent_model);
calcule_rent2 = @(indices, Uo, revenu_temp) calcule_rent(Uo, indices(2), indices(3), revenu_temp);


% Global score = score of all variables
% calc_score_vrai_fonction is just for rents! With spatial autocorrelation:
compute_score_amenities = @(sigma_amenities, error_amenities)...
         nansum(-log(2*pi)/2 - log(sigma_amenities) - 1/2*(error_amenities).^2 ./ sigma_amenities.^2);
compute_score_hous = @(sigma_hous, error_hous)...
         nansum(-log(2*pi)/2 - log(sigma_hous) - 1/2*(error_hous).^2 ./ sigma_hous.^2);
    

%% Iteration to determine likelihood for all combinations of parameters


f = waitbar(0,'Please wait...');
% Combvec gives all the combination of vectors
indexes = combvec(list_rho, list_beta, list_basic_q, list_Uo1, list_Uo2, list_Uo3, list_Uo4);
length_iterations = size(indexes, 2);
score_amenities = zeros(1, length_iterations);
score_hous = zeros(1, length_iterations);
parameters_amenities = zeros(length(which_regression) + 1, length_iterations);
residual_amenities_mat = zeros(sum(which_rent), length_iterations);
% error_income = zeros(1,length_iterations);

% Start iterations
tic
warning('off')
parfor index_iter = 1:length_iterations
        
        % Error on rents and amenities
        Uo = indexes(4:7, index_iter);
        beta = indexes(2, index_iter);
        residual_amenities = log(Uo) - log((1-beta).^(1-beta).*beta.^beta.*(revenu_rent - indexes(3,index_iter).*repmat(data_rent_which, size(which_rent_mat,1),1)) ./ repmat(data_rent_which, size(which_rent_mat,1),1).^beta);
        residual_amenities = residual_amenities(which_rent_mat);
        residual_amenities(abs(imag(residual_amenities)) > 0) = NaN;
                
        % Regression as matrix
        parameters_amenities(:,index_iter) = matrix_predict(~isnan(residual_amenities),:) \ real(residual_amenities(~isnan(residual_amenities)));
        error_amenities = real(residual_amenities(~isnan(residual_amenities))) - (matrix_predict(~isnan(residual_amenities),:)*parameters_amenities(:,index_iter));
        residual_amenities_mat(:,index_iter) = residual_amenities; 
        
        % fitglm is slow...
        % table_reg = table_reg_global;
        % table_reg.residu = real(residual_amenities);
        % model_amenity = fitglm(table_reg, model_spec);
        % parameters_amenities(:,index_iter) = model_amenity.Coefficients.Estimate;
        % parameters_amenities_SE(:,index_iter) = model_amenity.Coefficients.SE;
        % error_amenities = model_amenity.Residuals.Raw;

        
        % Estimated simulated rents
        % rent_model = calcule_rent2(indexes(index_iter), Uo, revenu_rent);
        % error_rent = log(rent_model(which_rent_mat)) - log(data_rent_which)';

        % Errors on allocation of income groups (to be deleted)
        % [~, income_model] = max(rent_model);
        % error_income(index_iter) = sum((income_model + 1 - income_class_which).^2);
        
        
        % Error on dwelling sizes
        hous = calcule_hous2(indexes(:,index_iter), revenu_hous, data_rent_which);
        hous = hous(which_hous_mat);
        error_hous = log(hous') - log(data_dwelling_size_which);
        
        % Scores
        % score_rent(index_iter) = compute_score_rent(sqrt(nansum(error_rent).^2./sum(~isnan(error_rent))), error_rent);
        score_amenities(index_iter) = compute_score_amenities(sqrt(nansum(error_amenities.^2)./sum(~isnan(error_amenities))), error_amenities);
        score_hous(index_iter) = compute_score_hous(sqrt(nansum(error_hous.^2)./sum(~isnan(error_hous))), error_hous);
        
        % Iterations and waitbar
        if floor(index_iter/50)*50 == index_iter
          time_left = toc/index_iter*(length_iterations - index_iter); %in seconds
          waitbar(index_iter/length_iterations, f, sprintf('%0.1f%% reste %0.0f m %0.0f s total %0.0f m', index_iter/length_iterations*100, floor(time_left/60), time_left - floor(time_left/60)*60, round(toc/60)));          
        end
                                        
end


close(f)
% [score_max,IX] = max(real(score_amenities + score_hous)); % Parameters with max likelihood
% indexes_output = indexes(:,IX)';
% indexes_output_names = which_indexes;

% Amenities
% parameters_amenities = parameters_amenities; %(:,IX);
% SE_amenities = parameters_amenities_SE; %(:,IX);


% Amenities 2
% Uo = indexes(4:7, IX);
% beta = indexes(2, IX);
% residual_amenities = log(Uo) - log((1-beta).^(1-beta).*beta.^beta.*(revenu_rent - indexes(3,IX).*repmat(data_rent_which, size(which_rent_mat,1),1)) ./ repmat(data_rent_which, size(which_rent_mat,1),1).^beta);
% residual_amenities = residual_amenities(which_rent_mat);
% residual_amenities(abs(imag(residual_amenities)) > 0) = NaN;
% table_reg_global.residu = real(residual_amenities);


% table_reg = table_reg_global;
% amenities = exp(parameters_amenities(2:8)' * [table_reg.distance_distr_parks, table_reg.distance_ocean, table_reg.distance_urban_herit, table_reg.airport_cone2, table_reg.slope, table_reg.distance_protected_envir, table_reg.RDP_proximity]') ;
% Uo = indexes(4:6, IX) .* exp(parameters_amenities(1));
% rent = calcule_rent2(indexes(:,IX), Uo./amenities, revenu_rent);
% rent_model = rent(which_rent_mat);
% 
% hous = calcule_hous2(indexes(:,IX), revenu_hous, rent);
% hous_model = hous(which_hous_mat);


% beta = indexes_output(2);
% basic_q = indexes_output(3);
% residual_amenities = log((1-beta).^(1-beta).*beta.^beta .* (revenu_rent - basic_q .* repmat(data_rent(which_rent), size(which_rent_mat,1), 1))./repmat(data_rent(which_rent), size(which_rent_mat,1), 1).^beta);
% residual_amenities = residual_amenities(which_rent_mat);
% residual_amenities(abs(imag(residual_amenities)) > 0) = NaN;
% table_reg = table_reg_global;
% table_reg.residu = real(residual_amenities);
% model_amenity = fitglm(table_reg, model_spec);
% parameters_amenities = model_amenity.Coefficients.Estimate;
% SE_amenities = model_amenity.Coefficients.SE;



disp('*** Estimation of beta and q0 done ***')
% interval = confidence_interval(indexes(:,IX), which_indexes, compute_score_global);

end


function [interval] = confidence_interval(indices_max,quoi_indices,compute_score)
%% Confidence interval

d_beta=1.05;

fprintf('\n');
beta_interval = zeros(size(indices_max));
for index = 1:size(indices_max,1),
    indices_ici = indices_max;
    score_tmp = compute_score(indices_ici);
    
    indices_ici = indices_max;
    score_tmp2 = compute_score(indices_ici);
    
    indices_ici = indices_max;
    indices_ici(index) = indices_ici(index) - indices_ici(index)*(d_beta-1);
    score_tmp3 = compute_score(indices_ici);
    
    indices_ici = indices_max;
    dd_l_beta = -(score_tmp2 + score_tmp3 - 2*score_tmp) / (indices_ici(index)*(d_beta-1))^2;
    beta_interval(index) = 1.96 / (sqrt( abs(dd_l_beta)));
    fprintf('%s\t\t%g (%g ; %g)\n',quoi_indices{index},indices_ici(index),indices_ici(index)-beta_interval(index),indices_ici(index)+beta_interval(index))
end
interval = [indices_ici-beta_interval,indices_ici+beta_interval];

end

function [ score ] = calc_score_vrai_fonction(model_rent, data_rent, sigma1, mu2, sigma2, rho, weighting)
% Log-likelihood for rents


deti = log(abs(det(eye(size(weighting)) - rho*weighting)));

score = log_likelihood_boxcox_SER_reduit_sigma(data_rent, model_rent, sigma1, mu2, sigma2, rho, weighting, deti);


end

function [log_likelihood] = log_likelihood_boxcox_SER_reduit_sigma(data, model, sigma1, mu2, sigma2, rho, weighting, deti)
% Estimates probability with an error term following a normal law (mu,
% sigma) summed with a box cox of parameter lambda

if rho~=0
    data = log(data) - rho*log(data)*weighting;
    model = (model) - rho*(model)*weighting;
else
    data = log(data);
    model = (model);
end

normal_law_log = @(x,mu,sigma_2) - 1/2*log(2*pi*sqrt(sigma_2)) - 1/2*(x-mu).^2/(sigma_2);

test2 = normal_law_log(data,mu2,sigma1^2 + sigma2^2)+...
    log((1 - normcdf(model,(sigma1^2*mu2 + sigma2^2*data) / (sigma1^2 + sigma2^2), sqrt((sigma1^2 + sigma2^2) / sigma1*sigma2))));
test1 = log(normcdf(model,mu2,sigma2)) + normal_law_log(model,data,sigma1.^2);

log_densite_proba = test1;
log_densite_proba(test2>(test1+3)) = test2(test2>(test1+3));
log_densite_proba((test2<=(test1+3))&(test2>=(test1-3))) = log(exp(test1((test2<=(test1+3))&(test2>=(test1-3))))+exp(test2((test2<=(test1+3))&(test2>=(test1-3)))));


log_jacobian2 = deti;

log_likelihood = sum(log_densite_proba) + log_jacobian2;

end

function [log_likelihood] = log_likelihood_just_normal(data,model,sigma,rho,weighting,deti)

%calcvule la proba avec un terme d'erreur suivant une loi normale (mu,
%sigma) additif avec une box cox de parametre lambda

if nargin<6,
    deti=log(abs(det(eye(size(weighting))-rho*weighting)));
end

% data(data==0)=min(data~=0)/100;
% model(model==0)=min(model~=0)/100;
% data=data+0.1;
% model=model+0.1;
if rho~=0
    data=log(data)-rho*log(data)*weighting;
    model=(model)-rho*(model)*weighting;
else
    data=log(data);
    model=(model);
end

% partie2=normpdf(data,mu2,sqrt(sigma1^2+sigma2^2)).*...
%     (1-normcdf(model,(sigma1^2*mu2+sigma2^2*data)/(sigma1^2+sigma2^2),sqrt((sigma1^2+sigma2^2)/sigma1*sigma2)));
% partie1=normcdf(model,mu2,sigma2).*normpdf(model,data,sigma1);
% densite_proba=partie1+partie2;
% %log_densite_proba=(-log(2*pi)/2-log(sigma)-1/2*(error-mu).^2./sigma.^2);
% log_densite_proba=log(densite_proba);

normal_law_log=@(x,mu,sigma_2) -1/2*log(2*pi*sqrt(sigma_2))-1/2*(x-mu).^2/(sigma_2);

log_densite_proba=normal_law_log(data,model,sigma^2);

log_jacobian2=deti;

log_likelihood=sum(log_densite_proba)+log_jacobian2;
end

function [ index_model,index_data ] = selectionne_points( X_data,Y_data,X_model,Y_model,pas )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

if size(X_data)~=size(Y_data),
    return;
end

index_data=true(size(X_data));

%put data as a column vector
if size(X_data,1)<size(X_data,2),
    X_data=X_data';
    Y_data=Y_data';
end

X_data2=repmat(X_data,1,size(X_model,2));
Y_data2=repmat(Y_data,1,size(Y_model,2));

X_model2=repmat(X_model,size(X_data,1),1);
Y_model2=repmat(Y_model,size(Y_data,1),1);

dist=(X_data2-X_model2).^2+(Y_data2-Y_model2).^2;
[rien,index_model]=min(dist,[],2);

%option: we remove points if they are too far from datapoints
if nargin>4,
    choix=(rien<2*pas);
    index_data(~choix)=false;
    index_model=index_model(choix);
end
end

function [ weighting ] = create_weighting( X_data,Y_data,pas )
%create_weighting creates weights matrix for spatial autocorrelation
%analysis
%   Detailed explanation goes here

if size(X_data)~=size(Y_data),
    return;
end

index_data = true(size(X_data));

%put data as a row vector
if size(X_data,1)>size(X_data,2),
    X_data=X_data';
    Y_data=Y_data';
end

X_data2=repmat(X_data,size(X_data,2),1);
Y_data2=repmat(Y_data,size(Y_data,2),1);

X_data3=repmat(X_data',1,size(X_data,2));
Y_data3=repmat(Y_data',1,size(X_data,2));

dist=(X_data2-X_data3).^2+(Y_data2-Y_data3).^2;
weighting=(dist<pas^2);
weighting=logical(weighting-eye(size(weighting)));
end

%% Utilities

function [ utili ] = utilite_q(q, revenu, basic_q, beta)
% Precalculation of the utility
%   [ utili ] = utilite( Ro,revenu,basic_q,param ) estimates
%   utility for a dwelling size "q" and an income (net of transport costs)
%   "revenu". 
% Without amenities 

utili = ((1-beta) .* revenu).^(1-beta)...
         .*(q - basic_q)...
         ./ (q - (1-beta).*basic_q).^(1-beta);
utili(q < basic_q) = 0;

% If the income is 0 then utility is 0 (we need to specify it to avoid
% NaNs later)
utili(revenu == 0) = 0;


end

function [ utili ] = utilite_R(Ro, revenu, basic_q, beta)
% Precalculation of the utility
%   [ utili ] = utilite( Ro,revenu,basic_q,param ) estimates
%   utility for a rent "Ro" and an income (net of transport costs)
%   "revenu". 
% Without amenities 

utili = (1-beta)^(1-beta)*beta^beta...
          .*sign(revenu - basic_q.*Ro)...
          .*abs(revenu - basic_q.*Ro)...
          ./(Ro.^beta);
utili((revenu - basic_q.*Ro) < 0) = 0;

% If the income is 0 then utility is 0 (we need to specify it to avoid
% NaNs later)
utili(revenu==0) = 0;


end
