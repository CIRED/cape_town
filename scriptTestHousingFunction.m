%% Script to test CES housing function and compare with Cobb-Douglas
% For now, simple simulations


% Assume a negative exponential function for rents
x = 0:1:60;
R = 100.*exp(-0.04*x);

% Parameters
b = 0.25;
a = 0.75;
K = 10; % scale, does not matter a lot here
rho_delta = 0.1;

% Cobb-Douglas
H_cd = K^(1/a) .* (b.*R./(rho_delta)).^(b/a);

% CES
b = 1.90;
a = 0.90;
sigma = 0.75; % sigma = 0.75; % elasticity of substitution
K_CES = K*10;
H_ces = K_CES.*a^(-sigma/(1-sigma)) .* (1 - b^sigma.*(rho_delta./(K_CES.*R)).^(1-sigma)).^(sigma/(1-sigma));

% CES 2
b = 0.90;
a = 0.1;
sigma = 0.75; % sigma = 0.75; % elasticity of substitution
K_CES = K/36;
H_ces_2 = K_CES.*a^(-sigma/(1-sigma)) .* (1 - b^sigma.*(rho_delta./(K_CES.*R)).^(1-sigma)).^(sigma/(1-sigma));


% Compare two functions
close all
plot(R,H_cd)
hold on 
plot(R,H_ces)
plot(R,H_ces_2)


% Compare two functions
% close all
% plot(x,H_cd)
% hold on 
% plot(x,H_ces)

