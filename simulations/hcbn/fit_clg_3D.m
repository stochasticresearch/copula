function [ x1_domain, x1_multinomial_est, mle_params ] = fit_clg_3D( x )
%FIT_CONDITIONAL_LINEAR_GAUSSIAN_DIST Fits a mixed distribution to a
%conditional linear gaussian distribution.
% Inputs:
%  x - 3-D data, with the first dimension being the discrete distribution
%      and the second and third being continuous
% 
% Outputs:
%  x1_domain  - the domain of the discrete random variable
%  mle_params - a vector of the MLE parameters for each Gaussian
%               corresponding to a value in the domain of the discrete RV.
%               x(1,1) = u_1, x(2,1) = sigma_1 ... 

x1 = x(:,1);
x2 = x(:,2);
x3 = x(:,3);

% estimate mle parameters for relationship between X_1 and X_2
[x1_domain, x1_multinomial_est, x1x2_mle_params] = fit_clg_2D([x1 x2]);

% estimate mle parameters for relationship between X_1 and X_3
[~, ~, x1x3_mle_params] = fit_clg_2D([x1 x3]);

% estimate mle_parameters for relationship between X_2 and X_3
[mu_x2, sigma_x2] = normfit(x2);
[mu_x3, sigma_x3] = normfit(x3);
rho = corrcoef([x2 x3]);

mle_params = {};
mle_params{1} = x1x2_mle_params;
mle_params{2} = x1x3_mle_params;
mle_params{3} = {};
mle_params{3}{1} = [mu_x2; sigma_x2];
mle_params{3}{2} = [mu_x3; sigma_x3];
mle_params{3}{3} = rho(1,2);

end

