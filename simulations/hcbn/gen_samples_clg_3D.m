function [ X1X2, X1X3, X2X3 ] = gen_samples_clg_3D( x1_domain, x1_multinomial_est, mle_params, N )
%GEN_SAMPLES_CONDITIONAL_LINEAR_GAUSSIAN Generates samples of a conditional
%linear gaussian model
% Input:
%  x1_domain - the domain of the discrete random variable
%  x1_multinomial_est - the estimated frequency of each value in the
%                       x1_domain occuring
%  mle_params - a cell array containing the parameters of the distributions
%  N - the number of samples to generate
% 
% Outputs:
%  X - the generated samples

% generate X1-X2
X1X2 = gen_samples_clg_2D( x1_domain, x1_multinomial_est, mle_params{1}, N );

% generate X1-X3
X1X3 = gen_samples_clg_2D( x1_domain, x1_multinomial_est, mle_params{2}, N );

% generate X2-X3
p = mle_params{3};
p_x2 = p{1};
p_x3 = p{2};
rho = p{3};
MU = [p_x2(1) p_x3(1)];
SIGMA = [p_x2(2).^2 rho*p_x2(2)*p_x3(2); rho*p_x2(2)*p_x3(2) p_x3(2).^2];
X2X3 = mvnrnd(MU,SIGMA,N);

end

