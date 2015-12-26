function [ x1_domain, x1_multinomial_est, x2_mle_params ] = fit_clg_2D( x )
%FIT_CONDITIONAL_LINEAR_GAUSSIAN_DIST Fits a mixed distribution to a
%conditional linear gaussian distribution.
% Inputs:
%  x - 2-D data, with the first dimension being the discrete distribution
%      and the second being continuous
% 
% Outputs:
%  x1_domain  - the domain of the discrete random variable
%  mle_params - a vector of the MLE parameters for each Gaussian
%               corresponding to a value in the domain of the discrete RV.
%               x(1,1) = u_1, x(2,1) = sigma_1 ... 

x1 = x(:,1);
x2 = x(:,2);

n = length(x1);

x1_domain = unique(x1);
h = histogram(x1,length(x1_domain));
x1_multinomial_est = h.Values/n;

x2_mle_params = zeros(2,length(x1_domain));

for ii=1:length(x1_domain)
    q_i = x1_domain(ii);
    % grab data associated with this q_i
    idx = x1==q_i;
    data = x2(idx);
    
    % estimate MLE parameters of gaussian for this
    [muhat, sigmahat] = normfit(data);
    x2_mle_params(1,ii) = muhat;
    x2_mle_params(2,ii) = sigmahat;
end

end

