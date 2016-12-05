function [pval] = ktauhatpval(ktauhatVal, M, varargin)
%RSDMPVAL - computes the p-value of a given ktauhat measure and the number of
%samples upon which that ktauhat measure was calculated against independence
%hypothesis
% Inputs:
%  rsdmVal - the rsdm measure
%  M - the sample size used to compute this rsdm measure
%  varargin{1} - type of data for which this RSDM metric was computed,
%                continuous, discrete, hybrid1, or hybrid2.  continuous is
%                default
% values are hard-coded, look at the script rsdm_runPower which has a
% execution-cell inside it that generates these vectors!
        
muVecContinuous = [-0.0006 -0.0010 0.0017 0.0003 -0.0003 0.0007 -0.0000 -0.0006 0.0005 0.0002];
muVecHybrid1 = [-0.0012 -0.0009 0.0016 0.0005 -0.0001 0.0005 -0.0001 -0.0008 0.0008 0.0002];
muVecHybrid2 = [-0.0002 -0.0013 0.0008 -0.0001 -0.0003 0.0008 -0.0004 -0.0007 0.0005 0.0002];
muVecDiscrete = [-0.0010 -0.0012 0.0009 0.0000 -0.0000 0.0008 -0.0007 -0.0012 0.0009 0.0002];

sigmaVecContinuous = [0.0656 0.0484 0.0383 0.0331 0.0314 0.0261 0.0258 0.0231 0.0223 0.0217];
sigmaVecHybrid1 = [0.0699 0.0510 0.0403 0.0348 0.0324 0.0271 0.0267 0.0242 0.0230 0.0224];
sigmaVecHybrid2 = [0.0708 0.0514 0.0394 0.0344 0.0329 0.0273 0.0269 0.0239 0.0231 0.0224];
sigmaVecDiscrete = [0.0808 0.0605 0.0471 0.0414 0.0391 0.0326 0.0322 0.0290 0.0276 0.0269];

if(length(varargin)>=1)
    selector = varargin{1};
else
    selector = 'continuous';
end

if(strcmpi(selector, 'continuous'))
    muVec = muVecContinuous;
    sigmaVec = sigmaVecContinuous;
elseif(strcmpi(selector, 'hybrid1'))
    muVec = muVecHybrid1;
    sigmaVec = sigmaVecHybrid1;
elseif(strcmpi(selector, 'hybrid2'))
    muVec = muVecHybrid2;
    sigmaVec = sigmaVecHybrid2;
else    % assume discrete
    muVec = muVecDiscrete;
    sigmaVec = sigmaVecDiscrete;
end

M_vec = 100:100:1000;
% iterpolate for each value of k, mu, sigma for the given M
mu = interp1(M_vec, muVec, M, 'spline');
sigma = interp1(M_vec, sigmaVec, M, 'spline');

pval = 1-normcdf(ktauhatVal, k, mu, sigma);

end