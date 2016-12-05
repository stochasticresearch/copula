function [pval] = rsdmpval(rsdmVal, M, varargin)
%RSDMPVAL - computes the p-value of a given RSDM measure and the number of
%samples upon which that RSDM measure was calculated against independence
%hypothesis
% Inputs:
%  rsdmVal - the rsdm measure
%  M - the sample size used to compute this rsdm measure
%  varargin{1} - type of data for which this RSDM metric was computed,
%                continuous, discrete, hybrid1, or hybrid2.  continuous is
%                default
% values are hard-coded, look at the script rsdm_runPower which has a
% execution-cell inside it that generates these vectors!
        
alphaVecContinuous = [9.8809 9.2044 9.9171 9.4418 9.0166 9.5859 9.3859 10.0273 8.8882 9.8018];
alphaVecHybrid1 = [6.7514 8.0427 8.5287 7.9632 7.7020 8.4284 8.6409 8.6185 8.1371 8.5404];
alphaVecHybrid2 = [7.1260 7.5886 7.6971 8.5812 8.0322 7.9201 8.3371 8.2508 8.0096 8.2058];
alphaVecDiscrete = [6.7738 6.7823 7.6465 7.4081 7.0405 7.3359 7.4724 7.7517 7.4365 7.8571];

betaVecContinuous = [64.9703 95.1447 132.4876 143.5421 155.5770 184.3349 196.4452 225.4885 210.0091 245.0919];
betaVecHybrid1 = [38.5076 76.5102 105.7296 113.3117 126.6973 150.3545 170.6229 182.6477 181.6230 202.5102];
betaVecHybrid2 = [40.8525 72.2452 95.5144 121.4174 131.7344 144.2803 165.9454 174.9780 181.2124 196.5054];
betaVecDiscrete = [39.4650 58.7468 83.6098 92.3611 102.3834 116.8612 130.9255 142.3449 145.5351 163.0671];

if(length(varargin)>=1)
    selector = varargin{1};
else
    selector = 'continuous';
end

if(strcmpi(selector, 'continuous'))
    alphaVec = alphaVecContinuous;
    betaVec = betaVecContinuous;
elseif(strcmpi(selector, 'hybrid1'))
    alphaVec = alphaVecHybrid1;
    betaVec = betaVecHybrid1;
elseif(strcmpi(selector, 'hybrid2'))
    alphaVec = alphaVecHybrid2;
    betaVec = betaVecHybrid2;
else    % assume discrete
    alphaVec = alphaVecDiscrete;
    betaVec = betaVecDiscrete;
end

M_vec = 100:100:1000;
% iterpolate for each value of k, mu, sigma for the given M
alpha = interp1(M_vec, alphaVec, M, 'spline');
beta = interp1(M_vec, betaVec, M, 'spline');

pval = 1-betacdf(rsdmVal, k, alpha, beta);

end