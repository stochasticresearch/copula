function [pval] = rsdmpval(rsdmVal, M)
%RSDMPVAL - computes the p-value of a given RSDM measure and the number of
%samples upon which that RSDM measure was calculated
% Inputs:
%  rsdmVal - the rsdm measure
%  M - the sample size used to compute this rsdm measure

% values are hard-coded, look at the 3rd and 4th execution cells to see how
% these vectors were generated
kVec = [-0.1294, -0.1382, -0.1167, -0.1524, -0.1166, -0.0958, -0.1346, ...
        -0.1260, -0.1729, -0.1102];
muVec = [0.1128, 0.0770, 0.0633, 0.0541, 0.0482, 0.0431, 0.0400, ...
         0.0368, 0.0358, 0.0330];
sigmaVec = [0.0357, 0.0242, 0.0190, 0.0178, 0.0151, 0.0131, 0.0127, 0.0117, ...
            0.0116, 0.0104];
M_vec = 100:100:1000;

% iterpolate for each value of k, mu, sigma for the given M

k = interp1(M_vec, kVec, M, 'spline');
mu = interp1(M_vec, muVec, M, 'spline');
sigma = interp1(M_vec, sigmaVec, M, 'spline');

pval = 1-gevcdf(rsdmVal, k, sigma, mu);

end