%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%* 
%**************************************************************************

% A script to see the performace of a simple BN w/ a discrete and
% continuous node w/ the copula construction.  Here, we use the full joint
% distribution formula given by Song, Li, Yuan.  See Joint Regression 
% Analysis of Correlated Data Using Gaussian Copulas -- Biometrics 2009,
% Eq 9

clear;
clc;

rng(12345);

aa = 1; bb = 2;
D = 2;
dag = zeros(D,D);
dag(aa,bb) = 1;
discreteNodes = [aa];
nodeNames = {'A', 'B'};
discreteNodeNames = {'A'};
bntPath = '../bnt'; addpath(genpath(bntPath));

% generate the data
M = 1000;

copulaType = 'Frank';
alpha = 8;
U = copularnd(copulaType, alpha, M);

a_probs = [0.25 0.25 0.25 0.25];
a_dist = makedist('Multinomial','Probabilities',a_probs);

xContinuous = [normrnd(-2,0.3,M/2,1); normrnd(2,0.8,M/2,1)];
xContinuous = xContinuous(randperm(M),:);
[fContinous,xiContinuous] = emppdf(xContinuous,0);
FContinuous = empcdf(xContinuous,0);
continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous,0);
X_hybrid(:,1) = a_dist.icdf(U(:,1));
for ii=1:M
    X_hybrid(ii,2) = continuousDistInfo.icdf(U(ii,2));
end

% estimate the copula density
X_hybrid_continued = X_hybrid;
X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));
% generate pseudo-observations
% U_hybrid_continued = pobs(X_hybrid_continued, 'ecdf', 100); %% HCBN
U_hybrid_continued = pobs(X_hybrid_continued);        %% not HCBN

K = 100; h = 0.01;
u = linspace(0,1,K);
[U1,U2] = ndgrid(u);

c_est = empcopulapdf(fliplr(U_hybrid_continued), h, K, 'betak');        % do a fliplr to emulate how hcbn.m processes the data
c_actual = reshape(copulapdf(copulaType, [U1(:) U2(:)], alpha), K, K);

% integrate the discrete dimnesion (u1).  This equates to taking the
% partial derivative of C w.r.t. u2 (the continuous dimension)
C_est_discreteIntegrate = cumtrapz(u,c_est,2);      % because we did a fliplr, the discrete dimension is dim=2 now
C_actual_discreteIntegrate = cumtrapz(u,c_actual,1);

% now compute the expression for f(y1,y2), where y1 is discrete, and y2 is
% continuous
% The expression for f(y) is:
% f(y1,y2) = g2(y2)*[C(G1(y1),G2(y2)) - C(G1(y1-),G2(y2))], where C is the
% partial derivative of the copula w.r.t. u2, b/c u1 is discrete dimension

isdiscrete = 0;
[f_y2_est, xi_y2] = emppdf(X_hybrid(:,2),isdiscrete);
F_y2_est = empcdf(X_hybrid(:,2),isdiscrete);
disty2Est = rvEmpiricalInfo(xi_y2,f_y2_est,F_y2_est,isdiscrete);

isdiscrete = 1;
[f_y1_est, xi_y1] = emppdf(X_hybrid(:,1),isdiscrete);
F_y1_est = empcdf(X_hybrid(:,1),isdiscrete);
disty1Est = rvEmpiricalInfo(xi_y1,f_y1_est,F_y1_est,isdiscrete);

mteObj = mtebn(X_hybrid, discreteNodes, dag);
clgObj = clgbn(X_hybrid, discreteNodes, dag);
hcbnObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, K, h, dag); 

for y1=1:4
    
    f_y2_given_y1_copula_est = zeros(1,length(xiContinuous));
    f_y2_given_y1 = zeros(1,length(xiContinuous));
    f_y2_given_y1_KDE = zeros(1,length(xiContinuous));
    f_y2_given_y1_clg = zeros(1,length(xiContinuous));
    f_y2_given_y1_mte = zeros(1,length(xiContinuous));
    f_y2_given_y1_hcbn = zeros(1,length(xiContinuous));
    
    % estimate KDE, MTE, and CLG
    X_continuous_subset = [];
    for jj=1:M
        if(X_hybrid(jj,1)==y1)
            X_continuous_subset = [X_continuous_subset; X_hybrid(jj,2)];
        end
    end
    isdiscrete = 0;
    [f_kde,xi_kde] = emppdf(X_continuous_subset, isdiscrete);
    F_kde = empcdf(X_continuous_subset, isdiscrete);
    conditionalKDE = rvEmpiricalInfo(xi_kde, f_kde, F_kde, isdiscrete);
        
    for xi_idx=1:length(xi_y2)
        xi = xi_y2(xi_idx);
        
        % notice, we flip the arguments for u2_est and u1_est to be the
        % continuous first, then discrete, this is b/c this is how hcbn.m
        % estiamtes the copula (see how c_est was generated above).  
        u2_est = [disty2Est.cdf(xi) disty1Est.cdf(y1)];
        u1_est = [disty2Est.cdf(xi) disty1Est.cdf(y1-1)];
        
        u2_actual = [a_dist.cdf(y1) continuousDistInfo.cdf(xi)];
        u1_actual = [a_dist.cdf(y1-1) continuousDistInfo.cdf(xi)];
        
        f_y2_given_y1(xi_idx) = (continuousDistInfo.pdf(xi)*(empcopulaval(C_actual_discreteIntegrate,u2_actual) - ...
                                                               empcopulaval(C_actual_discreteIntegrate,u1_actual) ))/ a_dist.pdf(y1);
        f_y2_given_y1_copula_est(xi_idx) = (disty2Est.pdf(xi)*(empcopulaval(C_est_discreteIntegrate,u2_est) - ...
                                                          empcopulaval(C_est_discreteIntegrate,u1_est)))/disty1Est.pdf(y1);
        f_y2_given_y1_KDE(xi_idx) = conditionalKDE.pdf(xi);
        
        f_y2_given_y1_clg(xi_idx) = normpdf(xi, clgObj.bnParams{2}{y1}.Mean, clgObj.bnParams{2}{y1}.Covariance);
        f_y2_given_y1_mte(xi_idx) = mteObj.bnParams{2}{y1}.mte_info.pdf(xi);
        f_y2_given_y1_hcbn(xi_idx) = hcbnObj.computeMixedConditionalProbability_([y1 xi], [bb aa], bb);
    end
    
    f = figure(1);
    plot(xiContinuous,f_y2_given_y1, 'b*-', ...
         xiContinuous, f_y2_given_y1_copula_est, 'rp--', ...
         xiContinuous, f_y2_given_y1_KDE, ...
         xiContinuous, f_y2_given_y1_clg, ...
         xiContinuous, f_y2_given_y1_mte, ...
         xiContinuous, f_y2_given_y1_hcbn, 'kh--'); grid on;
    legend('Generative Model', 'Copula Estimate_{RANK}', 'KDE Estimate', 'CLG', 'MTE', 'HCBN_{ECDF}'); 
    title(sprintf('Y_1 = %d',y1));
    pause;
    clf(f);
    
end