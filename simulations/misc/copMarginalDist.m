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

%% Understand the effect of INV Probability transform parametric w/ sample size

M = 2000;
x2 = [normrnd(-2,0.1,M/2,1); normrnd(2,0.8,M/2,1)];
x2 = x2(randperm(M),:);

[f2,xi2] = emppdf(x2,0);
F2 = empcdf(x2,0);
continuousMarg2 = rvEmpiricalInfo(xi2,f2,F2);

numPtsToTest = [1e2 1e3 1e4 1e5];
for numPtToTest=numPtsToTest
    u = unifrnd(0,1,numPtToTest,1);
    f_xform = zeros(1,numPtToTest);
    for ii=1:numPtToTest
        f_xform(ii) = continuousMarg2.invDistribution(u(ii));
    end
    [femp,xemp] = emppdf(f_xform, 0);
    fig1 = figure(1);
    plot(xi2,f2,xemp,femp); hold on; 
    isequal(xi2,xemp)
    title(sprintf('num pts=%d', numPtToTest));
    legend('True Marginal', 'Approximated Marginal w/ INV CDF');
    grid on;
    pause;
    clf(fig1);
end
close(fig1);

%%

% A script to understand copula conditional densities and why they may be
% skewed? -- in the HCBN context

clear;
clc;

% stuff to initialize the HCBN OBJ
aa = 1; bb = 2;
D =2;
dag = zeros(D,D);
dag(aa,bb) = 1;
discreteNodes = [aa];
nodeNames = {'A', 'B'};
bntPath = '../bnt'; addpath(genpath(bntPath));
discreteNodeNames = {'A'};

M = 2000;
D = 2; alpha = 10; 
U = copularnd('Gumbel', alpha, M);
X_hybrid = zeros(M,2);
probs = [0.15 0.35 0.4 0.1];
a_dist = makedist('Multinomial','Probabilities',probs);

% make both X1 and X2 multimodal distributions
x2 = [normrnd(-2,0.3,M/2,1); normrnd(2,0.8,M/2,1)];
x2 = x2(randperm(M),:);

[f2,xi2] = emppdf(x2,0);
F2 = empcdf(x2,0);
continuousMarg2 = rvEmpiricalInfo(xi2,f2,F2);
X_hybrid(:,1) = a_dist.icdf(U(:,1));
for ii=1:M
    X_hybrid(ii,2) = continuousMarg2.invDistribution(U(ii,2));
end

[f1est,x1est] = emppdf(X_hybrid(:,1),1);
F1est = empcdf(X_hybrid(:,1),1);
discreteEstDistInfo = rvEmpiricalInfo(x1est,f1est,F1est);

[f2est,x2est] = emppdf(X_hybrid(:,2),0); F2est = empcdf(X_hybrid(:,2),0);
continuousEstDistInfo = rvEmpiricalInfo(x2est,f2est,F2est);

hcbnObj = hcbn(bntPath, X_hybrid, nodeNames, discreteNodeNames, dag); 
mteObj = mte(X_hybrid, discreteNodes, dag);
clgObj = clg(X_hybrid, discreteNodes, dag);

X_hybrid_continued = X_hybrid;
X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));
% generate pseudo-observations
U_hybrid_continued = pseudoobs(X_hybrid_continued, 'ecdf', 100);
h = 0.05; K = 100; c_est = empcopulapdf(U_hybrid_continued, h, K, 'betak');

% fig1 = figure(1); 
% [uU1, uU2] = ndgrid(linspace(0,1,K));
% c_actual = reshape(frankcopulapdf([uU1(:) uU2(:)],alpha),K,K);
% subplot(1,2,1); surf(uU1, uU2, c_est); grid on; title('Estimated Copula');
% subplot(1,2,2); surf(uU1, uU2, c_actual); grid on; title('Actual Copula');
% pause;

% plot some conditional distributions from the known copula and known
% marginal distributions
x1_discrete_conditional_vec = 1:1:length(probs);

for x1_discrete_conditional=x1_discrete_conditional_vec
    
    X_continuous_subset = [];
    for jj=1:M
        if(X_hybrid(jj,1)==x1_discrete_conditional)
            X_continuous_subset = [X_continuous_subset; X_hybrid(jj,2)];
        end
    end
    % estimate the empirical PDF/CDF for this
    isdiscrete = 0;
    [f,xi] = emppdf(X_continuous_subset, isdiscrete);
    F = empcdf(X_continuous_subset, isdiscrete);
    conditionalKDE = rvEmpiricalInfo(xi, f, F);

    fx2_givenx1_copula = zeros(1,length(xi2));
    fx2_givenx1_copulaest = zeros(1,length(xi2));
    fx2_givenx1_hcbn = zeros(1,length(xi2));
    fx2_givenx1_clg = zeros(1,length(xi2));
    fx2_givenx1_mte = mteObj.bnParams{2}{x1_discrete_conditional}.mte_info.density;
    % compute the conditional density w/ the formula
    % f(x2|x1) = c(F(x1),F(x2))*f(x2)
    for ii=1:length(xi2)
        xi_val = xi2(ii);
        uuGenerative = [a_dist.cdf(x1_discrete_conditional) continuousMarg2.queryDistribution(xi_val)];
        uuGenerative = fixU(uuGenerative);
        
        fx2_givenx1_copula(ii) = copulapdf('Frank', uuGenerative, alpha)*f2(ii);
        uuEst = [discreteEstDistInfo.queryDistribution(x1_discrete_conditional) ...
                 continuousEstDistInfo.queryDistribution(xi_val)];
        uuEst = fixU(uuEst);
        fx2_givenx1_copulaest(ii) = empcopula_val(c_est, uuEst)*f2est(ii);
        
        uuHcbn = [hcbnObj.empInfo{1}.distribution(x1_discrete_conditional) ...
                  hcbnObj.empInfo{2}.queryDistribution(xi_val)];
        uuHcbn = fixU(uuHcbn);
        c_hcbn = empcopula_val(hcbnObj.copulaFamilies{2}.c, fliplr(uuHcbn));        % I think this should be a fliplr
        f_x2_hcbn = hcbnObj.empInfo{2}.density(ii);
        fx2_givenx1_hcbn(ii) = c_hcbn*f_x2_hcbn;
        
        fx2_givenx1_clg(ii) = normpdf(xi_val, clgObj.bnParams{2}{x1_discrete_conditional}.Mean, clgObj.bnParams{2}{x1_discrete_conditional}.Mean);
    end
    
    % plot the actual vs the copula version and compare the differences
    fig2 = figure(2);
    plot(xi2, fx2_givenx1_copula, ...
        xi2, fx2_givenx1_copulaest, ...
        xi2, fx2_givenx1_hcbn, ...
        xi, conditionalKDE.density, ...
        mteObj.bnParams{2}{x1_discrete_conditional}.mte_info.domain, fx2_givenx1_mte, ...
        xi2, fx2_givenx1_clg); 
    grid on; title(sprintf('Discrete << x1=%d',x1_discrete_conditional));
    legend('Generative Model PDF', 'Estimated Copula', 'HCBN Copula Estimate', ...
        'KDE Conditional', 'MTE', 'CLG');
    pause;
    clf(fig2);
end

close all;