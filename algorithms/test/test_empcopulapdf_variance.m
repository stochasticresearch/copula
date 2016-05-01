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

% Some experiments on copula density estimation and the variance associated
% with the density estimation

%% 2 discrete parents, 1 continuous child, test bias/var of continuous marg
clear;
clc;

load('mte_default_gaussian.mat');
MTE_DEFAULT_GAUSSIAN = mte_default_gaussian;

D = 3;
numMCSims = 25;
alpha = 5;
h = 0.05; 
K = 25;
M_vec = 100:100:500;

a_dist = makedist('Multinomial', 'probabilities', [0.1 0.2 0.3 0.05 0.05 0.15 0.05 0.1]);
b_dist = makedist('Multinomial', 'probabilities', [0.2 0.1 0.05 0.3 0.15 0.05 0.1 0.05]);
betapdf_a = 2;
betapdf_b = 4;
u = linspace(0,1,K);
domainC = u;   % beta-pdf is only between 0 and 1

[U1_3,U2_3,U3_3] = ndgrid(u);
[U1_2,U2_2] = ndgrid(u);
c_actual_all = reshape(gumbelcopulapdf([U1_3(:) U2_3(:) U3_3(:)], alpha),K,K,K);
c_actual_parents = reshape(gumbelcopulapdf([U1_2(:) U2_2(:)], alpha),K,K);
C_actual_all_discrete_integrate = cumtrapz(u,cumtrapz(u,c_actual_all,1),2);
C_actual_parents_discrete_integrate = cumtrapz(u,cumtrapz(u,c_actual_parents,1),2);

f_C_given_AB_actual_mat = zeros(length(a_dist.Probabilities), length(b_dist.Probabilities), length(domainC));
f_C_given_AB_estCopula_mat = zeros(numMCSims, length(a_dist.Probabilities), length(b_dist.Probabilities), length(domainC), length(M_vec));
f_C_given_AB_estCLG_mat = zeros(numMCSims, length(a_dist.Probabilities), length(b_dist.Probabilities), length(domainC), length(M_vec));
f_C_given_AB_estMTE_mat = zeros(numMCSims, length(a_dist.Probabilities), length(b_dist.Probabilities), length(domainC), length(M_vec));

for aVal=1:length(a_dist.Probabilities)
    for bVal=1:length(b_dist.Probabilities)
        cValIdx = 1;
        for cVal=domainC
            uuActual_ABC_1 = [a_dist.cdf(aVal) ...
                              b_dist.cdf(bVal) ...
                              betacdf(cVal, betapdf_a, betapdf_b)];
            uuActual_ABC_2 = [a_dist.cdf(aVal-1) ...
                              b_dist.cdf(bVal) ...
                              betacdf(cVal, betapdf_a, betapdf_b)];
            uuActual_ABC_3 = [a_dist.cdf(aVal) ...
                              b_dist.cdf(bVal-1) ...
                              betacdf(cVal, betapdf_a, betapdf_b)];
            uuActual_ABC_4 = [a_dist.cdf(aVal-1) ...
                              b_dist.cdf(bVal-1) ...
                              betacdf(cVal, betapdf_a, betapdf_b)];
            uuActual_AB_1 = uuActual_ABC_1(1:2);
            uuActual_AB_2 = uuActual_ABC_2(1:2);
            uuActual_AB_3 = uuActual_ABC_3(1:2);
            uuActual_AB_4 = uuActual_ABC_4(1:2);
            
            f_C = betapdf(cVal, betapdf_a, betapdf_b);
                    
            C_actual_ABC = empcopulaval(C_actual_all_discrete_integrate, uuActual_ABC_1) - ...
                           empcopulaval(C_actual_all_discrete_integrate, uuActual_ABC_2) - ...
                           empcopulaval(C_actual_all_discrete_integrate, uuActual_ABC_3) + ...
                           empcopulaval(C_actual_all_discrete_integrate, uuActual_ABC_4);
            C_actual_AB = empcopulaval(C_actual_parents_discrete_integrate, uuActual_AB_1) - ...
                          empcopulaval(C_actual_parents_discrete_integrate, uuActual_AB_2) - ...
                          empcopulaval(C_actual_parents_discrete_integrate, uuActual_AB_3) + ...
                          empcopulaval(C_actual_parents_discrete_integrate, uuActual_AB_4);
            
            f_C_given_AB_actual = f_C*C_actual_ABC/C_actual_AB;
            f_C_given_AB_actual_mat(aVal,bVal,cValIdx) = f_C_given_AB_actual;
            cValIdx = cValIdx + 1;
        end
    end
end
                    
mVecIdx = 1;
for M=M_vec
    for mcSimNum=1:numMCSims
        fprintf('Running MCSim=%d\n', mcSimNum);
        U = gumbelcopularnd(M, D, alpha);
        
        X = zeros(size(U));
        X(:,1) = a_dist.icdf(U(:,1));
        X(:,2) = b_dist.icdf(U(:,2));
        X(:,3) = betainv(U(:,3), betapdf_a, betapdf_b);
        
        % generate estimates of A, B, and C
        a_samps = random(a_dist, M, 1);
        b_samps = random(b_dist, M, 1);
        c_samps = betarnd(betapdf_a, betapdf_b, M, 1);
        
        isdiscrete = 1;
        [fAest,xAest] = emppdf(a_samps,isdiscrete);
        FAest = empcdf(a_samps,isdiscrete);
        distAEst = rvEmpiricalInfo(xAest,fAest,FAest);

        [fBest,xBest] = emppdf(b_samps,isdiscrete);
        FBest = empcdf(b_samps,isdiscrete);
        distBEst = rvEmpiricalInfo(xBest,fBest,FBest);

        isdiscrete = 0;
        [fCest,xCest] = emppdf(c_samps,isdiscrete); 
        FCest = empcdf(c_samps,isdiscrete);
        distCEst = rvEmpiricalInfo(xCest,fCest,FCest);
        
        % continue A and B
        X(:,1) = continueRv(X(:,1));
        X(:,2) = continueRv(X(:,2));
        
        % generate pseudo observations
        U_est = pseudoobs(X);
        
        % compute empirical copula
        c_est_all = empcopulapdf(U_est, h, K, 'betak');
        c_est_parents = empcopulapdf(U_est(:,1:2), h, K, 'betak');
        
        % integrate discrete part out to compute the PDF's
        C_est_all_discrete_integrate = cumtrapz(u,cumtrapz(u,c_est_all,1),2);
        C_est_parents_discrete_integrate = cumtrapz(u,cumtrapz(u,c_est_parents,1),2);
        
        % compute the marginal continuous distribution (C)
        for aVal=1:length(a_dist.Probabilities)
            for bVal=1:length(b_dist.Probabilities)
                
                % find the subset of values which have this combination of
                % parents
                combo = [aVal bVal];
                X_continuous_subset = [];
                for jj=1:M
                    if(isequal(combo, X(jj,1:2)))
                        X_continuous_subset = [X_continuous_subset; X(jj,3)];
                    end
                end
                
                % fit gaussian to the continuous subset
                if(isempty(X_continuous_subset))
                    gaussFit_mu = 0;
                    gaussFit_sigma = 1;
                else
                    [gaussFit_mu, gaussFit_sigma] = normfit(X_continuous_subset);
                end
                
                % fit MTE to the continuous subset
                if(isempty(X_continuous_subset))
                    mteFit = MTE_DEFAULT_GAUSSIAN;
                else
                    mteFit = estMteDensity(X_continuous_subset);
                end
                
                cValIdx = 1;
                for cVal=domainC
                    uuEst_ABC_1 = [distAEst.cdf(aVal) ...
                                   distBEst.cdf(bVal) ...
                                   distCEst.cdf(cVal)];
                    uuEst_ABC_2 = [distAEst.cdf(aVal-1) ...
                                   distBEst.cdf(bVal) ...
                                   distCEst.cdf(cVal)];
                    uuEst_ABC_3 = [distAEst.cdf(aVal) ...
                                   distBEst.cdf(bVal-1) ...
                                   distCEst.cdf(cVal)];
                    uuEst_ABC_4 = [distAEst.cdf(aVal-1) ...
                                   distBEst.cdf(bVal-1) ...
                                   distCEst.cdf(cVal)];
                    uuEst_AB_1 = uuEst_ABC_1(1:2);
                    uuEst_AB_2 = uuEst_ABC_2(1:2);
                    uuEst_AB_3 = uuEst_ABC_3(1:2);
                    uuEst_AB_4 = uuEst_ABC_4(1:2);
                    
                    f_C_est = distCEst.pdf(cVal);
                    
                    C_est_ABC = empcopulaval(C_est_all_discrete_integrate, uuEst_ABC_1, 'interpolate', 1) - ...
                                empcopulaval(C_est_all_discrete_integrate, uuEst_ABC_2, 'interpolate', 1) - ...
                                empcopulaval(C_est_all_discrete_integrate, uuEst_ABC_3, 'interpolate', 1) + ...
                                empcopulaval(C_est_all_discrete_integrate, uuEst_ABC_4, 'interpolate', 1);
                    C_est_AB = empcopulaval(C_est_parents_discrete_integrate, uuEst_AB_1, 'interpolate', 1) - ...
                               empcopulaval(C_est_parents_discrete_integrate, uuEst_AB_2, 'interpolate', 1) - ...
                               empcopulaval(C_est_parents_discrete_integrate, uuEst_AB_3, 'interpolate', 1) + ...
                               empcopulaval(C_est_parents_discrete_integrate, uuEst_AB_4, 'interpolate', 1);
                    
                    f_C_given_AB_est_copula = f_C_est*C_est_ABC/C_est_AB;
                    f_C_given_AB_est_clg = normpdf(cVal,gaussFit_mu,gaussFit_sigma);
                    f_C_given_AB_est_mte = mteFit.pdf(cVal);
                    
                    if(isnan(f_C_given_AB_est_copula))
                        1;  % interactive debugging :D
                        f_C_given_AB_est_copula = 0;       % TODO: is this correct?
                    end
                    
                    f_C_given_AB_estCopula_mat(mcSimNum, aVal, bVal, cValIdx, mVecIdx) = f_C_given_AB_est_copula;
                    f_C_given_AB_estCLG_mat(mcSimNum, aVal, bVal, cValIdx, mVecIdx) = f_C_given_AB_est_clg;
                    f_C_given_AB_estMTE_mat(mcSimNum, aVal, bVal, cValIdx, mVecIdx) = f_C_given_AB_est_mte;
                    
                    cValIdx = cValIdx + 1;
                end
            end
        end 
    end
    mVecIdx = mVecIdx + 1;
end

f_C_given_AB_actual_mat_repmat = repmat(f_C_given_AB_actual_mat,[ones(1,D) length(mVec)]);

f_C_given_AB_estCopula_mat_mcAvg = squeeze(mean(f_C_given_AB_estCopula_mat,1));
f_C_given_AB_estCLG_mat_mcAvg = squeeze(mean(f_C_given_AB_estCLG_mat,1));
f_C_given_AB_estMTE_mat_mcAvg = squeeze(mean(f_C_given_AB_estMTE_mat,1));

% calculate bias of these estimators for each possible parent combination
copulaEstBias = f_C_given_AB_actual_mat_repmat-f_C_given_AB_estCopula_mat_mcAvg;
clgEstBias = f_C_given_AB_actual_mat_repmat-f_C_given_AB_estCLG_mat_mcAvg;
mteEstBias = f_C_given_AB_actual_mat_repmat-f_C_given_AB_estMTE_mat_mcAvg;

% calculate variance of these estimators for each possible parent
% combination
copulaEstVar = squeeze(mean(f_C_given_AB_estCopula_mat.^2,1))-f_C_given_AB_estCopula_mat_mcAvg.^2;
clgEstVar = squeeze(mean(f_C_given_AB_estCLG_mat.^2,1))-f_C_given_AB_estCLG_mat_mcAvg.^2;
mteEstVar = squeeze(mean(f_C_given_AB_estMTE_mat.^2,1))-f_C_given_AB_estMTE_mat_mcAvg.^2;
