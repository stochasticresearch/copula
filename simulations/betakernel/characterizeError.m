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

% Script which generates multivariate copulas of different types with discrete 
% marginal data, and then reconstructs the copula using beta kernels. 
% It then characterizes the error between the constructed copula and 
%  1.) the actual copula
%  2.) c(F(x1) ... f(Xn)) and f(x1 ... xn) / [ f(x1) * ... * f(xn) ]

clear;
clc;


% define all the parametrics we would like to test over
copulaTypes = {'Gaussian', 'T', 'Frank', 'Gumbel', 'Clayton'};
K_vec = 25:25:100;
M_vec = 250:250:1000;
N_vec = 2:4;
h_vec = [0.01 0.05 0.1 0.2];

% TODO: parametrize these also!
Rho_2D = [1 0.3; 0.3 1];
Rho_3D = [1 0.3 0.6; 0.3 1 0.7; 0.6 0.7 1];
Rho_4D = [1 0.3 0.6 0.4; 0.3 1 0.7 0.1; 0.6 0.7 1 0.2; 0.4 0.1 0.2 1];
nu = 1;
alpha = 6;
numECDFPts = 100;

% define the marginal distributions to be multinomial distributions, you
% should also parametrize this with different probabilities and for
% different number of discrete outcomes
discreteMarginalDist_1 = makedist('Multinomial','Probabilities',[0.2 0.2 0.2 0.2 0.2]);
discreteMarginalDist_2 = makedist('Multinomial','Probabilities',[0.5 0.1 0.1 0.1 0.2]);
discreteMarginalDist_3 = makedist('Multinomial','Probabilities',[0.9 0.01 0.01 0.04 0.04]);
discreteMarginalDist_4 = makedist('Multinomial','Probabilities',[0.2 0.3 0.4 0.05 0.05]);
discreteMarginalDistributions = {discreteMarginalDist_1, ...
    discreteMarginalDist_2, ...
    discreteMarginalDist_3, ...
    discreteMarginalDist_4};

e_dens_dobsMat = zeros(length(copulaTypes), ...
                 length(K_vec),...
                 length(M_vec), ...
                 length(N_vec),...
                 length(h_vec));
e_dens_bestcaseMat = zeros(size(e_dens_dobsMat));
e_dens_trueMat = zeros(size(e_dens_dobsMat));

e_probMat = zeros(size(e_dens_dobsMat));    % contains the error between 
                                            % c(F(x1) ... F(xn)) and
                                            % f(x1 ... xn) / [f(x1) * ... * f(xn) ]

cTypeIdx = 1; 
for copulaType=copulaTypes
    KvecIdx = 1;
    for K=K_vec
        MvecIdx = 1;
        for M=M_vec
            NvecIdx = 1;
            for N=N_vec
                hvecIdx = 1;
                for h=h_vec
                    u = linspace(0,1,K);
                    if(N==2)
                        Rho = Rho_2D;
                        [U1,U2] = ndgrid(u);
                        copPdf_U = [U1(:) U2(:)];
                        reshapeVec = [K,K];
                    elseif(N==3)
                        Rho = Rho_3D;
                        [U1,U2,U3] = ndgrid(u);
                        copPdf_U = [U1(:) U2(:) U3(:)];
                        reshapeVec = [K,K,K];
                    elseif(N==4)
                        Rho = Rho_4D;
                        [U1,U2,U3,U4] = ndgrid(u);
                        copPdf_U = [U1(:) U2(:) U3(:) U4(:)];
                        reshapeVec = [K,K,K,K];
                    end
                    % generate the copula random variates to set the
                    % dependency, and the actual copula densities
                    if(strcmp(copulaType, 'Gaussian'))
                        U = copularnd(copulaType, Rho, M);
                        c = copulapdf(copulaType, copPdf_U, Rho);
                    elseif(strcmp(copulaType, 'T'))
                        U = copularnd(copulaType, Rho, nu, M);
                        c = copulapdf(copulaType, copPdf_U, Rho, nu);
                    elseif(strcmp(copulaType, 'Frank'))
                        U = frankcopularnd(M, N, alpha);
                        c = frankcopulapdf(copPdf_U, alpha);
                    elseif(strcmp(copulaType, 'Gumbel'))
                        U = gumbelcopularnd(M, N, alpha);
                        c = gumbelcopulapdf(copPdf_U, alpha);
                    elseif(strcmp(copulaType, 'Clayton'))
                        U = claytoncopularnd(M, N, alpha);
                        c = claytoncopulapdf(copPdf_U, alpha);
                    else
                        error('Unknown copula type!');
                    end
                    c = reshape(c, reshapeVec);
                    
                    % apply icdf function to generate discrete random
                    % variates w/ the defined dependency structure
                    X = zeros(M,N);
                    for nn=1:N
                        X(:,nn) = icdf(discreteMarginalDistributions{nn}, U(:,nn));
                    end
                    
                    % continue the discrete observations
                    X_continued = continueRv( X );
                    
                    % generate pseudo-observations from the empirical
                    % marginal distribution functions for the continued
                    % random variables, and estimate the empirical copula
                    % density from those
                    U_pseudoObs = pobs(X_continued);
                    
                    % estimate the copula density from the continued 
                    % observations - call this c_hat_dobs
                    c_hat_dobs = empcopulapdf(U_pseudoObs, h, K, 'betak');
                    
                    % compute error between estimated copula (c_hat_dobs)
                    % density and actual copula density c, call this e_dobs
                    e_dobs = hyperFunctionError(c,c_hat_dobs);
                    
                    % estimate the copula density w/ the copula random 
                    % variates directly - call this c_hat
                    c_hat = empcopulapdf(U, h, K, 'betak');
                    
                    % compute error between c_hat and c, call this
                    % e_bestcase
                    e_bestcase = hyperFunctionError(c,c_hat);
                    
                    % compute error between c_hat and c_hat_dobs, call this
                    % e_true
                    e_true = hyperFunctionError(c_hat,c_hat_dobs);
                    
                    e_dens_dobsMat(cTypeIdx,KvecIdx,MvecIdx,NvecIdx,hvecIdx) = e_dobs;
                    e_dens_bestcaseMat(cTypeIdx,KvecIdx,MvecIdx,NvecIdx,hvecIdx) = e_bestcase;
                    e_dens_trueMat(cTypeIdx,KvecIdx,MvecIdx,NvecIdx,hvecIdx) = e_true;
                    
                    % compute the error between  c(F(x1) ... F(xn)) and
                    % f(x1 ... xn) / [f(x1) * ... * f(xn) ], we use
                    % c_hat_dobs for the copula calculation
                    [prob_using_density, combos] = computeEmpiricalDiscretProb( X );
                    prob_using_c_hat_obs = zeros(size(prob_using_density));
                    U_psuedoObs = zeros(1,length(combo));
                    for comboIdx=1:size(combos,1)
                        combo = combos(comboIdx,:);
                        
                        % change the combo to pseudo-observations
                        for jj=1:length(combo)
                            U_pseudoObs(jj) = cdf(discreteMarginalDistributions{jj},combo(jj));
                        end
                        
                        % query the copula value
                        prob_using_c_hat_obs(comboIdx) = empcopulaval(c_hat_dobs, U_pseudoObs);
                    end
                    
                    hvecIdx = hvecIdx + 1;
                end
                NvecIdx = NvecIdx + 1;
            end
            MvecIdx = MvecIdx + 1;
        end
        KvecIdx = KvecIdx + 1;
    end
    cTypeIdx = cTypeIdx + 1;
end