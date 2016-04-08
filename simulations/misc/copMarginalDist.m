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
continuousDistInfo = rvEmpiricalInfo(xi2,f2,F2);

numPtsToTest = [1e2 1e3 1e4 1e5];
for numPtToTest=numPtsToTest
    u = unifrnd(0,1,numPtToTest,1);
    f_xform = zeros(1,numPtToTest);
    for ii=1:numPtToTest
        f_xform(ii) = continuousDistInfo.invDistribution(u(ii));
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
% skewed? -- in the HCBN context (all 2D tests)

clear;
clc;

% stuff to initialize the HCBN OBJ
aa = 1; bb = 2;
D = 2;
dag = zeros(D,D);
dag(aa,bb) = 1;
discreteNodes = [aa];
nodeNames = {'A', 'B'};
bntPath = '../bnt'; addpath(genpath(bntPath));
discreteNodeNames = {'A'};

% parametrization variables
mVec = 500:500:2000;    % TODO: more finer grained after you like results
copulaTypeVec = {'Frank', 'Gumbel', 'Clayton'};
alphaVec = 1:3:10;
x2DistTypeVec = {'Multimodal', 'Uniform', 'Gaussian', 'ThickTailed'};
numKLDivCalculated = 7;
probs = [0.15 0.35 0.4 0.1];        % TODO: parametrize this also
a_dist = makedist('Multinomial','Probabilities',probs);

klDivMat = zeros(length(copulaTypeVec),...
                 length(x2DistTypeVec),...
                 length(alphaVec),...
                 length(mVec),...
                 numKLDivCalculated, ...
                 length(probs));

numMCSims = 25;
klDivMCMat = zeros(numKLDivCalculated,length(probs),numMCSims);

fid = fopen('/home/kiran/ownCloud/PhD/sim_results/copMarginalDist.txt', 'a');
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
numTotalLoops = length(copulaTypeVec)*length(x2DistTypeVec)*length(alphaVec)*length(mVec)*numMCSims;
progressIdx = 1;
for copulaTypeVecIdx=1:length(copulaTypeVec)
    for x2DistTypeVecIdx=1:length(x2DistTypeVec)
        for alphaVecIdx=1:length(alphaVec)
            for mVecIdx=1:length(mVec)
                copulaType = copulaTypeVec{copulaTypeVecIdx};
                x2DistType = x2DistTypeVec{x2DistTypeVecIdx};
                alpha = alphaVec(alphaVecIdx);
                M = mVec(mVecIdx);
                
                %%%%%%%%%%% MAIN SIMULATION CODE %%%%%%%%%%
                for mcSimNum=1:numMCSims
                    progressAmt = progressIdx/numTotalLoops*100;
                    progressStr = sprintf('copulaType=%s x2DistType=%s alpha=%d M=%d MC Sim# = %d || Progress=%0.04f\n', ...
                        copulaType, x2DistType, alpha, M, mcSimNum, progressAmt);
                    dispstat(progressStr,'timestamp');
                    fprintf(fid, progressStr);
                    
                    U = copularnd(copulaType, alpha, M);
                    X_hybrid = zeros(M,2);

                    % make both X1 and X2 multimodal distributions
                    if(strcmp(x2DistType, 'Multimodal'))
                        x2 = [normrnd(-2,0.3,M/2,1); normrnd(2,0.8,M/2,1)];
                    elseif(strcmp(x2DistType, 'Uniform'))
                        x2 = unifrnd(-2,2,M,1);
                    elseif(strcmp(x2DistType, 'UnimodalSkewed'))
                        x2 = betarnd(2,5,M,1);
                    elseif(strcmp(x2DistType, 'Gaussian'))
                        x2 = normrnd(2,0.5,M,1);
                    elseif(strcmp(x2DistType, 'ThickTailed'))
                        x2 = trnd(1, M, 1);
                    else
                        error('Unknown X2 Dist Type!');
                    end
                    x2 = x2(randperm(M),:);     % permute for evenness of samples

                    [f2,xi2] = emppdf(x2,0);
                    F2 = empcdf(x2,0);
                    continuousDistInfo = rvEmpiricalInfo(xi2,f2,F2);
                    X_hybrid(:,1) = a_dist.icdf(U(:,1));
                    for ii=1:M
                        X_hybrid(ii,2) = continuousDistInfo.invDistribution(U(ii,2));
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
                    % U_hybrid_continued = pseudoobs(X_hybrid_continued, 'ecdf', 100); %% HCBN  implements this version, so we compare here w/ non-ecdf
                    U_hybrid_continued = pseudoobs(X_hybrid_continued);
                    h = 0.05; K = 100; c_est = empcopulapdf(U_hybrid_continued, h, K, 'betak');

                    % c_estTest = empcopulapdf(fliplr(U_hybrid_continued), h, K, 'betak');
                    % surf( (c_estTest-hcbnObj.copulaFamilies{2}.c).^2 ); pause;

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
                        [f_kde,xi_kde] = emppdf(X_continuous_subset, isdiscrete);
                        F_kde = empcdf(X_continuous_subset, isdiscrete);
                        conditionalKDE = rvEmpiricalInfo(xi_kde, f_kde, F_kde);

                        fx2_givenx1_copula = zeros(1,length(xi2));
                        fx2_givenx1_copulaestf2est = zeros(1,length(xi2));
                        fx2_givenx1_copulaestf2Actual = zeros(1,length(xi2));
                        fx2_givenx1_copulaActualf2est = zeros(1,length(xi2));
                        fx2_givenx1_hcbn = zeros(1,length(xi2));
                        fx2_givenx1_clg = zeros(1,length(xi2));
                        fx2_givenx1_mte = zeros(1,length(xi2));
                        fx2_givenx1_conditionalKDE = zeros(1,length(xi2));
                        % compute the conditional density w/ the formula
                        % f(x2|x1) = c(F(x1),F(x2))*f(x2)
                        for ii=1:length(xi2)
                            xi2_val = xi2(ii);
                            uuGenerative = [a_dist.cdf(x1_discrete_conditional) continuousDistInfo.queryDistribution(xi2_val)];
                            uuGenerative = fixU(uuGenerative);

                            fx2_givenx1_copula(ii) = copulapdf(copulaType, uuGenerative, alpha)*continuousDistInfo.queryDensity(xi2_val);
                            uuEst = [discreteEstDistInfo.queryDistribution(x1_discrete_conditional) ...
                                     continuousEstDistInfo.queryDistribution(xi2_val)];
                            uuEst = fixU(uuEst);
                            fx2_givenx1_copulaestf2est(ii) = empcopula_val(c_est, uuEst)*continuousEstDistInfo.queryDensity(xi2_val);
                            fx2_givenx1_copulaestf2Actual(ii) = empcopula_val(c_est, uuEst)*continuousDistInfo.queryDensity(xi2_val);
                            fx2_givenx1_copulaActualf2est(ii) = copulapdf(copulaType, uuGenerative, alpha)*continuousEstDistInfo.queryDensity(xi2_val);

                            uuHcbn = [hcbnObj.empInfo{1}.distribution(x1_discrete_conditional) ...
                                      hcbnObj.empInfo{2}.queryDistribution(xi2_val)];
                            uuHcbn = fixU(uuHcbn);
                            c_hcbn = empcopula_val(hcbnObj.copulaFamilies{2}.c, fliplr(uuHcbn));        % I think this should be a fliplr b/c copula
                                                                                                        % was estimated w/ [u2 u1]
                            f_x2_hcbn = hcbnObj.empInfo{2}.queryDensity(xi2_val);
                            fx2_givenx1_hcbn(ii) = c_hcbn*f_x2_hcbn;       
                            fx2_givenx1_clg(ii) = normpdf(xi2_val, clgObj.bnParams{2}{x1_discrete_conditional}.Mean, clgObj.bnParams{2}{x1_discrete_conditional}.Covariance);
                            fx2_givenx1_mte(ii) = mteObj.bnParams{2}{x1_discrete_conditional}.mte_info.queryDensity(xi2_val);
                            fx2_givenx1_conditionalKDE(ii) = conditionalKDE.queryDensity(xi2_val);
                        end

                        % compute kl_divergences
                        div_cEstf2Est = kldivergence(fx2_givenx1_copula, fx2_givenx1_copulaestf2est, xi2);
                        div_cEstf2Actual = kldivergence(fx2_givenx1_copula, fx2_givenx1_copulaestf2Actual, xi2);
                        div_cActualf2Est = kldivergence(fx2_givenx1_copula, fx2_givenx1_copulaActualf2est, xi2);
                        div_hcbn = kldivergence(fx2_givenx1_copula, fx2_givenx1_hcbn, xi2);
                        div_kde = kldivergence(fx2_givenx1_copula, fx2_givenx1_conditionalKDE, xi2);
                        div_mte = kldivergence(fx2_givenx1_copula, fx2_givenx1_mte, xi2);
                        div_clg = kldivergence(fx2_givenx1_copula, fx2_givenx1_clg, xi2);

                        klDivMCMat(1,x1_discrete_conditional,mcSimNum) = div_cEstf2Est;
                        klDivMCMat(2,x1_discrete_conditional,mcSimNum) = div_cEstf2Actual;
                        klDivMCMat(3,x1_discrete_conditional,mcSimNum) = div_cActualf2Est;
                        klDivMCMat(4,x1_discrete_conditional,mcSimNum) = div_hcbn;
                        klDivMCMat(5,x1_discrete_conditional,mcSimNum) = div_kde;
                        klDivMCMat(6,x1_discrete_conditional,mcSimNum) = div_mte;
                        klDivMCMat(7,x1_discrete_conditional,mcSimNum) = div_clg;

                %         % plot the actual vs the copula version and compare the differences
                %         fig1 = figure(1);
                %         plot(xi2, fx2_givenx1_copula, 'b*-', ...
                %              xi2, fx2_givenx1_copulaestf2est, ...
                %              xi2, fx2_givenx1_copulaestf2Actual, ...
                %              xi2, fx2_givenx1_copulaActualf2est, ...
                %              xi2, fx2_givenx1_hcbn, 'o-',...
                %              xi2, fx2_givenx1_conditionalKDE, ...
                %              xi2, fx2_givenx1_mte, ...
                %              xi2, fx2_givenx1_clg); 
                %         grid on; title(sprintf('x_1=%d',x1_discrete_conditional));
                %         h_legend = legend('$c*f(x_2)$', ...
                %             strcat('$\hat{c}_{RANK}*\hat{f}(x_2)$', sprintf(' KLdiv=%0.02f', div_cEstf2Est)), ...
                %             strcat('$\hat{c}_{RANK}*f(x_2)$', sprintf(' KLdiv=%0.02f', div_cEstf2Actual)), ...
                %             strcat('$c*\hat{f}(x_2)$', sprintf(' KLdiv=%0.02f', div_cActualf2Est)), ...
                %             strcat('$\hat{c}_{HCBN-ECDF}$', sprintf(' KLdiv=%0.02f', div_hcbn)), ...
                %             strcat('KDE', sprintf(' KLdiv=%0.02f', div_kde)), ...
                %             strcat('MTE', sprintf(' KLdiv=%0.02f', div_mte)), ...
                %             strcat('CLG', sprintf(' KLdiv=%0.02f', div_clg)) );
                %         set(h_legend,'FontSize',10);
                %         set(h_legend,'Interpreter','latex')
                %         pause(1);
                %         clf(fig1);
                    end
                    progressIdx = progressIdx + 1;
                end
                meanKLDivMCMat = mean(klDivMCMat,3);
                for x1_discrete_conditional=x1_discrete_conditional_vec
                    klDivMat(copulaTypeVecIdx, ...
                             x2DistTypeVecIdx,...
                             alphaVecIdx,...
                             mVecIdx,...
                             1, x1_discrete_conditional) = meanKLDivMCMat(1,x1_discrete_conditional);
                    klDivMat(copulaTypeVecIdx, ...
                             x2DistTypeVecIdx,...
                             alphaVecIdx,...
                             mVecIdx,...
                             2, x1_discrete_conditional) = meanKLDivMCMat(2,x1_discrete_conditional);
                    klDivMat(copulaTypeVecIdx, ...
                             x2DistTypeVecIdx,...
                             alphaVecIdx,...
                             mVecIdx,...
                             3, x1_discrete_conditional) = meanKLDivMCMat(3,x1_discrete_conditional);
                    klDivMat(copulaTypeVecIdx, ...
                             x2DistTypeVecIdx,...
                             alphaVecIdx,...
                             mVecIdx,...
                             4, x1_discrete_conditional) = meanKLDivMCMat(4,x1_discrete_conditional);
                    klDivMat(copulaTypeVecIdx, ...
                             x2DistTypeVecIdx,...
                             alphaVecIdx,...
                             mVecIdx,...
                             5, x1_discrete_conditional) = meanKLDivMCMat(5,x1_discrete_conditional);
                    klDivMat(copulaTypeVecIdx, ...
                             x2DistTypeVecIdx,...
                             alphaVecIdx,...
                             mVecIdx,...
                             6, x1_discrete_conditional) = meanKLDivMCMat(6,x1_discrete_conditional);
                    klDivMat(copulaTypeVecIdx, ...
                             x2DistTypeVecIdx,...
                             alphaVecIdx,...
                             mVecIdx,...
                             7, x1_discrete_conditional) = meanKLDivMCMat(7,x1_discrete_conditional);
                end
                %%%%%%%%%%% END OF MAIN SIMULATION CODE %%%%%%%%%%
            end
        end
    end
end

dispstat('Finished.','keepprev');
fclose(fid);
save('/home/kiran/ownCloud/PhD/sim_results/klDivMat.mat', 'klDivMat');