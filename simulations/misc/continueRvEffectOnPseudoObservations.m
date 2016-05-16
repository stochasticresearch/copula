% A script which shows the effect of continuing a RV on
% pseudo-observations, which are then used to estimate the copula density

%% experiments on pseudo-observations of marginal distributions
M = 250;
X1 = normrnd(0,1,M,1);
X2 = exprnd(2, M, 1);
X3 = betarnd(2,3,M,1);
X4 = trnd(5, M, 1);

U1 = pseudoobs(X1);
U2 = pseudoobs(X2);
U3 = pseudoobs(X3);
U4 = pseudoobs(X4);
subplot(2,2,1); ksdensity(X1); hold on; ksdensity(U1); legend('X1', 'U1'); grid on;
subplot(2,2,2); ksdensity(X2); hold on; ksdensity(U2); legend('X1', 'U1'); grid on;
subplot(2,2,3); ksdensity(X3); hold on; ksdensity(U3); legend('X1', 'U1'); grid on;
subplot(2,2,4); ksdensity(X4); hold on; ksdensity(U4); legend('X1', 'U1'); grid on;

%% D=2
clear;
clc;

M = 1000;
alpha = 3;
U = copularnd('Gumbel', alpha, M);

% prob = [0.5 0.5];
% prob = [0.5 0.15 0.15 0.05 0.05 0.05 0.025 0.025];
% prob = .1*ones(1,10);
% prob = [0.25 0.25 0.25 0.25];
prob = [0.6 0.25 0.1 0.05];
a_dist = makedist('Multinomial','Probabilities',prob);
X_hybrid(:,1) = a_dist.icdf(U(:,1));
X_hybrid(:,2) = norminv(U(:,2), 0, 1);

X_hybrid_continued = X_hybrid;
X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));
% 
% for mm=1:M
%     if(X_hybrid_continued(mm,1)==1 || X_hybrid_continued(mm,1)==3)
%         X_hybrid_continued(mm,1) = X_hybrid_continued(mm,1)+rand()-1;
%     end
% end

U_hybrid = pseudoobs(X_hybrid);
U_hybrid_continued = pseudoobs(X_hybrid_continued);
U_hybrid_corrected = pseudoobs(X_hybrid_continued, 'correction', 50);

tau1 = corr(U_hybrid(:,1),U_hybrid(:,2),'type','spearman');
tau2 = corr(U_hybrid_continued(:,1),U_hybrid_continued(:,2),'type','spearman');

xx = linspace(0,1,M);
tau1_rescaled = tau1/(12*sqrt(prod(var(U_hybrid))));
tau2_rescaled = tau2/(12*sqrt(prod(var(U_hybrid_continued))));
yy1 = tau1_rescaled*xx;
yy2 = tau2_rescaled*xx;
yy3 = tau2*xx;

subplot(3,2,1);
scatter(X_hybrid(:,1),X_hybrid(:,2)); grid on; title('X Hybrid');
subplot(3,2,2);
scatter(X_hybrid_continued(:,1),X_hybrid_continued(:,2)); grid on; title('X Hybrid Continued');
subplot(3,2,3);
scatter(U_hybrid(:,1),U_hybrid(:,2)); grid on; title(sprintf('U Hybrid - tau=%0.04f', tau1_rescaled));
hold on;
scatter(U(:,1), U(:,2), 'r')
plot(xx,yy1, 'g+')
plot(xx,yy2, 'k+')
plot(xx,yy3, 'm--')
subplot(3,2,4);
scatter(U_hybrid_continued(:,1),U_hybrid_continued(:,2)); grid on; title('U Hybrid Continued');
subplot(3,2,6);
scatter(U_hybrid_corrected(:,1),U_hybrid_corrected(:,2)); 
axis([0 1 0 1]);
grid on; title('U Hybrid Continued & Fixed');

%% simulation for testing difference in tau's as a metric for determining
%  if we want to do this idea
clear;
clc;

copulaTypeVec = {'Frank', 'Gumbel', 'Clayton', 'Gaussian'};
alphaVec = 1:20;
rhoVec = 0:.05:.95;
M = 1000;
continuousDistTypeVec = {'Multimodal', 'Uniform', 'Gaussian', 'ThickTailed'};
discreteDistCfgVec = [1 2 3];
sRhoRescaledVec = zeros(length(copulaTypeVec), length(discreteDistCfgVec), length(continuousDistTypeVec), length(alphaVec));
sRhoActualVec = zeros(length(copulaTypeVec), length(alphaVec));
numMCSims = 25;

for copulaTypeVecIdx=1:length(copulaTypeVec)
    for alphaVecIdx=1:length(alphaVec)
        for discreteDistCfg=discreteDistCfgVec
            for continuousDistTypeVecIdx=1:length(continuousDistTypeVec)
                copulaType = copulaTypeVec{copulaTypeVecIdx};
                if(strcmpi(copulaType, 'Gaussian'))
                    alpha = rhoVec(alphaVecIdx);
                else
                    alpha = alphaVec(alphaVecIdx);
                end

                continuousDistType = continuousDistTypeVec{continuousDistTypeVecIdx};

                % apply the icdf transform to U(:,2)
                % make both X1 and X2 multimodal distributions
                if(strcmp(continuousDistType, 'Multimodal'))
                    xContinuous = [normrnd(-2,0.3,M/2,1); normrnd(2,0.8,M/2,1)];
                elseif(strcmp(continuousDistType, 'Uniform'))
                    xContinuous = unifrnd(-2,2,M,1);
                elseif(strcmp(continuousDistType, 'Gaussian'))
                    xContinuous = normrnd(2,0.5,M,1);
                elseif(strcmp(continuousDistType, 'ThickTailed'))
                    xContinuous = trnd(1, M, 1);
                else
                    error('Unknown X2 Dist Type!');
                end
                xContinuous = xContinuous(randperm(M),:);     % permute for evenness of samples
                [fContinous,xiContinuous] = emppdf(xContinuous,0);
                FContinuous = empcdf(xContinuous,0);
                continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous);

                % manual toggling between uniform and skewed
                if(discreteDistCfg==1)
                    prob = [0.5 0.5];
                    prob = [0.8 0.2];
                elseif(discreteDistCfg==2)
                    prob = .25*ones(1,4);
                    prob = [0.6 0.25 0.1 0.05];
                elseif(discreteDistCfg==3)
                    prob = .1*ones(1,10);
                    prob = [0.5 0.15 0.15 0.05 0.05 0.05 0.025 0.025];
                end
                a_dist = makedist('Multinomial','Probabilities',prob);

                sRhoRescaled = 0;
                sRhoActual = copulastat(copulaType, alpha, 'type', 'Spearman');
                for mcSimNum=1:numMCSims
                    U = copularnd(copulaType, alpha, M);
                    X_hybrid(:,1) = a_dist.icdf(U(:,1));
                    for ii=1:M
                        X_hybrid(ii,2) = continuousDistInfo.icdf(U(ii,2));
                    end
                    X_hybrid_continued = X_hybrid;
                    X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));
                    
                    U_hybrid = pseudoobs(X_hybrid);
                    tau1 = corr(U_hybrid(:,1),U_hybrid(:,2),'type','spearman');
                    sRhoRescaled = sRhoRescaled + tau1/(12*sqrt(prod(var(U_hybrid))));
                end
                sRhoRescaled = sRhoRescaled/numMCSims;

                sRhoRescaledVec(copulaTypeVecIdx, discreteDistCfg, continuousDistTypeVecIdx, alphaVecIdx) = sRhoRescaled;
                sRhoActualVec(copulaTypeVecIdx, alphaVecIdx) = sRhoActual;
            end
        end
    end
end

titleAppendStr = 'Left-Skew';
figure(1);
subplot(2,2,1);
plot(alphaVec, squeeze(sRhoRescaledVec(1,1,1,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,1,1,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,1,1,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,1,1,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Multimodal Continuous');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,2);
plot(alphaVec, squeeze(sRhoRescaledVec(1,1,2,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,1,2,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,1,2,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,1,2,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Uniform Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,3);
plot(alphaVec, squeeze(sRhoRescaledVec(1,1,3,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,1,3,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,1,3,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,1,3,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Gaussian Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,4);
plot(alphaVec, squeeze(sRhoRescaledVec(1,1,4,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,1,4,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,1,4,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,1,4,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Thick-Tailed Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');
mtit(sprintf('CFG=1 -- %s', titleAppendStr));

figure(2);
subplot(2,2,1);
plot(alphaVec, squeeze(sRhoRescaledVec(1,2,1,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,2,1,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,2,1,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,2,1,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+'); 
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Multimodal Continuous');
xlabel('\alpha');
ylabel('\Delta \tau');

subplot(2,2,2);
plot(alphaVec, squeeze(sRhoRescaledVec(1,2,2,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,2,2,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,2,2,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,2,2,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Uniform Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,3);
plot(alphaVec, squeeze(sRhoRescaledVec(1,2,3,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,2,3,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,2,3,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,2,3,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Gaussian Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,4);
plot(alphaVec, squeeze(sRhoRescaledVec(1,2,4,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,2,4,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,2,4,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,2,4,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+'); 
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Thick-Tailed Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');
mtit(sprintf('CFG=2 -- %s', titleAppendStr));

figure(3);
subplot(2,2,1);
plot(alphaVec, squeeze(sRhoRescaledVec(1,3,1,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,3,1,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,3,1,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,3,1,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+'); 
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Multimodal Continuous');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,2);
plot(alphaVec, squeeze(sRhoRescaledVec(1,3,2,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,3,2,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,3,2,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,3,2,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Uniform Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,3);
plot(alphaVec, squeeze(sRhoRescaledVec(1,3,3,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,3,3,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,3,3,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,3,3,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+');
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Gaussian Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');

subplot(2,2,4);
plot(alphaVec, squeeze(sRhoRescaledVec(1,3,4,:)), 'b', ...
     alphaVec, sRhoActualVec(1,:), 'b+', ...
     alphaVec, squeeze(sRhoRescaledVec(2,3,4,:)), 'r', ...
     alphaVec, sRhoActualVec(2,:), 'r+', ...
     alphaVec, squeeze(sRhoRescaledVec(3,3,4,:)), 'k', ...
     alphaVec, sRhoActualVec(3,:), 'k+', ...
     alphaVec, squeeze(sRhoRescaledVec(4,3,4,:)), 'g', ...
     alphaVec, sRhoActualVec(4,:), 'g+'); 
grid on; 
legend('Frank-Rescaled', 'Frank-Actual', ...
       'Gumbel-Rescaled', 'Gumbel-Actual', ...
       'Clayton-Rescaled', 'Clayton-Actual', ...
       'Gaussian-Rescaled', 'Gaussian-Actual'); 
title('Thick-Tailed Continuous Distribution');
xlabel('\alpha');
ylabel('\tau_R');
mtit(sprintf('CFG=3 -- %s', titleAppendStr));
%% D = 2
clear;
clc;

K = 100; h = 0.05;
copulaTypeVec = {'Frank', 'Gumbel', 'Clayton', 'Gaussian'};
alphaVec = [1 8 20];
continuousDistTypeVec = {'Multimodal', 'Uniform', 'Gaussian', 'ThickTailed'};
mVec = 1000;
RhoVecs_2D = cell(1,length(alphaVec)); 
RhoVecs_2D{1} = [1 -0.9; -0.9 1]; RhoVecs_2D{2} = [1 -0.65; -0.65 1];
RhoVecs_2D{3} = [1 0.35; 0.35 1]; RhoVecs_2D{4} = [1 0.1; 0.1 1];

discreteDistTypeDescriptionVec = cell(1,9);
discreteDistTypeVec = cell(1,9);
prob = [0.25 0.25 0.25 0.25];   % uniform
discreteDistTypeVec{1} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{1} = '4-Uniform';
prob = [0.6 0.2 0.15 0.05];     % left skew
discreteDistTypeVec{2} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{2} = '4-SkewLeft';
prob = fliplr([0.6 0.2 0.15 0.05]);     % right skew
discreteDistTypeVec{3} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{3} = '4-SkewRight';
prob = [0.5 0.5];
discreteDistTypeVec{4} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{4} = '2-Uniform';
prob = [0.8 0.2];               % left skew
discreteDistTypeVec{5} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{5} = '2-SkewLeft';
prob = fliplr([0.8 0.2]);               % right skew
discreteDistTypeVec{6} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{6} = '2-SkewRight';
prob = [0.125 0.125 0.125 0.125 0.125 0.125 0.125 0.125];  % uniform
discreteDistTypeVec{7} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{7} = '8-Uniform';
prob = [0.5 0.15 0.15 0.05 0.05 0.05 0.025 0.025];  % skew left
discreteDistTypeVec{8} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{8} = '8-SkewLeft';
prob = fliplr([0.5 0.15 0.15 0.05 0.05 0.05 0.025 0.025]);  % skew right
discreteDistTypeVec{9} = makedist('Multinomial','Probabilities',prob); discreteDistTypeDescriptionVec{9} = '8-SkewRight';

fig_num = 1;
numTotalLoops = length(copulaTypeVec)*length(continuousDistTypeVec)*length(alphaVec)*length(discreteDistTypeVec)*length(mVec);
fid = fopen('/Users/kiran/ownCloud/PhD/sim_results/idea/progress.txt', 'a');
for copulaTypeVecIdx=1:length(copulaTypeVec)
    for continuousDistTypeVecIdx=1:length(continuousDistTypeVec)
        for alphaVecIdx=1:length(alphaVec)
            for discreteDistTypeVecIdx=1:length(discreteDistTypeVec)
                for mVecIdx=1:length(mVec)
                    copulaType = copulaTypeVec{copulaTypeVecIdx};
                    continuousDistType = continuousDistTypeVec{continuousDistTypeVecIdx};
                    discreteDist = discreteDistTypeVec{discreteDistTypeVecIdx};
                    discreteDistType = discreteDistTypeDescriptionVec{discreteDistTypeVecIdx};
                    if(strcmp(copulaType, 'Gaussian'))
                        alpha = RhoVecs_2D{alphaVecIdx};        % alpha is Rho here
                        alphaPrint = alpha(1,2);
                    else
                        alpha = alphaVec(alphaVecIdx);
                        alphaPrint = alpha;
                    end
                    M = mVec(mVecIdx);
                    
                    fprintf(fid, 'Progress=%0.02f %%\n', (fig_num/numTotalLoops)*100);
                    fprintf('Progress=%0.02f %%\n', (fig_num/numTotalLoops)*100);
                    
                    % Generate U
                    U = copularnd(copulaType, alpha, M);

                    % apply the icdf transform to U(:,1)
                    X_hybrid(:,1) = discreteDist.icdf(U(:,1));

                    % apply the icdf transform to U(:,2)
                    % make both X1 and X2 multimodal distributions
                    if(strcmp(continuousDistType, 'Multimodal'))
                        xContinuous = [normrnd(-2,0.3,M/2,1); normrnd(2,0.8,M/2,1)];
                    elseif(strcmp(continuousDistType, 'Uniform'))
                        xContinuous = unifrnd(-2,2,M,1);
                    elseif(strcmp(continuousDistType, 'UnimodalSkewed'))
                        xContinuous = betarnd(2,5,M,1);
                    elseif(strcmp(continuousDistType, 'Gaussian'))
                        xContinuous = normrnd(2,0.5,M,1);
                    elseif(strcmp(continuousDistType, 'ThickTailed'))
                        xContinuous = trnd(1, M, 1);
                    else
                        error('Unknown X2 Dist Type!');
                    end
                    xContinuous = xContinuous(randperm(M),:);     % permute for evenness of samples
                    [fContinous,xiContinuous] = emppdf(xContinuous,0);
                    FContinuous = empcdf(xContinuous,0);
                    continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous);
                    for ii=1:M
                        X_hybrid(ii,2) = continuousDistInfo.icdf(U(ii,2));
                    end

                    X_hybrid_continued = X_hybrid;
                    X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));
                    U_hybrid_continued = pseudoobs(X_hybrid_continued);
                    U_hybrid = pseudoobs(X_hybrid);

                    % estimate the copula for U_hybrid_continued
                    c_est = empcopulapdf(U_hybrid_continued, h, K, 'betak');
                    [U1,U2] = ndgrid(linspace(0,1,K));
                    
%                     fig1=figure(1);
                    set(gcf, 'Position', get(0,'Screensize'));
                    subplot(2,6,[1 2]); scatter(U(:,1),U(:,2)); 
                    title(sprintf('%s alpha=%0.02f', copulaType, alphaPrint)); grid on; xlabel('u_1'); ylabel('u_2');
                    subplot(2,6,[3 4]); scatter(U_hybrid(:,1),U_hybrid(:,2)); 
                    title(sprintf('%s | %s Pseudo-Obs', continuousDistType, discreteDistType)); grid on; xlabel('u_1'); ylabel('u_2');
                    subplot(2,6,[5 6]); scatter(U_hybrid_continued(:,1),U_hybrid_continued(:,2));
                    title('Continued Pseudo-Observations'); grid on; xlabel('u_1'); ylabel('u_2');
                    
                    subplot(2,6,[7 8 9]); surf(U1,U2,reshape( copulapdf(copulaType, [U1(:), U2(:)], alpha), K, K)); title('Actual Copula'); xlabel('u_1'); ylabel('u_2'); grid on;
                    subplot(2,6,[10 11 12]); surf(U1,U2,c_est); title('Estimated Copula'); xlabel('u_1'); ylabel('u_2'); grid on;                    
                    figureFilename = sprintf('/Users/kiran/ownCloud/PhD/sim_results/idea/fig_%d', fig_num);
                    fig = gcf;
                    fig.PaperPositionMode = 'auto';
                    print(figureFilename,'-dpng','-r0')
                    fig_num = fig_num + 1;
%                     pause(1);
                    pause;
%                     clf(fig1);
                end
            end
        end
    end
end
fclose(fid);