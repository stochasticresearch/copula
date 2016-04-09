% A script which shows the effect of continuing a RV on
% pseudo-observations, which are then used to estimate the copula density

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
fid = fopen('/home/kiran/ownCloud/PhD/sim_results/idea/progress.txt', 'a');
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
                    X_hybrid(:,1) = discreteDist.icdf(U(:,1));
                    for ii=1:M
                        X_hybrid(ii,2) = continuousDistInfo.invDistribution(U(ii,2));
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
                    figureFilename = sprintf('/home/kiran/ownCloud/PhD/sim_results/idea/fig_%d', fig_num);
                    fig = gcf;
                    fig.PaperPositionMode = 'auto';
                    print(figureFilename,'-dpng','-r0')
                    fig_num = fig_num + 1;
                    pause(1);
%                     clf(fig1);
                end
            end
        end
    end
end
fclose(fid);