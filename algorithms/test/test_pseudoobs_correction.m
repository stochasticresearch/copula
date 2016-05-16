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
%**************************************************************************

% generate 2-D hybrid data w/ different dependency, and show how our
% correction improves the copula density estimate
alphaVec = [1 10];
copulaTypeVec = {'Clayton', 'Gaussian'};
RhoVecs_2D = [-0.1 -0.9];  
M = 1000;
% continuousDistTypeVec = {'Multimodal', 'Uniform', 'Gaussian'};
continuousDistTypeVec = {'Multimodal'};
discreteDistCfgVec = [1 2];
skewVec = [1 2];

correctionFactor = .75;
h = 0.05;
K = 25;
numMCSims = 1;

for copulaTypeVecIdx=1:length(copulaTypeVec)
    for discreteDistCfg=discreteDistCfgVec
        for skewVecCfg=skewVec
            for continuousDistTypeVecIdx=1:length(continuousDistTypeVec)
                for alphaVecIdx=1:length(alphaVec)
                    copulaType = copulaTypeVec{copulaTypeVecIdx};
                    if(strcmpi(copulaType, 'Gaussian'))
                        alpha = RhoVecs_2D(alphaVecIdx);
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
                        if(skewVecCfg==1)
                            prob = [0.5 0.5];
                        else
                            prob = [0.8 0.2];
                        end
                    elseif(discreteDistCfg==2)
                        if(skewVecCfg==1)
                            prob = .1*ones(1,10);
                        else
                            prob = [0.5 0.15 0.15 0.05 0.05 0.05 0.025 0.025];
                        end
                    end
                    a_dist = makedist('Multinomial','Probabilities',prob);

                    [U1,U2] = ndgrid(linspace(0,1,K));
                    c_actual = copulapdf(copulaType, [U1(:) U2(:)], alpha);
                    for mcSimNum=1:numMCSims
                        U = copularnd(copulaType, alpha, M);
                        X_hybrid(:,1) = a_dist.icdf(U(:,1));
                        for ii=1:M
                            X_hybrid(ii,2) = continuousDistInfo.icdf(U(ii,2));
                        end
                        X_hybrid_continued = X_hybrid;
                        X_hybrid_continued(:,1) = continueRv(X_hybrid(:,1));

                        U_hybrid = pobs(X_hybrid);
                        U_hybrid_continued = pobs(X_hybrid_continued);
                        sRho = ndcorr(U_hybrid, 'spearman');
                        U_hybrid_continued_corrected = pobs(X_hybrid_continued, 'correction', sRho, correctionFactor);

                        % estimate the copula w/ U_hybrid and
                        % U_hybrid_corrected
                        c_est = empcopulapdf(U_hybrid_continued, h, K, 'betak');
                        c_est_corrected = empcopulapdf(U_hybrid_continued_corrected, h, K, 'betak');

                        % compare the MSE between these 2 estimators
                        c_est_mse = sum( (c_est(:)-c_actual).^2 );
                        c_est_corrected_mse = sum( (c_est_corrected(:)-c_actual).^2 );

                        % TODO: store these
                    end
                    
                    if(sRho>=0)
                        jamiexx = linspace(0,1,K);
                        jamieyy = sRho*jamiexx;
                    else
                        jamiexx = linspace(0,1,K);
                        jamieyy = sRho*jamiexx + 1;
                    end
                    
                    
                    f = figure(1);
                    subplot(2,3,1);
                    scatter(U_hybrid(:,1), U_hybrid(:,2)); grid on;  hold on;
                    scatter(jamiexx,jamieyy, 'g+');
                    title(sprintf('U Hybrid - %s -- %0.02f', continuousDistType, sRho));
                    
                    subplot(2,3,2);
                    scatter(U_hybrid_continued(:,1), U_hybrid_continued(:,2)); grid on; 
                    if(skewVecCfg==1)
                        title('U Hybrid Continued -- No Skew');
                    else
                        title('U Hybrid Continued -- Left Skew');
                    end
                    hold on;
                    scatter(U(:,1),U(:,2), 'r');
                    
                    subplot(2,3,3);
                    scatter(U_hybrid_continued_corrected(:,1), U_hybrid_continued_corrected(:,2)); grid on;
                    title('U Hybrid Corrected');
                    
                    subplot(2,3,4);
                    surf( U1, U2, reshape(c_actual, K, K) ); grid on; 
                    title(sprintf('%s | theta=%d', copulaType, alpha));
                    xlabel('U_1'); ylabel('U_2');
                    
                    subplot(2,3,5);
                    surf( U1, U2, c_est ); grid on; 
                    title(sprintf('C Est - MSE=%0.03f', c_est_mse));
                    xlabel('U_1'); ylabel('U_2');
                    
                    subplot(2,3,6);
                    surf( U1, U2, c_est_corrected ); grid on; 
                    title(sprintf('C Est Corrected - MSE=%0.03f', c_est_corrected_mse));
                    xlabel('U_1'); ylabel('U_2');

                    pause;
                    clf(f);
                end
            end
        end
    end
end