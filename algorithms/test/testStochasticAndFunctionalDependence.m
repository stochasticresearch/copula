%**************************************************************************
%* 
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>
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

%% Try to emulate a Clayton Copula w/ function Dependence + Noise
clear;
clc;

M = 500;
minVal = 0;
maxVal = 1;

isdiscrete = 0;
numMCSims = 100;

X = rand(M,1)*(maxVal-minVal)+minVal; X(X<=0) = eps;
[xPdf, xi] = emppdf(X, isdiscrete);
xCdf = empcdf(X, isdiscrete);
xEmpInfo = rvEmpiricalInfo(xi, xPdf, xCdf, isdiscrete);

u = 0.05:0.05:0.95;
v = 0.05:0.05:0.95;
C_generated = zeros(length(u),length(v));

C_clayton = zeros(length(u),length(v)); alpha = 5;
C_frank = zeros(length(u),length(v));
C_gumbel = zeros(length(u),length(v));
C_gaussian = zeros(length(u),length(v)); rho = 0.5;

noiseVec = 0.05:0.05:.5;
distanceVec = zeros(4,length(noiseVec));
for kk = 1:length(noiseVec)
    noisePower = noiseVec(kk);
    fprintf('Computing for Noise=%0.02f\n', noisePower);
    
    computedDistanceClayton = 0;
    computedDistanceFrank = 0;
    computedDistanceGumbel = 0;
    computedDistanceGaussian = 0;
    for ll=1:numMCSims
        Y = log(X) + randn(M,1)*noisePower;
        % compute the empirical info of Y
        [yPdf, yi] = emppdf(Y, isdiscrete);
        yCdf = empcdf(Y, isdiscrete);
        yEmpInfo = rvEmpiricalInfo(yi, yPdf, yCdf, isdiscrete);

        % compute the empirical 3-D CDF
        [F_empirical,xi_empirical] = ksdensity([X Y],'function','cdf');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % changed internal ksdensity/mvksdensity.m function to compute 100x100
        % grid for 2-D, rather than the default 30.  
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %     % plot for sanity check
    %     x1_empirical = xi_empirical(:,1);
    %     x2_empirical = xi_empirical(:,2);
    %     xx = linspace(min(x1_empirical),max(x1_empirical));
    %     yy = linspace(min(x2_empirical),max(x2_empirical));
    %     [xq,yq] = meshgrid(xx,yy);
    %     orig_state = warning;
    %     warning('off','all');
    %     Z = griddata(x1_empirical,x2_empirical,F_empirical,xq,yq);
    %     warning(orig_state);
    %     surf(xq,yq,Z);
    %     xlabel('x'); ylabel('y');
    %     pause;

        % apply the M copula (b/c we have monotonic increasing)
        for ii=1:length(u)
            for jj=1:length(v)
                uu = u(ii);
                vv = v(jj);
                % compute inverse's
                H_inv = xEmpInfo.icdf(uu);
                G_inv = yEmpInfo.icdf(vv);

                % find closest value in F_empirical & store as our "generated'
                % copula (This is 2nd version of Sklar's Theorem)
                idx = knnsearch(xi_empirical,[H_inv G_inv]);
                C_generated(ii,jj) = F_empirical(idx);

                C_clayton(ii,jj) = copulacdf('Clayton', [uu vv], alpha);
                C_frank(ii,jj) = copulacdf('Frank', [uu vv], alpha);
                C_gumbel(ii,jj) = copulacdf('Gumbel', [uu vv], alpha);
                C_gaussian(ii,jj) = copulacdf('Gaussian', [uu vv], rho);
            end
        end

        % compute distance metric between these functions
        computedDistanceClayton = computedDistanceClayton + sum(sum((C_generated-C_clayton).^2));
        computedDistanceFrank = computedDistanceFrank + sum(sum((C_generated-C_frank).^2));
        computedDistanceGumbel = computedDistanceGumbel + sum(sum((C_generated-C_gumbel).^2));
        computedDistanceGaussian = computedDistanceGaussian + sum(sum((C_generated-C_gaussian).^2));
    end
    % store the average error
    distanceVec(1,kk) = computedDistanceClayton/numMCSims;
    distanceVec(2,kk) = computedDistanceFrank/numMCSims;
    distanceVec(3,kk) = computedDistanceGumbel/numMCSims;
    distanceVec(4,kk) = computedDistanceGaussian/numMCSims;
    
%     subplot(1,2,1);
%     surf(C_generated);
%     title(sprintf('Noise = %0.02f', noisePower));
%     subplot(1,2,2);
%     surf(C_clayton);
%     title(sprintf('Distance = %0.02f', computedDistanceClayton));
%     pause;
    
end

plot(noiseVec,distanceVec(1,:),...
     noiseVec,distanceVec(2,:),...
     noiseVec,distanceVec(3,:),...
     noiseVec,distanceVec(4,:));
legend('Clayton (\alpha=5)','Frank (\alpha=5)','Gumbel (\alpha=5)','Gaussian(\rho=0.5)');
grid on;
xlabel('\sigma^2', 'FontSize',20)
ylabel('Error','FontSize',20);

%% A Script which generates pseudo-samples using functional dependence, and 
% tries to fit to different copulas, and computes the fitting error.
clear;
clc;

M = 500;

% Gaussian Copula = linear dependence + AWGN
x = rand(M,1); x(x<=0) = eps;
y = log(x) + randn(M,1)*.35;
uv = pobs([x y]);
[paramhat_frank,paramci_frank] = copulafit('frank',uv)
[paramhat_gumbel,paramci_gumbel] = copulafit('gumbel',uv)
[paramhat_clayton,paramci_clayton] = copulafit('clayton',uv)