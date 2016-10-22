%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

% Find the optimal hyper-parameters for RSDM, which equates to 

%% Script which determines the optimal parameters for the RSDM metric
clear;
clc;

depsToOptimize = {'independence', 'linear', 'quadratic', 'cubic', ...
                  'sine-1', 'sine-2', 'fourth-root', 'circular', 'step'};

nsim = 500;
% for independence, we find the optimal parameters which minimize the
% metric, for everything else we find the parameters which maximize the
% metric (this should have the effect of maximizing the statistical power
% of the RSDM metric);

% CFG-1 OPTIMIZATION TEST
cfg = 1;
minscanincrVec = [0.025 0.05 0.1];
diffthreshVec = 20:20:200;
alphaVec = 0.01:0.01:0.1;

results = zeros(length(depsToOptimize), length(minscanincrVec), length(diffthreshVec), length(alphaVec) );

M = 500;
dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...\n'),'keepthis','timestamp');

numLoopsPerDep = length(minscanincrVec)*length(diffthreshVec)*length(alphaVec)*nsim;
jj = 1;
for depToOptimize=depsToOptimize
    dispstat(sprintf('Simulating Dependency Type %s -- %0.02f', char(depToOptimize), jj/length(depsToOptimize)*100 ),'keepthis','timestamp');
    kk = 1;
    progressIdx = 1;
    for minscanincr=minscanincrVec
        ll = 1;
        for diffthresh=diffthreshVec
            mm = 1;
            for alpha=alphaVec
                nn = 1;
                metricVec = zeros(1,nsim);
                parfor ii=1:nsim
%                         dispstat(sprintf('Monte-Carlo -- %0.02f', progressIdx/numLoopsPerDep*100),'timestamp');

                    % generate X & Y
                    x = rand(M,1);
                    noise = 1.5*randn(M,1);
                    switch(char(depToOptimize))
                        case 'independence'
                            y = rand(M,1) + noise;
                        case 'linear'
                            y = x + noise; 
                        case 'quadratic'
                            y = 4*(x-.5).^2 + noise;
                        case 'cubic'
                            y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3) + noise;
                        case 'sine-1'
                            y = sin(4*pi*x) + 2*noise;      % noise scaling factor same as above .. not 100% sure why though
                        case 'sine-2'
                            y = sin(16*pi*x) + noise;
                        case 'fourth-root'
                            y = x.^(1/4) + noise;
                        case 'circular'
                            y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + 3/4*15/30*randn(M,1);  % noise=3, num_noise=30, l=15
                        case 'step'
                            y = (x > 0.5) + 3*5*15/30*randn(M,1); % noise=3, num_noise=30, l=15
                        otherwise
                            error('unrecognized dep!');
                    end
                    % compute a running average of the copula metric
                    metricVec(ii) = rsdm_2(x,y, minscanincr, diffthresh, alpha);
                    progressIdx = progressIdx + 1;
                end
                switch(char(depToOptimize))
                    case 'independence'
                        metric = quantile(metricVec, 0.95);
                    otherwise
                        metric = quantile(metricVec, 0.05);
                end
                results(jj,kk,ll,mm) = metric;
                mm = mm + 1;
            end
            ll = ll + 1;
        end
        kk = kk + 1;
    end
    jj = jj + 1;
end

% find the configuration which maximizes the *overall* power for all
% different types of dependencies

% define some convenience variables
linearResults = squeeze(results(2,:,:,:)-results(1,:,:,:));
quadraticResults = squeeze(results(3,:,:,:)-results(1,:,:,:));
cubicResults = squeeze(results(4,:,:,:)-results(1,:,:,:));
sine1Results = squeeze(results(5,:,:,:)-results(1,:,:,:));
sine2Results = squeeze(results(6,:,:,:)-results(1,:,:,:));
fourthRootResults = squeeze(results(7,:,:,:)-results(1,:,:,:));
circularResults = squeeze(results(8,:,:,:)-results(1,:,:,:));
stepResults = squeeze(results(9,:,:,:)-results(1,:,:,:));

optimfound = 0;
thresh = 0;
threshdecr = 0.0125;
minthresh = -0.02;
while(~optimfound || thresh>=minthresh)
    lIdx = find(linearResults>thresh);
    qIdx = find(quadraticResults>thresh);
    cIdx = find(cubicResults>thresh);
    s1Idx = find(sine1Results>thresh);
    s2Idx = find(sine2Results>thresh);
    frIdx = find(fourthRootResults>thresh);
    ciIdx = find(circularResults>thresh);
    stIdx = find(stepResults>thresh);
    
    bestIdxs = intersect(intersect(intersect(intersect(intersect(intersect(intersect(lIdx,qIdx),cIdx),s1Idx),s2Idx),frIdx),ciIdx),stIdx);
    
    if(isempty(bestIdxs))
        thresh = thresh - threshdecr;
    else
        optimfound = 1;
    end
end

if(ispc)
    save(sprintf('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rsdm2_optimal_params_cfg%d.mat',cfg));
elseif(ismac)
    save(sprintf('/Users/kiran/ownCloud/PhD/sim_results/independence/rsdm2_optimal_params_cfg%d.mat',cfg));
else
    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/independence/rsdm2_optimal_params_cfg%d.mat',cfg));
end

for ii=1:length(bestIdxs)
    [I1,I2] = ind2sub(size(squeeze(results(1,:,:,:))),bestIdxs(ii));
    fprintf('Best Idx = [%d %d]\n', I1, I2);
    % print the expected "power" for each of these configurations so we can
    % choose the best one manually
    fprintf('Linear = %0.03f\n', results(2,1,I1,I2)-results(1,1,I1,I2));
    fprintf('Quadratic = %0.03f\n', results(3,1,I1,I2)-results(1,1,I1,I2));
    fprintf('Cubic = %0.03f\n', results(4,1,I1,I2)-results(1,1,I1,I2));
    fprintf('Sine-1 = %0.03f\n', results(5,1,I1,I2)-results(1,1,I1,I2));
    fprintf('Sine-2 = %0.03f\n', results(6,1,I1,I2)-results(1,1,I1,I2));
    fprintf('FourthRoot = %0.03f\n', results(7,1,I1,I2)-results(1,1,I1,I2));
    fprintf('Circular = %0.03f\n', results(8,1,I1,I2)-results(1,1,I1,I2));
    fprintf('Step = %0.03f\n', results(9,1,I1,I2,I3)-results(1,1,I1,I2));
    fprintf('------------------------------------------------------\n\n');
end