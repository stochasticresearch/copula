%% Generates the Type I error curves for testing conditional independence
clear;
clc;

rng(315);
dbstop if error;

nsim = 500;

num_noise = 30;
noise = 3;

M = 500;

numDepTests = 6;        % the number of different dependency tests we will conduct
                        % TODO: add copula dependencies as well

rscdmTypeINull = zeros(1,nsim);
rscdmTypeIAlt  = zeros(1,nsim);
rscdmTypeIPower = zeros(numDepTests, num_noise);

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
minVal = 0;
maxVal = 1;

testAlpha = 0.05;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the Type I Error simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 1;
num_noise_test_max = 20;
for l=num_noise_test_min:num_noise_test_max
    for typ=1:numDepTests
        dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        parfor ii=1:nsim
            dispstat(sprintf('Simulating -- %0.02f', ii/nsim*100),'timestamp');
            % Generate data from     Y-->X<--Z
            y = rand(M,1)*(maxVal-minVal)+minVal;
            z = rand(M,1)*(maxVal-minVal)+minVal;
            switch(typ)
                case 1
                    % linear
                    x = y + z + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    x = 4*(y-.5).^2 + 4*(z-.5).^2 + ...
                        noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    x = 128*(y-1/3).^3-48*(y-1/3).^3-12*(y-1/3)+ ...
                        128*(z-1/3).^3-48*(z-1/3).^3-12*(z-1/3)+ ...
                        10* noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    x = sin(2*pi*y) + sin(2*pi*z) + ...
                        2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % fourth root
                    x = y.^(1/4) + z.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % circle
                    binomailVals = (2*binornd(1,0.5,M,1)-1);
                    x = binomailVals .* (sqrt(1 - (2*y - 1).^2)) + ...
                        binomailVals .* (sqrt(1 - (2*z - 1).^2)) + ...
                        noise/4*l/num_noise*randn(M,1);
                % TODO: add copula based tests
                
                % TODO: decide whether to add post-nonlinear based tests
                % also?
                    
                otherwise
                    error('unknown dep type!');
            end
            
            % calculate the metric
            metricConditioned = rscdm(y, z, x);
            metricUnconditioned = rsdm(y, z);
            
            rscdmTypeIAlt(ii) = metricConditioned;
            rscdmTypeINull(ii) = metricUnconditioned;
        end
        
        % compute the cutoff
        rscdm_cut = quantile(rscdmTypeINull, 1-testAlpha);
        % here we look for >, b/c we want to see that once it is
        % conditioned upon the variable, we get a lower value, meaning we
        % would declare conditional independence
        rscdmTypeIPower(typ, l) = sum(rscdmTypeIAlt > rscdm_cut)/nsim;
    end
end

% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rscdmPower_TypeI.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rscdmPower_TypeI.mat');
end

% TODO: do the inlet for the plots instead of a title

noiseVec = (num_noise_test_min:num_noise_test_max)/10;
h1 = subplot(3,2,1);
hh1 = plot(noiseVec, rscdmTypeIPower(1,num_noise_test_min:num_noise_test_max), 'o-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Type I Power'); grid on;
title('Linear');
h1.FontSize = 20; 

h2 = subplot(3,2,2);
hh2 = plot(noiseVec, rscdmTypeIPower(2,num_noise_test_min:num_noise_test_max), 'o-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Type I Power'); grid on;
title('Parabolic')
h2.FontSize = 20; 

h3 = subplot(3,2,3); 
hh3 = plot(noiseVec, rscdmTypeIPower(3,num_noise_test_min:num_noise_test_max), 'o-.');  
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Type I Power'); grid on;
title('Cubic');
h3.FontSize = 20; 

h4 = subplot(3,2,4); 
hh4 = plot(noiseVec, rscdmTypeIPower(4,num_noise_test_min:num_noise_test_max), 'o-.');  
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Type I Power'); grid on;
title('Sine')
h4.FontSize = 20; 

h5 = subplot(3,2,5); 
hh5 = plot(noiseVec, rscdmTypeIPower(5,num_noise_test_min:num_noise_test_max), 'o-.');
xlabel('Noise Level'); ylabel('Type I Power'); grid on;
title('Fourth-Root');
h5.FontSize = 20; 

h6 = subplot(3,2,6); 
hh6 = plot(noiseVec, rscdmTypeIPower(6,num_noise_test_min:num_noise_test_max), 'o-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Type I Power'); grid on;
title('Circle');
h6.FontSize = 20; 

%% Generates the Type II error curves for testing conditional independence
clear;
clc;

rng(315);
dbstop if error;

nsim = 500;

num_noise = 30;
noise = 3;

M = 500;

numDepTests = 6;        % the number of different dependency tests we will conduct
                        % TODO: add copula dependencies as well

rscdmTypeIINull = zeros(1,nsim);
rscdmTypeIIAlt  = zeros(1,nsim);
rscdmTypeIIPower = zeros(numDepTests, num_noise);

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
xMin = 0;
xMax = 1;

testAlpha = 0.05;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the Type II Error simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 1;
num_noise_test_max = 20;
for l=num_noise_test_min:num_noise_test_max
    for typ=1:numDepTests
        dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        parfor ii=1:nsim
            dispstat(sprintf('Simulating -- %0.02f', ii/nsim*100),'timestamp');
            % Generate data from     Y<--X-->Z
            x = rand(M,1)*(xMax-xMin)+xMin;
            switch(typ)
                case 1
                    % linear
                    y = x + noise*(l/num_noise)*randn(M,1); 
                    z = x + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    y = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                    z = 4*(x-.5).^2 + noise*(l/num_noise)*randn(M,1);
                case 3
                    % cubic
                    y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                    z = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3)+10* noise*(l/num_noise)*randn(M,1);
                case 4
                    % low-freq sin
                    y = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                    z = sin(4*pi*x) + 2*noise*(l/num_noise)*randn(M,1);
                case 5
                    % fourth root
                    y = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                    z = x.^(1/4) + noise*(l/num_noise)*randn(M,1);
                case 6
                    % circle
                    y = (2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                    z = (2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2)) + noise/4*l/num_noise*randn(M,1);
                
                % TODO: add copula based tests
                
                % TODO: decide whether to add post-nonlinear based tests
                % also?
                    
                otherwise
                    error('unknown dep type!');
            end
            % calculate the metric
            metricConditioned = rscdm(y, z, x);
            metricUnconditioned = rsdm(y, z);
            
            rscdmTypeIIAlt(ii) = metricConditioned;
            rscdmTypeIINull(ii) = metricUnconditioned;
        end
        
        % compute the cutoff - we do alpha, not 1-alpha here b/c we are
        % looking to see when the test will fail on the lower-end
        rscdm_cut = quantile(rscdmTypeIINull, testAlpha);
        % here we look for <, b/c we want to see that once it is
        % conditioned upon the variable, we get a lower value, meaning we
        % would declare conditional independence
        rscdmTypeIIPower(typ, l) = sum(rscdmTypeIIAlt < rscdm_cut)/nsim;
    end
end


% save the data
if(ispc)
    save('C:\\Users\\Kiran\\ownCloud\\PhD\\sim_results\\independence\\rscdmPower_TypeII.mat');
else
    save('/home/kiran/ownCloud/PhD/sim_results/independence/rscdmPower_TypeII.mat');
end

% TODO: do the inlet for the plots instead of a title

noiseVec = (num_noise_test_min:num_noise_test_max)/10;
h1 = subplot(3,2,1);
hh1 = plot(noiseVec, rscdmTypeIIPower(1,num_noise_test_min:num_noise_test_max), 'o-.'); 
xlabel('Noise Level'); ylabel('Type II Power'); grid on;
title('Linear');
h1.FontSize = 20; 

h2 = subplot(3,2,2);
hh2 = plot(noiseVec, rscdmTypeIIPower(2,num_noise_test_min:num_noise_test_max), 'o-.'); 
xlabel('Noise Level'); ylabel('Type II Power'); grid on;
title('Parabolic')
h2.FontSize = 20; 

h3 = subplot(3,2,3); 
hh3 = plot(noiseVec, rscdmTypeIIPower(3,num_noise_test_min:num_noise_test_max), 'o-.');  
xlabel('Noise Level'); ylabel('Type II Power'); grid on;
title('Cubic');
h3.FontSize = 20; 

h4 = subplot(3,2,4); 
hh4 = plot(noiseVec, rscdmTypeIIPower(4,num_noise_test_min:num_noise_test_max), 'o-.');  
xlabel('Noise Level'); ylabel('Type II Power'); grid on;
title('Sine')
h4.FontSize = 20; 

h5 = subplot(3,2,5); 
hh5 = plot(noiseVec, rscdmTypeIIPower(5,num_noise_test_min:num_noise_test_max), 'o-.');
xlabel('Noise Level'); ylabel('Type II Power'); grid on;
title('Fourth-Root');
h5.FontSize = 20; 

h6 = subplot(3,2,6); 
hh6 = plot(noiseVec, rscdmTypeIIPower(6,num_noise_test_min:num_noise_test_max), 'o-.'); 
axis([min(noiseVec) max(noiseVec) 0 1]);
xlabel('Noise Level'); ylabel('Type II Power'); grid on;
title('Circle');
h6.FontSize = 20; 

%% Do a conditionally independent test

clear;
clc;

rng(123);

M = 500;
noise = 0.01;
alpha = 0.05;

% Generate data from     Y<--X-->Z
x = rand(M,1);
% y = 4*(x-0.5).^2 + noise*randn(M,1);
% z = 4*(x-0.5).^2 + noise*randn(M,1);
% y = x + noise*randn(M,1);
% z = x + noise*randn(M,1);
y = sin(2*pi*x) + noise*randn(M,1);
z = 4*(x-0.5).^2 + noise*randn(M,1);

rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% RSCDM conditions X on Y and Z, and sees how related Y and Z
% are to each other ... in this graphical model, they should be UNRELATED
% after the effect of X is removed... i.e. close to independent.  To see
% why, look at the graphical model, Y indep Z | X according to
% d-separation.  So if we condition upon X (i.e. remove teh effect of X on
% Y and Z separately), then we should get independence.
[rscdmVal, RxStacked, RyStacked, RxPtsStacked, RyPtsStacked] = rscdm(y,z,x);
rscdmVal_Z = sqrt(M-4)*abs(rscdmVal);
testVal = norminv(1-alpha/2);

subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm1), 'FontSize', 20);

subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm3), 'FontSize', 20);

subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm2), 'FontSize', 20);

subplot(3,12,13:17);
scatter(RxPtsStacked, RxStacked); grid on;
xlabel('x', 'FontSize', 20); ylabel('r_y', 'FontSize', 20);

subplot(3,12,20:24);
scatter(RyPtsStacked, RyStacked); grid on;
xlabel('x', 'FontSize', 20); ylabel('r_z', 'FontSize', 20);

subplot(3,12,25:29);
scatter(RxStacked,RyStacked); grid on; 
xlabel('r_y', 'FontSize', 20); ylabel('r_z', 'FontSize', 20);  
title(sprintf('rscdm=%0.2f', rscdmVal), 'FontSize', 20);

subplot(3,12,32:36);
scatter(pobs(RxStacked),pobs(RyStacked), 'r'); grid on; 
xlabel('F_{r_y}', 'FontSize', 20); ylabel('F_{r_z}', 'FontSize', 20);


%% Do a conditionally dependent test

% Tests Conditional Independence w/ RSDM and the residuals processing
% algorithm.

clear;
clc;

rng(123);

M = 500;
noise = 0.5;
alpha = 0.05;

% Generate data from     Y-->X<--Z
y = rand(M,1);
z = rand(M,1);
x = y.^2+z.^3 + noise;      % TODO: fix the bug :D
% x = y + z + noise;


rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% In this graphical model, Y and Z are independent of each other, but when
% conditioned upon X, they become dependent.  Refer to the rules of
% d-separation to see why, this is a V-Structure!
[rscdmVal, RxStacked, RyStacked, RxPtsStacked, RyPtsStacked] = rscdm(y,z,x);
rscdmVal_Z = sqrt(M-4)*abs(rscdmVal);
testVal = norminv(1-alpha/2);

subplot(3,12,1:3);
scatter(x,y); grid on; xlabel('x', 'FontSize', 20); ylabel('y', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm1), 'FontSize', 20);

subplot(3,12,5:8);
scatter(y,z); grid on; xlabel('y', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm3), 'FontSize', 20);

subplot(3,12,10:12);
scatter(x,z); grid on; xlabel('x', 'FontSize', 20); ylabel('z', 'FontSize', 20); 
title(sprintf('RSDM=%0.2f', rsdm2), 'FontSize', 20);

subplot(3,12,13:17);
scatter(RxPtsStacked, RxStacked); grid on;
xlabel('x', 'FontSize', 20); title('r_y', 'FontSize', 20);

subplot(3,12,20:24);
scatter(RyPtsStacked, RyStacked); grid on;
xlabel('x', 'FontSize', 20); title('r_z', 'FontSize', 20);

subplot(3,12,25:29);
scatter(RxStacked,RyStacked); grid on; xlabel('r_y', 'FontSize', 20); ylabel('r_z', 'FontSize', 20);  
title(sprintf('rscdm=%0.2f', rscdmVal), 'FontSize', 20);

subplot(3,12,32:36);
scatter(pobs(RxStacked),pobs(RyStacked), 'r'); grid on; xlabel('F_{r_y}'); ylabel('F_{r_z}');
