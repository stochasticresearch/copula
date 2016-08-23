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

% TODO: insert optimal parameters here ... these are not optimal
rsdm_rectWidths = [1 0.75 0.5 0.25 0.125];
rsdm_overlapIncr = 0.055;
rsdm_alpha = 0.03;
rsdm_MA_size = 2;
rsdm_peakFinderSel = 8;
rsdm_boostFactor = 1;

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
minVal = 0;
maxVal = 1;

testAlpha = 0.05;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the Type II Error simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 5;
num_noise_test_max = 20;
for l=num_noise_test_min:num_noise_test_max
    for typ=1:numDepTests
        dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        for ii=1:nsim
            dispstat(sprintf('Simulating Null -- %0.02f', ii/nsim*100),'timestamp');
            % Generate data from     Y-->X<--Z
            y = rand(M,1)*(maxVal-minVal)+minVal;
            z = rand(M,1)*(maxVal-minVal)+minVal;
            switch(typ)
                case 1
                    % linear
                    x = y + z + noise*(l/num_noise)*randn(M,1); 
                case 2
                    % parabolic
                    x = 4*(y-.5).^2 + 4*(x-.5).^2 + ...
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
                    x = binomialVals .* (sqrt(1 - (2*y - 1).^2)) + ...
                        binomialVals .* (sqrt(1 - (2*z - 1).^2)) + ...
                        noise/4*l/num_noise*randn(M,1);
                % TODO: add copula based tests
                
                % TODO: decide whether to add post-nonlinear based tests
                % also?
                    
                otherwise
                    error('unknown dep type!');
            end
            
            % calculate the metric
            metricConditioned = rscdm(y, z, x, 'kendall', rsdm_rectWidths, ...
                rsdm_overlapIncr, rsdm_alpha, rsdm_MA_size, rsdm_peakFinderSel, rsdm_boostFactor);
            metricUnconditioned = copuladeptest10(y, z, 'kendall', rsdm_rectWidths, ...
                rsdm_overlapIncr, rsdm_alpha, rsdm_MA_size, rsdm_peakFinderSel, rsdm_boostFactor);
            
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

% TODO: insert optimal parameters here ... these are not optimal
rsdm_rectWidths = [1 0.75 0.5 0.25 0.125];
rsdm_overlapIncr = 0.055;
rsdm_alpha = 0.03;
rsdm_MA_size = 2;
rsdm_peakFinderSel = 8;
rsdm_boostFactor = 1;

% Simon & Tibshirani use xMin=0, xMax=1 for performing their analysis ...
xMin = 0;
xMax = 1;

testAlpha = 0.05;

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the Type II Error simulation...\n'),'keepthis','timestamp');
num_noise_test_min = 5;
num_noise_test_max = 20;
for l=num_noise_test_min:num_noise_test_max
    for typ=1:numDepTests
        dispstat(sprintf('Computing for noise level=%d Dependency Test=%d',l, typ),'keepthis', 'timestamp');
        % simulate data under the null w/ correct marginals
        for ii=1:nsim
            dispstat(sprintf('Simulating Null -- %0.02f', ii/nsim*100),'timestamp');
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
            metricConditioned = rscdm(y, z, x, 'kendall', rsdm_rectWidths, ...
                rsdm_overlapIncr, rsdm_alpha, rsdm_MA_size, rsdm_peakFinderSel, rsdm_boostFactor);
            metricUnconditioned = copuladeptest10(y, z, 'kendall', rsdm_rectWidths, ...
                rsdm_overlapIncr, rsdm_alpha, rsdm_MA_size, rsdm_peakFinderSel, rsdm_boostFactor);
            
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

% Plot the results

%% Do a conditionally independent test

clear;
clc;

M = 500;
noise = 0.2;

% Generate data from     Y<--X-->Z
x = rand(M,1)*2-1;
y = x + noise*randn(M,1);
z = x.^2 + noise*randn(M,1);

rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% RSCDM conditions X on Y and Z, and sees how related Y and Z
% are to each other ... in this graphical model, they should be UNRELATED
% after the effect of X is removed... i.e. close to independent.  To see
% why, look at the graphical model, Y indep Z | X according to
% d-separation.  So if we condition upon X (i.e. remove teh effect of X on
% Y and Z separately), then we should get independence.
[rscdmVal, Rx_aligned, Ry_aligned] = rscdm(y,z,x);

subplot(3,3,1);
scatter(x,y); grid on; xlabel('x'); ylabel('y'); title(sprintf('RSDM=%0.2f', rsdm1));

subplot(3,3,2);
scatter(x,z); grid on; xlabel('x'); ylabel('z'); title(sprintf('RSDM=%0.2f', rsdm2));

subplot(3,3,3);
scatter(y,z); grid on; xlabel('y'); ylabel('z'); title(sprintf('RSDM=%0.2f', rsdm3));

subplot(3,3,4);
scatter(pobs(x),pobs(y), 'r'); grid on; xlabel('F_X(x)'); ylabel('F_Y(y)');

subplot(3,3,5);
scatter(pobs(x),pobs(z), 'r'); grid on; xlabel('F_X(x)'); ylabel('F_Z(z)');

subplot(3,3,6);
scatter(pobs(y),pobs(z), 'r'); grid on; xlabel('F_Y(y)'); ylabel('F_Z(z)');

subplot(3,3,7);
scatter(Rx_aligned,Ry_aligned); grid on; xlabel('R_x'); ylabel('R_y');  title(sprintf('RSCDM=%0.2f', rscdmVal));

subplot(3,3,8);
scatter(pobs(Rx_aligned),pobs(Ry_aligned), 'r'); grid on; xlabel('F_{R_x}'); ylabel('F_{R_y}');


%% Do a conditionally dependent test

% Tests Conditional Independence w/ RSDM and the residuals processing
% algorithm

clear;
clc;

M = 500;
noise = 0.2;

% Generate data from     Y-->X<--Z
y = rand(M,1);
z = rand(M,1);
x = y.^2+z.^3;


rsdm1 = rsdm(x, y);
rsdm2 = rsdm(x, z);
rsdm3 = rsdm(y, z);
% In this graphical model, Y and Z are independent of each other, but when
% conditioned upon X, they become dependent.  Refer to the rules of
% d-separation to see why, this is a V-Structure!
[rscdmVal, Rx_aligned, Ry_aligned] = rscdm(y,z,x);

subplot(3,3,1);
scatter(x,y); grid on; xlabel('x'); ylabel('y'); title(sprintf('RSDM=%0.2f', rsdm1));

subplot(3,3,2);
scatter(x,z); grid on; xlabel('x'); ylabel('z'); title(sprintf('RSDM=%0.2f', rsdm2));

subplot(3,3,3);
scatter(y,z); grid on; xlabel('y'); ylabel('z'); title(sprintf('RSDM=%0.2f', rsdm3));

subplot(3,3,4);
scatter(pobs(x),pobs(y), 'r'); grid on; xlabel('F_X(x)'); ylabel('F_Y(y)');

subplot(3,3,5);
scatter(pobs(x),pobs(z), 'r'); grid on; xlabel('F_X(x)'); ylabel('F_Z(z)');

subplot(3,3,6);
scatter(pobs(y),pobs(z), 'r'); grid on; xlabel('F_Y(y)'); ylabel('F_Z(z)');

subplot(3,3,7);
scatter(Rx_aligned,Ry_aligned); grid on; xlabel('R_x'); ylabel('R_y');  title(sprintf('RSCDM=%0.2f', rscdmVal));

subplot(3,3,8);
scatter(pobs(Rx_aligned),pobs(Ry_aligned), 'r'); grid on; xlabel('F_{R_x}'); ylabel('F_{R_y}');
