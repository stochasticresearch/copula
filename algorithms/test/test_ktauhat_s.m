%**********************************************************************
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
%* along with this program.  If not, see <http://www.gnu.org/licenses/>
%* 
%**********************************************************************

%% Test the ktauhat "streaming" mode
clear;
clc;
dbstop if error;

M = 100; numDiscreteIntervals = 4;

% x = [1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9];
% y = [5.1 4.2 3.3 2.4 1.5 2.6 3.7 4.8 5.9];

% x = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25];
% y = [5 4 3 2 1 2 3 4 5 4  3  2  1  2  3  4  5  4  3  2  1  2  3  4  5];

% x = [1 2 3 3 4 4 5 5 6 7 8];
% y = [1 1 1 2 1 2 1 2 2 2 2];

% x = [1 1 1 1 2 2 2 2];
% y = [1 1 1 1 2 2 2 2];

rng(1234);

testContinuous = 1;
testHybrid1 = 1;
testHybrid2 = 1;
testDiscrete = 1;

tol = 0.02;
numTests = 9;
CORRECTION_FACTOR_SETTING = 5;      % this corresponds to NO correction factor
for testNum=1:numTests
    x = rand(M,1);
    x = sort(x);
    x_discrete = discretizeRv(x,numDiscreteIntervals)';
    
    switch(testNum)
        case 1
            y = x;
            y_discrete = x_discrete;
        case 2
            y = (x-0.5).^2;
            y_discrete = (x_discrete-0.5).^2;
        case 3
            y = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3);
            y_discrete = 128*(x_discrete-1/3).^3-48*(x_discrete-1/3).^3-12*(x_discrete-1/3);
        case 4
            y = sin(4*pi*x);
            y_discrete = sin(4*pi*x_discrete);
        case 5
            y = sin(16*pi*x);
            y_discrete = sin(16*pi*x_discrete);
        case 6
            y = x.^(1/4);
            y_discrete = x_discrete.^(1/4);
        case 7
            y=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2));
            y_discrete=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x_discrete - 1).^2));
        case 8
            y = (x > 0.5);
            y_discrete = (x_discrete > 0.5);
        case 9
            y = rand(M,1);
            y_discrete = discretizeRv(y,numDiscreteIntervals)';
    end
    y_hybrid2 = discretizeRv(y,numDiscreteIntervals)';
    
    y_hybrid1 = sortSubblocks(x_discrete, y);
    y_hybrid2 = sortSubblocks(x, y_hybrid2);
    y_discrete = sortSubblocks(x_discrete, y_discrete);

    fprintf('***** Testing DepType = %d *****\n', testNum);
    
    % test equality for for all continuous/hybrid/all discrete cases
    kso_continuous = ktauhat_s(x, y);
    kso_hybrid1 = ktauhat_s(x_discrete,y_hybrid1);
    kso_hybrid2 = ktauhat_s(x,y_hybrid2);
    kso_discrete = ktauhat_s(x_discrete,y_discrete);
    for ii=2:M
%         fprintf('*******************************\n');
        x_continuous_subset = x(1:ii)';
        x_discrete_subset = x_discrete(1:ii);
        y_continuous_subset = y(1:ii);
        y_hybrid1_subset = y_hybrid1(1:ii);
        y_hybrid2_subset = y_hybrid2(1:ii);
        y_discrete_subset = y_discrete(1:ii);
        
        
        % sort the hybrid data and the discrete data such that ktauhat and
        % ktauhat_s *should* match-up
        
        if(testContinuous)
            tau1_c = ktauhat(x_continuous_subset, y_continuous_subset, CORRECTION_FACTOR_SETTING);
            tau2_c = kso_continuous.consume(1);
            
            if(abs(tau1_c-tau2_c)>tol)
                warning('Continuous Error: ii=%d', ii);
                warning('\t tau1=%0.04f tau2=%0.04f\n', tau1_c, tau2_c);
            end
        end
        
        if(testHybrid1)
            tau1_h1 = ktauhat(x_discrete_subset, y_hybrid1_subset, CORRECTION_FACTOR_SETTING);
            tau2_h1 = kso_hybrid1.consume(1);
            if(abs(tau1_h1-tau2_h1)>tol)
                warning('Hybrid-1 Error: ii=%d', ii);
                warning('\t tau1=%0.04f tau2=%0.04f\n', tau1_h1, tau2_h1);
            end
        end
        
        if(testHybrid2)
            tau1_h2 = ktauhat(x_continuous_subset, y_hybrid2_subset, CORRECTION_FACTOR_SETTING);
            tau2_h2 = kso_hybrid2.consume(1);
            if(abs(tau1_h2-tau2_h2)>tol)
                warning('Hybrid-2 Error: ii=%d', ii);
                warning('\t tau1=%0.04f tau2=%0.04f\n', tau1_h2, tau2_h2);
            end
        end
        
        if(testDiscrete)      % circular isn't correct for discrete
            tau1_d = ktauhat(x_discrete_subset, y_discrete_subset, CORRECTION_FACTOR_SETTING);
            tau2_d = kso_discrete.consume(1);
            if(abs(tau1_d-tau2_d)>tol)
                warning('Discrete Error: ii=%d', ii);
                warning('\t tau1=%0.04f tau2=%0.04f\n', tau1_d, tau2_d);
            end
        end
        
    end
%     fprintf('*******************************\n');
end

%% test ktauhat_s restart mode
clear;
clc;
dbstop if error;

M = 500;
rng(1234);

x = rand(M,1);
x = sort(x);
y = (x-0.5).^2;

tol = 1e-3;

kso = ktauhat_s(x, y);
iiBegin = 1;
iiEnd = M/5;
for ii=iiBegin+1:iiEnd
    x_continuous_subset = x(iiBegin:ii);
    y_continuous_subset = y(iiBegin:ii);
    
    tau1_c = ktauhat(x_continuous_subset, y_continuous_subset);
    tau2_c = kso.consume(1);

    if(abs(tau1_c-tau2_c)>tol)
        warning('Continuous Error: ii=%d', ii);
        warning('\t tau1=%0.04f tau2=%0.04f\n', tau1_c, tau2_c);
    end
end

% perform a reset state
kso.clearState();
iiBegin = iiEnd+1;
iiEnd = M;
for ii=iiBegin+1:iiEnd
    x_continuous_subset = x(iiBegin:ii);
    y_continuous_subset = y(iiBegin:ii);
    
    tau1_c = ktauhat(x_continuous_subset, y_continuous_subset);
    tau2_c = kso.consume(1);

    if(abs(tau1_c-tau2_c)>tol)
        warning('Continuous Error: ii=%d', ii);
        warning('\t tau1=%0.04f tau2=%0.04f\n', tau1_c, tau2_c);
    end
end
