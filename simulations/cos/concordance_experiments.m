%******************************************************************************
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

% Script explores what different concordance measures do for bivariate
% generated distributions with different dependencies.

clear;
clc;

% Generate the data
range = 4;
n = 1000;
x = -range + 2*range*rand(n,1);       % generate uniform RV's

figtitlestr = sprintf('X ~ Uniform(-%d,%d)',range,range);

% affine dependence
fprintf('Processing Affine Dependence\n');
m = -20:.1:20;
ktau_affine = zeros(1,length(m));
srho_affine = zeros(1,length(m));
ii = 1;
for mm=m
    y = mm*x + 1;
    
    % calculate kendall's tau
    ktau = corr(x,y,'type','Kendall');
    ktau_affine(ii) = ktau;
    
    % calculate spearman's rho
    srho = corr(x,y,'type','Spearman');
    srho_affine(ii) = srho;
        
    ii = ii + 1;
end

figure;
subplot(2,3,1);
plot(m,ktau_affine,'b',m,srho_affine,'ro');
grid on
xlabel('\alpha')
title('Y = \alpha * X + 1')
legend('Kendall \tau', 'Spearman \rho')

% quadratic dependence
fprintf('Processing Quadratic Dependence\n');
ktau_quadratic = zeros(1,length(m));
srho_quadratic = zeros(1,length(m));
ii = 1;
for mm=m
    y = mm*x.^2;
    
    % calculate kendall's tau
    ktau = corr(x,y,'type','Kendall');
    ktau_quadratic(ii) = ktau;
    
    % calculate spearman's rho
    srho = corr(x,y,'type','Spearman');
    srho_quadratic(ii) = srho;
        
    ii = ii + 1;
end

subplot(2,3,2);
plot(m,ktau_quadratic,'b',m,srho_quadratic,'ro');
grid on
xlabel('\alpha')
title('Y = \alpha * X^2')
legend('Kendall \tau', 'Spearman \rho')

% 4th order Poly dependence
fprintf('Processing 4th order polynomial dependence\n');
ktau_fourthOrder = zeros(1,length(m));
srho_fourthOrder = zeros(1,length(m));
ii = 1;
for mm=m
    y = mm*(x.^4-2*x.^2);
    
    % calculate kendall's tau
    ktau = corr(x,y,'type','Kendall');
    ktau_fourthOrder(ii) = ktau;
    
    % calculate spearman's rho
    srho = corr(x,y,'type','Spearman');
    srho_fourthOrder(ii) = srho;
        
    ii = ii + 1;
end

subplot(2,3,3);
plot(m,ktau_fourthOrder,'b',m,srho_fourthOrder,'ro');
grid on
xlabel('\alpha')
title('Y=\alpha * (X^4-2X^2)')
legend('Kendall \tau', 'Spearman \rho')

% sinusoidal dependence
fprintf('Processing Sinusoidal Dependence\n');

ktau_sinusoidal_amplitude = zeros(1,length(m));
srho_sinusoidal_amplitude = zeros(1,length(m));
ktau_sinusoidal_frequency = zeros(1,length(m));
srho_sinusoidal_frequency = zeros(1,length(m));
ktau_sinusoidal_phase = zeros(1,length(m));
srho_sinusoidal_phase = zeros(1,length(m));

ii = 1;
for mm=m
    y = mm.*sin(x);
    % calculate kendall's tau
    ktau = corr(x,y,'type','Kendall');
    ktau_sinusoidal_amplitude(ii) = ktau;
    % calculate spearman's rho
    srho = corr(x,y,'type','Spearman');
    srho_sinusoidal_amplitude(ii) = srho;
    
    y = sin(mm*x);
    % calculate kendall's tau
    ktau = corr(x,y,'type','Kendall');
    ktau_sinusoidal_frequency(ii) = ktau;
    % calculate spearman's rho
    srho = corr(x,y,'type','Spearman');
    srho_sinusoidal_frequency(ii) = srho;
    
    y = sin(x+mm);
    % calculate kendall's tau
    ktau = corr(x,y,'type','Kendall');
    ktau_sinusoidal_phase(ii) = ktau;
    % calculate spearman's rho
    srho = corr(x,y,'type','Spearman');
    srho_sinusoidal_phase(ii) = srho;
    
    ii = ii + 1;
end

subplot(2,3,4);
plot(m,ktau_sinusoidal_amplitude,'b',m,srho_sinusoidal_amplitude,'ro');
grid on
xlabel('\alpha')
title('Y=\alpha * sin(X)')
legend('Kendall \tau', 'Spearman \rho')

subplot(2,3,5);
plot(m,ktau_sinusoidal_frequency,'b',m,srho_sinusoidal_frequency,'ro');
grid on
xlabel('\alpha')
title('Y=sin( \alpha * X)')
legend('Kendall \tau', 'Spearman \rho')

subplot(2,3,6);
plot(m,ktau_sinusoidal_phase,'b',m,srho_sinusoidal_phase,'ro');
grid on
xlabel('\alpha')
title('Y=sin(X + \alpha)')
legend('Kendall \tau', 'Spearman \rho')

h = figtitle(figtitlestr);