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

%% Single Modal Gaussian Distribution test - low variance
clear;
clc;

% Test script for estMteDensity
numPts = 1000;
x = normrnd(0,0.3,numPts,1);
[f,xi] = ksdensity(x);
mte_estimate = estMteDensity(x);

plot(xi,f,mte_estimate.domain,mte_estimate.density);
grid on;
xlabel('x')
legend('KDE', 'MTE')

refLikelihood = -1*normlike([0,1],x);

mteLLVal = 0;
for ii=1:numPts
    mteLLVal = mteLLVal + log(mte_estimate.pdf(x(ii)));
end
title(sprintf('Ref LL=%0.02f, MTE LL=%0.02f', refLikelihood, mteLLVal));

%% Single Modal Uniform Distribution test
clear;
clc;

% Test script for estMteDensity
numPts = 1000;
x = unifrnd(0,5,numPts,1);
[f,xi] = ksdensity(x);
mte_estimate = estMteDensity(x);

plot(xi,f,mte_estimate.domain,mte_estimate.density);
grid on;
xlabel('x')
legend('KDE', 'MTE')

mteLLVal = 0;
for ii=1:numPts
    mteLLVal = mteLLVal + log(mte_estimate.pdf(x(ii)));
end
title(sprintf('MTE LL=%0.02f', mteLLVal));

%% Multimodal Gaussian Distrubtion Test
clear;
clc;

numPts = 100;
x = [normrnd(-2,0.3,numPts,1); normrnd(2,0.5,numPts,1)];
[f,xi] = ksdensity(x);
mte_estimate = estMteDensity(x);

plot(xi,f,mte_estimate.domain,mte_estimate.density);
grid on;
xlabel('x')
legend('KDE', 'MTE')

mteLLVal = 0;
for ii=1:numPts
    mteLLVal = mteLLVal + log(mte_estimate.pdf(x(ii)));
end
title(sprintf('MTE LL=%0.02f', mteLLVal));
