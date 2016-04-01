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

% A script to understand which functions are the best way to estimate the
% PDF and CDF for discrete and continuous random variables

%% Discrete RV Test
clear;
clc;

pd = makedist('Multinomial', 'probabilities',[0.1 0.2 0.5 0.2]);
% generate samples from this distribution
x_data = random(pd, 250, 1);

[F_ecdf, x_ecdf] = empcdf(x_data, 1);
f_ecdf = emppdf(x_data, 1);

[f_ksdensity, x_ksdensity] = emppdf(x_data, 0);
F_ksdensity = empcdf(x_data, 0);

subplot(2,2,1); plot(x_ecdf,F_ecdf); grid on; title('ECDF');
subplot(2,2,2); plot(x_ecdf,f_ecdf); grid on; title('PDF - Manual');
subplot(2,2,3); plot(x_ksdensity, F_ksdensity); grid on; title('CDF - ksdensity');
subplot(2,2,4); plot(x_ksdensity, f_ksdensity); grid on; title('PDF - ksdensity');

%% Continuous RV Test
clear;
clc;

x_data = normrnd(0, 0.5, 250, 1);

[F_ecdf, x_ecdf] = empcdf(x_data, 1);
f_ecdf = emppdf(x_data, 1);

fprintf('f_ecdf integral = %f\n', trapz(x_data, f_ecdf));

[f_ksdensity, x_ksdensity] = emppdf(x_data, 0);
F_ksdensity = empcdf(x_data, 0);

subplot(2,2,1); plot(x_ecdf,F_ecdf); grid on; title('ECDF');
subplot(2,2,2); plot(x_ecdf,f_ecdf); grid on; title('PDF - Manual');
subplot(2,2,3); plot(x_ksdensity, F_ksdensity); grid on; title('CDF - ksdensity');
subplot(2,2,4); plot(x_ksdensity, f_ksdensity); grid on; title('PDF - ksdensity');