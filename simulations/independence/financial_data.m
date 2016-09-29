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

% Perform bivariate and multivarate statistical dependency tests on
% finanical data.  

clear;
clc;

%% Load all the data and preprocess necessary items

% % % load Data_GlobalIdx1;   % Import daily index closings -- this contains the global indices
% % %                         % Date Ranges = Business days, April 27, 1993 to July 14, 2003  


% now, we use yahoo data to fetch the individual closing prices of each
% constituent within the indices
tickers = {'TSX', 'CAC', 'DAX', 'NIK', 'FTSE', 'SP'};

startDate = '01052013';
endDate = '01062013';

% tsxConstituentData = hist_stock_data(startDate, endDate, 'tsx_tickers.txt');

% we get 40/40 for CAC40.  cac40_tickers.txt is verified to be correct.
cac40ConstituentData = hist_stock_data(startDate, endDate, 'cac40_tickers.txt');      % NOT WORKING :(

% we get 29/30 for DAX.  dax_tickers.txt is verified to be correct.
% daxConstituentData = hist_stock_data(startDate, endDate, 'dax_tickers.txt');

% we get 93/100 for FTSE.  ftse100_tickers.txt is verified to be correct
% ftse100ConstituentData = hist_stock_data(startDate, endDate, 'ftse100_tickers.txt');

% we get 488/500 for S&P500.  sp500_tickers.txt is verified to be correct
% sp500ConstituentData = hist_stock_data(startDate, endDate, 'sp500_tickers.txt');