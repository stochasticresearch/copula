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

% shows how to use the genSynthData function
clear;
clc;

discreteType = {};
nodeA = [0.4 0.3 0.2 0.1]; discreteType{1} = nodeA;
nodeB = [0.6 0.1 0.05 0.25]; discreteType{2} = nodeB;
continuousType = 'Gaussian';
% continuousType = 'other';
M = 1000;
X = genSynthData(discreteType, continuousType, M);

% make a pairs style plot
subplot(5,5,1); hist(X(:,1)); xlabel('X_1'); grid on;
subplot(5,5,2); scatter(X(:,1),X(:,2)); xlabel('X_1'); ylabel('X_2'); grid on;
subplot(5,5,3); scatter(X(:,1),X(:,3)); xlabel('X_1'); ylabel('X_3'); grid on;
subplot(5,5,4); scatter(X(:,1),X(:,4)); xlabel('X_1'); ylabel('X_4'); grid on;
subplot(5,5,5); scatter(X(:,1),X(:,5)); xlabel('X_1'); ylabel('X_5'); grid on;

subplot(5,5,6); scatter(X(:,2),X(:,1)); xlabel('X_2'); ylabel('X_1'); grid on;
subplot(5,5,7); hist(X(:,2)); xlabel('X_2'); grid on;
subplot(5,5,8); scatter(X(:,2),X(:,3)); xlabel('X_2'); ylabel('X_3'); grid on;
subplot(5,5,9); scatter(X(:,2),X(:,4)); xlabel('X_2'); ylabel('X_4'); grid on;
subplot(5,5,10); scatter(X(:,2),X(:,5)); xlabel('X_2'); ylabel('X_5'); grid on;

subplot(5,5,11); scatter(X(:,3),X(:,1)); xlabel('X_3'); ylabel('X_1'); grid on;
subplot(5,5,12); scatter(X(:,3),X(:,2)); xlabel('X_3'); ylabel('X_2'); grid on;
subplot(5,5,13); hist(X(:,3)); xlabel('X_3'); grid on;
subplot(5,5,14); scatter(X(:,3),X(:,4)); xlabel('X_3'); ylabel('X_4'); grid on;
subplot(5,5,15); scatter(X(:,3),X(:,5)); xlabel('X_3'); ylabel('X_5'); grid on;

subplot(5,5,16); scatter(X(:,4),X(:,1)); xlabel('X_4'); ylabel('X_1'); grid on;
subplot(5,5,17); scatter(X(:,4),X(:,2)); xlabel('X_4'); ylabel('X_2'); grid on;
subplot(5,5,18); scatter(X(:,4),X(:,3)); xlabel('X_4'); ylabel('X_3'); grid on;
subplot(5,5,19); hist(X(:,4)); xlabel('X_4'); grid on;
subplot(5,5,20); scatter(X(:,4),X(:,5)); xlabel('X_4'); ylabel('X_5'); grid on;

subplot(5,5,21); scatter(X(:,5),X(:,1)); xlabel('X_5'); ylabel('X_1'); grid on;
subplot(5,5,22); scatter(X(:,5),X(:,2)); xlabel('X_5'); ylabel('X_2'); grid on;
subplot(5,5,23); scatter(X(:,5),X(:,3)); xlabel('X_5'); ylabel('X_3'); grid on;
subplot(5,5,24); scatter(X(:,5),X(:,4)); xlabel('X_5'); ylabel('X_4'); grid on;
subplot(5,5,25); hist(X(:,5)); xlabel('X_5'); grid on;