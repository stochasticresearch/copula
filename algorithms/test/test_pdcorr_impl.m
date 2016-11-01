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

%% A script which tests the differences between the R and Matlab implementations
% of pdcorr

clear;
clc;

nsim = 50;
M = 500;
thresh = 0.001;
pval_thresh = 0.01;
numR = 2000;

pdcov_mse_vec = zeros(9, nsim);
pdcov_pval_mse_vec = zeros(9, nsim);

for ii=1:nsim
    % generate uncorrelated random variables
    x = rand(M,1); y = rand(M,1); z = rand(M,1);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(1,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(1,ii) = (pval_r-pval_matlab).^2;
%     [pdcorr_r, pval_r] = pdcorr_R(x,y,z);
%     [pdcorr_matlab, pval_matlab] = pdcorr(x,y,z);
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in uncorrelated test');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in uncorrelated test');
%         pause;
%     end
        
    % generate correlated random variables
    xyz = mvnrnd([0 0 0], [1 0 .3; 0 1 -.1; .3 -.1 1], M);
    x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(2,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(2,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(xyz(:,1), xyz(:,2), xyz(:,3));
%     pdcorr_matlab = pdcorr(xyz(:,1), xyz(:,2), xyz(:,3));
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in correlated test 1');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in correlated test');
%         pause;
%     end
    
    xyz = mvnrnd([0 0 0], [1 0 .2; 0 1 -.8; .2 -.8 1], M);
    x = xyz(:,1); y = xyz(:,2); z = xyz(:,3);
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(3,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(3,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(xyz(:,1), xyz(:,2), xyz(:,3));
%     pdcorr_matlab = pdcorr(xyz(:,1), xyz(:,2), xyz(:,3));
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in correlated test 2');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in correlated test 2');
%         pause;
%     end
    
    % generate ^ structure of various types
    x = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*x + (1-gamma)*eps;
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(4,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(4,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(x,y,z);
%     pdcorr_matlab = pdcorr(x,y,z);
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in ^ test 1');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in ^ test 1');
%         pause;
%     end
    
    x = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*x + (1-gamma)*eps;
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(5,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(5,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(x,y,z);
%     pdcorr_matlab = pdcorr(x,y,z);
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in ^ test 2');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in ^ test 2');
%         pause;
%     end
    
    x = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    y = gamma*x + (1-gamma)*eps;
    z = gamma*x + (1-gamma)*eps;
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(6,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(6,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(x,y,z);
%     pdcorr_matlab = pdcorr(x,y,z);
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in ^ test 3');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in ^ test 3');
%         pause;
%     end
    
    % generate v structure of various types
    y = rand(M,1); z = rand(M,1);
    gamma = 0.2; eps = randn(M,1);
    x = gamma*(y+z) + (1-gamma)*eps;
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(7,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(7,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(x,y,z);
%     pdcorr_matlab = pdcorr(x,y,z);
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in v test 1');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in v test 1');
%         pause;
%     end
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.5; eps = randn(M,1);
    x = gamma*(y+z) + (1-gamma)*eps;
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(8,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(8,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(x,y,z);
%     pdcorr_matlab = pdcorr(x,y,z);
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in v test 2');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in v test 2');
%         pause;
%     end
    
    y = rand(M,1); z = rand(M,1);
    gamma = 0.8; eps = randn(M,1);
    x = gamma*(y+z) + (1-gamma)*eps;
%     [pdcov_r, pval_r] = pdcov_R(x,y,z);
%     [pdcov_matlab, pval_matlab] = pdcov(x,y,z);
    [pdcov_r, pval_r] = pdcov_R(x,y,z, numR);
    [pdcov_matlab, pval_matlab] = pdcov(x,y,z, numR);
    fprintf('pdcov_r=%0.06f pdcov_matlab=%0.06f\n', pdcov_r, pdcov_matlab);
    fprintf('pval_r=%0.06f pval_matlab=%0.06f\n', pval_r, pval_matlab);
    pdcov_mse_vec(9,ii) = (pdcov_r-pdcov_matlab).^2;
    pdcov_pval_mse_vec(9,ii) = (pval_r-pval_matlab).^2;
%     pdcorr_r = pdcorr_R(x,y,z);
%     pdcorr_matlab = pdcorr(x,y,z);
%     if( abs(pdcorr_r-pdcorr_matlab)>=thresh)
%         fprintf('metric fail in v test 3');
%         pause;
%     end
%     if( abs(pval_r - pval_matlab)>=thresh)
%         fprintf('pval fail in v test 3');
%         pause;
%     end
end

pdcov_mse = mean(pdcov_mse_vec,2);
pdcov_pval_mse = mean(pdcov_pval_mse_vec,2);