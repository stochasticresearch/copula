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

%% 2-D test for empcopulaval

clear;
clc;

fid = fopen('/home/kiran/ownCloud/PhD/sim_results/test_empcopulaval.diary', 'w');

% generate samples from Frank copula pdf
alpha = 4;
M = 1000;
D = 2;
K = 25;

numMCSims = 10;
mse_empcopulaval_query_cellarr = cell(1,numMCSims);
mse_est_cellarr = cell(1,numMCSims);

dispstat('','init'); % One time only initialization
dispstat(sprintf('Begining the simulation...'),'keepthis','timestamp');
hVec = .01:.03:.1;
Kvec = 25:25:100;
progressIdx = 1;
for h=hVec
    for K=Kvec
        for mcSimNum=1:numMCSims
            
            progress = progressIdx*100/(length(hVec)*length(Kvec)*numMCSims);
            progressStr = sprintf('Progress: h=%0.02f K=%d %0.02f%%',h,K,progress);
            dispstat(progressStr,'timestamp');
            fprintf(fid, '%s\n', progressStr);

            % generate data for the simulation
            U = frankcopularnd(M, D, alpha);
            h = 0.05;

            % estimate copula pdf empirically
            c_est = empcopulapdf(U, h, K, 'betak');

            multiplyFactor=2;
            u = linspace(.001,.999,K*multiplyFactor);
            [U1,U2] = ndgrid(u);
            UU = [U1(:) U2(:)];
            c_est_query = zeros(size(c_est)*multiplyFactor);
            c_actual_query = zeros(size(c_est)*multiplyFactor);
            for ii=1:size(UU,1)
                uu = UU(ii,:);
                c_est_query(ii) = empcopula_val(c_est, uu);         % because matlab stores column wise, this works as linear indexing
                c_actual_query(ii) = frankcopulapdf(uu, alpha);     % because matlab stores column wise, this works as linear indexing
            end

            se1 = (c_actual_query-c_est_query).^2;
            mse_empcopulaval_query_cellarr{mcSimNum} = se1;
            
            progressIdx = progressIdx + 1;
        end

        mse_empcopulaval_query_mat = mse_empcopulaval_query_cellarr{1};
        for ii=2:numMCSims
            mse_empcopulaval_query_mat = mse_empcopulaval_query_mat + mse_empcopulaval_query_cellarr{ii};
        end
        mse_empcopulaval_query_ave = mse_empcopulaval_query_mat/numMCSims;
        fig1 = figure;
        h1 = subplot(1,2,1); surf(U1,U2,mse_empcopulaval_query_ave); xlabel('u_1'); ylabel('u_2'); 
        grid on; 
        title(sprintf('\\mu MSE - K=%d h=%0.02f MC=%d', K, h, numMCSims));
        % calculate the variance
        var_empcopulaval_query = (mse_empcopulaval_query_cellarr{1}-mse_empcopulaval_query_mat).^2;
        for ii=2:numMCSims
            var_empcopulaval_query = var_empcopulaval_query + (mse_empcopulaval_query_cellarr{ii}-mse_empcopulaval_query_mat).^2;
        end
        var_empcopulaval_query_ave = var_empcopulaval_query/numMCSims;
        h2 = subplot(1,2,2); surf(U1,U2,var_empcopulaval_query_ave); xlabel('u_1'); ylabel('u_2'); 
        grid on; 
        title(sprintf('\\sigma^2 MSE - K=%d h=%0.02f MC=%d', K, h, numMCSims));

        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        print(sprintf('/home/kiran/ownCloud/PhD/sim_results/empcopulaval_K_%d_h_%0.02f', K, h),'-dpng')
        close(fig1);
    end
end
dispstat('Finished.','keepprev');
fclose(fid);