clear;
clc;

M = 1000;
D = 2;

Rho = [1 0.2; 0.2 1];
Z = mvnrnd([0 0], Rho, M);
U = normcdf(Z,0,1);

X = [unidinv(U(:,1),3) ...
     unidinv(U(:,2),4)];

X_xform = X+(rand(size(X))-1);
X_in = [X_xform(:,1) X_xform(:,2)];

K = 200;
[ C, ~, c ] = empcopula(X_in, K);
% subplot(2,2,1); contour(C); title('C empirical'); grid on;
subplot(2,3,1); contour(c{end}); title('c empirical'); grid on;
c_emp_smooth = smoothn(c{end}');
subplot(2,3,2); contour(c_emp_smooth); title('c empirical smooth'); grid on;
cc1 = cumsum(c_emp_smooth,1);
cc2 = cumsum(cc1,2);
subplot(2,3,3); contour(cc2); title('C empirical smooth'); grid on;


% calculate C & c actual
u = linspace(0,1,200);
[u1,u2] = meshgrid(u,u);
C_actual = copulacdf('Gaussian',[u1(:),u2(:)],Rho);
c_actual = copulapdf('Gaussian',[u1(:),u2(:)],Rho);
subplot(2,3,4); contour(reshape(C_actual,200,200)); title('C actual'); grid on;
subplot(2,3,5); contour(reshape(c_actual,200,200)); title('c actual'); grid on;
CC_actual = reshape(C_actual, 200, 200);


figure;
subplot(2,2,1);
contour(reshape(c_actual,200,200)); title('c actual'); grid on;
subplot(2,2,2);
smoothing = [0.3 0.3];
[xg, ~, ~] = gradients_xy(CC_actual, smoothing);
[~, yg, ~] = gradients_xy(xg, smoothing);
contour(yg); grid on; title('Smoothed Gradient of C actual')
subplot(2,2,3);
smoothing = [10 10];
[xg, ~, ~] = gradients_xy(C, smoothing);
[~, yg, ~] = gradients_xy(xg, smoothing);
contour(yg); grid on; title('Smoothed Gradient of C empirical')

% lets try a multivariate kernel density estimate - described in Copula
% Density Estimation book - technique page 9
h = .1;        % kernel bandwidth
rankedObs = tiedrank(X_in)/(M+1);
c_dens_kde_23 = zeros(200,200);
for uu=1:200
    for vv=1:200
        u = uu/200;
        v = vv/200;
        
        sumVal = sum( mvnpdf([ (u-rankedObs(:,1))/h (v-rankedObs(:,2))/h ]) );
        sumVal = sumVal/(M*h.^2);
        c_dens_kde_23(uu,vv) = sumVal;
    end
end
figure;
subplot(3,2,1);
contour(c_dens_kde_23)
title('Page 9 Technique c')

% page 9 technique cdf
subplot(3,2,2);
t1 = cumsum(c_dens_kde_23,1);
t2 = cumsum(t1,2);
contour(t2);
title('Page 9 Technique C')

% lets try a multivariate kernel density estimate - described in Copula
% Density Estimation book - mirror image technique page 14
c_dens_kde_24 = zeros(200,200);
for uu=1:200
    for vv=1:200
        u = uu/200;
        v = vv/200;
        
        K1 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
        K2 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
        K3 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );
        K4 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );
        K5 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v-2+rankedObs(:,2))/h );
        K6 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v-2+rankedObs(:,2))/h );
        K7 = normpdf( (u-2+rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
        K8 = normpdf( (u-2+rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );
        K9 = normpdf( (u-2+rankedObs(:,1))/h ) .* normpdf( (v-2+rankedObs(:,2))/h );
        
        sumVal = sum( K1 + K2 + K3 + K4 + K5 + K6 + K7 + K8 + K9 );
        sumVal = sumVal/(M*h.^2);
        c_dens_kde_24(uu,vv) = sumVal;
    end
end
subplot(3,2,3)
contour(c_dens_kde_24)
title('Page 14 Technique c')

t1 = cumsum(c_dens_kde_24,1);
t2 = cumsum(t1,2);
subplot(3,2,4)
contour(t2);
title('Page 14 Technique C')


% try beta-kernel approach
c_dens_betak = zeros(200,200);
for uu=1:200
    for vv=1:200
        u = uu/200;
        v = vv/200;
        
        K1 = betapdf(X_in(:,1),u/h+1,((1-u)/h) + 1);
        K2 = betapdf(X_in(:,2),v/h+1,((1-v)/h) + 1);
        
        c_dens_betak(uu,vv) = sum(K1.*K2)/(M*h.^2);
    end
end
subplot(3,2,5)
contour(c_dens_betak)
title('Beta Kernel Approach to c')

t1 = cumsum(c_dens_betak,1);
t2 = cumsum(t1,2);
subplot(3,2,6);
contour(t2);
title('Beta Kernel Approach C')

%%
% this dude's version: http://www.mathworks.com/matlabcentral/fileexchange/41187-fast-kernel-density-estimator--multivariate-
pdf = precluster( X_in' ) ;
H = estimateBandwidth( pdf ) ;
kde = constructKDE(pdf.Mu, pdf.w, H, pdf.Cov) ;
kde_construct = zeros(200,200);
for uu=1:200
    for vv=1:200
        u = uu/200;
        v = vv/200;
        sumVal = 0;
        for ii=1:length(kde.pdf.Mu)
            sumVal = sumVal + mvnpdf([u v], kde.pdf.Mu(:,ii)', kde.pdf.Cov{ii});
        end
        kde_construct(uu,vv) = sumVal;
    end
end

figure;
contour(kde_construct)

% kde = constructKDE(pdf.Mu, pdf.w, H, pdf.Cov) ;
% 
% figure;
% % plot the kdes as tabulated distributions
% subplot(1,2,1) ;
% visualizeKDE('kde', kde, 'tabulated', 0) ; 
% I_k = visualizeKDE('kde', kde, 'tabulated', 1) ; 
% title('Tabulated pdf: Kristan bw') ;
% subplot(1,2,2) ; 
% imagesc(I_k); title('Tabulated pdf: Kristan bw') ; 
% axis equal ; axis tight ; colormap gray ;

%%
% In this section, we generate 2-D and 3-D multivariate discrete
% distributions, and then compare the copula ratio's to the actual
% conditional density values to see which smoothing method we should use
% for the HCBN
clear;

% 2-D test
M = 1000;
D = 2;
Rho = [1 0.2; 0.2 1];
Z = mvnrnd([0 0], Rho, M);
U = normcdf(Z,0,1);
X = [unidinv(U(:,1),3) ...
     unidinv(U(:,2),4)];
X_xform = X+(rand(size(X))-1);
X_in = [X_xform(:,1) X_xform(:,2)];     % assume NodeIdx 2 is the parent

[ pdf_X1X2, combos_X1X2 ] = hist_discrete( X );

% create empirical info objects for X1 and X2 for easy querying of the
% distribution function
[F1,x1] = ecdf(X(:,1));
F1 = F1(2:end);
x1 = x1(2:end);
f1 = zeros(1,length(x1));
idx = 1;
for jj=1:length(x1)
    f1(idx) = sum(X(:,1)==x1(jj))/M;
    idx = idx + 1;
end
empInfoX1 = rvEmpiricalInfo(x1,f1,F1);

[F2,x2] = ecdf(X(:,2));
F2 = F2(2:end);
x2 = x2(2:end);
f2 = zeros(1,length(x2));
idx = 1;
for jj=1:length(x2)
    f2(idx) = sum(X(:,2)==x2(jj))/M;
    idx = idx + 1;
end
empInfoX2 = rvEmpiricalInfo(x2,f2,F2);


% Estimate copula density w/ empcoupla, and smoothing
mses = zeros(5,10,5);   % zz1 = K, zz2 = h, zz3 = [mse_smooth mse_p9 mse_p14 mse_p14v2 mse_betak]
zz1 = 1;
zz2 = 1;
% for K=[100,200,250,500]
for K=[25]
    [ C, ~, c ] = empcopula(X_in, K);
    c_emp_smooth = smoothn(c{end});
    % rescale to ensure it is a density
%     c_emp_smooth = c_emp_smooth/sum(c_emp_smooth(:));
%     fprintf('sum(c_emp_smooth)=%f\n', sum(c_emp_smooth(:)));

    % Estimate the copula density w/ p9 technique
    % h = .5;        % kernel bandwidth

    for h=linspace(0.01,10,10)

        rankedObs = tiedrank(X_in)/(M+1);
        c_dens_kde_p9 = zeros(K,K);
        for uu=1:K
            for vv=1:K
                u = uu/K;
                v = vv/K;

                sumVal = sum( mvnpdf([ (u-rankedObs(:,1))/h (v-rankedObs(:,2))/h ]) );
                sumVal = sumVal/(M*h.^2);
                c_dens_kde_p9(uu,vv) = sumVal;
            end
        end
%         fprintf('sum(c_dens_kde_p9)=%f\n', sum(c_dens_kde_p9(:)));

        % Estimate the copula density w/ p14 technique
        c_dens_kde_p14 = zeros(K,K);
        for uu=1:K
            for vv=1:K
                u = uu/K;
                v = vv/K;

                K1 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
                K2 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
                K3 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );
                K4 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );
                K5 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v-2+rankedObs(:,2))/h );
                K6 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v-2+rankedObs(:,2))/h );
                K7 = normpdf( (u-2+rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
                K8 = normpdf( (u-2+rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );
                K9 = normpdf( (u-2+rankedObs(:,1))/h ) .* normpdf( (v-2+rankedObs(:,2))/h );

                sumVal = sum( K1 + K2 + K3 + K4 + K5 + K6 + K7 + K8 + K9 );
                sumVal = sumVal/(M*h.^2);
                c_dens_kde_p14(uu,vv) = sumVal;
            end
        end
%         fprintf('sum(c_dens_kde_p14)=%f\n', sum(c_dens_kde_p14(:)));
        
        c_dens_kde_p14_v2 = zeros(K,K);
        for uu=1:K
            for vv=1:K
                u = uu/K;
                v = vv/K;

                K1 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
                K2 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v-rankedObs(:,2))/h );
                K3 = normpdf( (u-rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );
                K4 = normpdf( (u+rankedObs(:,1))/h ) .* normpdf( (v+rankedObs(:,2))/h );

                sumVal = sum( K1 + K2 + K3 + K4 );
                sumVal = sumVal/(M*h.^2);
                c_dens_kde_p14_v2(uu,vv) = sumVal;
            end
        end
%         fprintf('sum(c_dens_kde_p14v2)=%f\n', sum(c_dens_kde_p14_v2(:)));

        % estimate the copula density w/ beta kernels
        c_dens_betak = empcopdens_betak_v2(rankedObs(:,1), rankedObs(:,2), h, K);
%         fprintf('sum(c_dens_betak)=%f\n', sum(c_dens_betak(:)));
        
        % now compare the actual Rc w/ the expected Rc
        Rc_vals = zeros(6,M);     % row 1 = actual
                                    % row 2 = c_emp_smooth
                                    % row 3 = c_dens_kde_p9
                                    % row 4 = c_dens_kde_p14
                                    % row 5 = c_dens_kde_p14v2
                                    % row 6 = c_dens_betak
        ii = 1;
        for m=1:M
            yy = [X(m,1) X(m,2)];
            numIdx = find(yy(1)==combos_X1X2(1,:) & yy(2)==combos_X1X2(2,:));
            Rc_expect_num = pdf_X1X2(numIdx);
            Rc_expect_den = empInfoX1.queryDensity(X(m,1))*empInfoX2.queryDensity(X(m,2));
            Rc_expect = Rc_expect_num/Rc_expect_den;

            uu = [empInfoX1.queryDistribution(X(m,1)) empInfoX2.queryDistribution(X(m,2))];
            %%% NOTE: all den's will be the same b/c only 1 parent
            % actual copula ratio for c_emp_smooth
            [~,Rc_c_emp_smooth_num] = empcopula_val([], c_emp_smooth, uu);
            Rc_c_emp_smooth = Rc_c_emp_smooth_num/1;

            % actual copula ratio for c_dens_kde_p9
            [~,Rc_c_dens_kde_p9_num] = empcopula_val([], c_dens_kde_p9, uu);
            Rc_c_dens_kde_p9 = Rc_c_dens_kde_p9_num/1;

            % actual copula ratio for c_dens_kde_p14
            [~,Rc_c_dens_kde_p14_num] = empcopula_val([], c_dens_kde_p14, uu);
            Rc_c_dens_kde_p14 = Rc_c_dens_kde_p14_num/1;
            
            % actual copula ratio for c_dens_kde_p14_v2
            [~,Rc_c_dens_kde_p14v2_num] = empcopula_val([], c_dens_kde_p14_v2, uu);
            Rc_c_dens_kde_p14v2 = Rc_c_dens_kde_p14v2_num/1;

            % actual copula ratio for c_dens_betak
            [~,Rc_c_dens_betak_num] = empcopula_val([], c_dens_betak, uu);
            Rc_c_dens_betak = Rc_c_dens_betak_num/1;

            Rc_vals(:,ii) = [Rc_expect Rc_c_emp_smooth Rc_c_dens_kde_p9 Rc_c_dens_kde_p14 Rc_c_dens_kde_p14v2 Rc_c_dens_betak]';
            ii = ii + 1;
        end
%         Rc_vals = log(Rc_vals);
        % subplot(3,2,1); plot(Rc_vals(1,:)); grid on; title('Rc Expect');
        % subplot(3,2,2); plot(Rc_vals(2,:)); grid on; title('Rc Smooth');
        % subplot(3,2,3); plot(Rc_vals(3,:)); grid on; title('Rc P9 Technique');
        % subplot(3,2,4); plot(Rc_vals(4,:)); grid on; title('Rc P14 Technique');
        % subplot(3,2,5); plot(Rc_vals(5,:)); grid on; title('Rc Beta Kernels');

        % Calculate MSE for expect versus different density estimates
        mse_smooth = mean((Rc_vals(1,:)-Rc_vals(2,:)).^2);
        mse_p9 = mean((Rc_vals(1,:)-Rc_vals(3,:)).^2);
        mse_p14 = mean((Rc_vals(1,:)-Rc_vals(4,:)).^2);
        mse_p14v2 = mean((Rc_vals(1,:)-Rc_vals(5,:)).^2);
        mse_betak = mean((Rc_vals(1,:)-Rc_vals(6,:)).^2);
        
        mses(zz1,zz2,1) = mse_smooth;
        mses(zz1,zz2,2) = mse_p9;
        mses(zz1,zz2,3) = mse_p14;
        mses(zz1,zz2,4) = mse_p14v2;
        mses(zz1,zz2,5) = mse_betak;
        
        fprintf('2D-->Kernel K=%d BW=%f MSE_SMOOTH=%f MSE_P9=%f MSE_P14=%f MSE_P14_V2=%f MSE_BETAK=%f\n', ...
            K, h, mse_smooth, mse_p9, mse_p14, mse_p14v2, mse_betak);
        zz2 = zz2 + 1;
    end
    zz1 = zz1 + 1;
end

save('mses_2d.mat', 'mses')

%%
% In this section, we generate 3-D multivariate discrete
% distributions, and then compare the copula ratio's to the actual
% conditional density values to see which smoothing method we should use
% for the HCBN
clear;

% 3-D test
M = 1000;
D = 3;
Rho = [1 .4 .2; .4 1 -.8; .2 -.8 1];
Z = mvnrnd([0 0 0], Rho, M);
U = normcdf(Z,0,1);
X = [unidinv(U(:,1),3) ...
     unidinv(U(:,2),4) ...
     unidinv(U(:,2),5)];
X_xform = X+(rand(size(X))-1);
X_in = [X_xform(:,1) X_xform(:,2) X_xform(:,3)];     % assume NodeIdx 2,3 are the parents

[ pdf_X1X2X3, combos_X1X2X3 ] = hist_discrete( X );
[ pdf_X2X3, combos_X2X3 ] = hist_discrete( [X(:,2) X(:,3)] );

% create empirical info objects for X1 and X2 for easy querying of the
% distribution function
[F1,x1] = ecdf(X(:,1));
F1 = F1(2:end);
x1 = x1(2:end);
f1 = zeros(1,length(x1));
idx = 1;
for jj=1:length(x1)
    f1(idx) = sum(X(:,1)==x1(jj))/M;
    idx = idx + 1;
end
empInfoX1 = rvEmpiricalInfo(x1,f1,F1);

[F2,x2] = ecdf(X(:,2));
F2 = F2(2:end);
x2 = x2(2:end);
f2 = zeros(1,length(x2));
idx = 1;
for jj=1:length(x2)
    f2(idx) = sum(X(:,2)==x2(jj))/M;
    idx = idx + 1;
end
empInfoX2 = rvEmpiricalInfo(x2,f2,F2);

[F3,x3] = ecdf(X(:,3));
F3 = F3(2:end);
x3 = x3(2:end);
f3 = zeros(1,length(x3));
idx = 1;
for jj=1:length(x3)
    f3(idx) = sum(X(:,3)==x3(jj))/M;
    idx = idx + 1;
end
empInfoX3 = rvEmpiricalInfo(x3,f3,F3);

% for K=[100,200,250,500]
for K=[25]
    [ C, ~, c ] = empcopula(X_in, K);
    c_emp_smooth_num = smoothn(c{end});
    c_emp_smooth_num = c_emp_smooth_num/sum(c_emp_smooth_num(:));
    [ C_den, ~, c_den ] = empcopula(X_in(:,2:3), K);
    c_emp_smooth_den = smoothn(c_den{end});
    c_emp_smooth_den = c_emp_smooth_den/sum(c_emp_smooth_den(:));

    rankedObs = tiedrank(X_in)/(M+1);
    for h=linspace(0.01,10,10);
        c_dens_betak_num = empcopdens_betak_3d_v2(rankedObs(:,1), rankedObs(:,2), rankedObs(:,3), h, K);
        c_dens_betak_den = empcopdens_betak_v2(rankedObs(:,2), rankedObs(:,3), h, K);

        % now compare the actual Rc w/ the expected Rc
        Rc_vals = zeros(2,M);     % row 1 = actual
                                    % row 2 = c_dens_betak
        ii = 1;
        for m=1:M
            yy = [X(m,1) X(m,2) X(m,3)];
            numIdx = find(yy(1)==combos_X1X2X3(1,:) & yy(2)==combos_X1X2X3(2,:) & yy(3)==combos_X1X2X3(3,:));
            Rc_expect_num = pdf_X1X2X3(numIdx);
            den1Idx = find(yy(2)==combos_X2X3(1,:) & yy(3)==combos_X2X3(2,:));
            Rc_expect_den = empInfoX1.queryDensity(X(m,1))*pdf_X2X3(den1Idx);
            Rc_expect = Rc_expect_num/Rc_expect_den;

            uu_num = [empInfoX1.queryDistribution(X(m,1)) empInfoX2.queryDistribution(X(m,2)) empInfoX2.queryDistribution(X(m,3))];
            uu_den = [empInfoX2.queryDistribution(X(m,2)) empInfoX2.queryDistribution(X(m,3))];
            
            % actual copula ratio for c_dens_betak
            [~,Rc_c_dens_betak_num] = empcopula_val([], c_dens_betak_num, uu_num);
            [~,Rc_c_dens_betak_den] = empcopula_val([], c_dens_betak_den, uu_den);
            Rc_c_dens_betak = Rc_c_dens_betak_num/Rc_c_dens_betak_den;

            Rc_vals(:,ii) = [Rc_expect Rc_c_dens_betak]';
            ii = ii + 1;
        end
%         Rc_vals = log(Rc_vals);

        % Calculate MSE for expect versus different density estimates
        mse_betak = mean((Rc_vals(1,:)-Rc_vals(2,:)).^2);
        
        fprintf('3D-->Kernel K=%d BW=%f MSE_BETAK=%f\n', ...
            K, h, mse_betak);
        
    end
end