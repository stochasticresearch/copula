% Script which generates multivariate copulas of different types with discrete 
% marginal data, and then reconstructs the copula using beta kernels. 
% It then characterizes the error between the constructed copula and 
%  1.) the actual copula
%  2.) c(F(x1) ... f(Xn)) and f(x1 ... xn) / [ f(x1) * ... * f(xn) ]

clear;
clc;


% define all the parametrics we would like to test over
copulaTypes = {'Gaussian', 'T', 'Frank', 'Gumbel', 'Clayton'};
K_vec = 25:25:100;
M_vec = 250:250:1000;
N_vec = 2:4;
h_vec = [0.01 0.05 0.1 0.2];

% TODO: parametrize these also!
Rho_2D = [1 0.3; 0.3 1];
Rho_3D = [1 0.3 0.6; 0.3 1 0.7; 0.6 0.7 1];
Rho_4D = [1 0.3 0.6 0.4; 0.3 1 0.7 0.1; 0.6 0.7 1 0.2; 0.4 0.1 0.2 1];
nu = 1;
alpha = 6;

% define the marginal distributions to be multinomial distributions, you
% should also parametrize this with different probabilities and for
% different number of discrete outcomes
F_X1 = makedist('Multinomial','Probabilities',[0.2 0.2 0.2 0.2 0.2]);
F_X2 = makedist('Multinomial','Probabilities',[0.5 0.1 0.1 0.1 0.2]);
F_X3 = makedist('Multinomial','Probabilities',[0.9 0.01 0.01 0.04 0.4]);
F_X4 = makedist('Multinomial','Probabilities',[0.2 0.3 0.4 0.05 0.05]);
F_X = {F_X1, F_X2, F_X3, F_X4};

for copulaType=copulaTypes
    for K=K_vec
        for M=M_vec
            for N=N_vec
                for h=h_vec
                    if(N==2)
                        Rho = Rho_2D;
                        [U1,U2] = ndgrid(u);
                    elseif(N==3)
                        Rho = Rho_3D;
                        [U1,U2,U3] = ndgrid(u);
                    elseif(N==4)
                        Rho = Rho_4D;
                        [U1,U2,U3,U4] = ndgrid(u);
                    end
                    % generate the copula random variates to set the
                    % dependency
                    if(strcmp(copulaType, 'Gaussian'))
                        U = copularnd(copulaType, Rho, M);
                    elseif(strcmp(copulaType, 'T'))
                        U = copularnd(copulaType, Rho, nu, M);
                    elseif(strcmp(copulaType, 'Frank'))
                        U = frankcopularnd(M, N, alpha);
                    elseif(strcmp(copulaType, 'Gumbel'))
                        U = gumbelcopularnd(M, N, alpha);
                    elseif(strcmp(copulaType, 'Clayton'))
                        U = claytoncopularnd(M, N, alpha);
                    else
                        error('Unknown copula type!');
                    end
                    
                    % apply icdf function to generate discrete random
                    % variates w/ the defined dependency structure
                    X = zeros(M,N);
                    for nn=1:N
                        X(:,nn) = icdf(F_X{nn}, U(:,nn));
                    end
                    
                    % continue the discrete observations
                    X_continued = continueRv( X );
                    
                    % generate the actual copula density c
                    
                    
                    % estimate the copula density from the continued 
                    % observations - call this c_hat_dobs
                    c_hat_dobs = empcopuladensity(X_continued, h, K, 'betak');
                    
                    % compute error between estimated copula (c_hat_dobs)
                    % density and actual copula density c, call this e_dobs
                    
                    % estimate the copula density w/ the copula random 
                    % variates directly - call this c_hat
                    
                    % compute error between c_hat and c, call this
                    % e_bestcase
                    
                    % compute error between c_hat and c_hat_dobs, call this
                    % e_true
                end
            end
        end
    end
end

K = 25;
u = linspace(0,1,K);
[U1,U2] = ndgrid(u);
c2 = copulapdf(copulaType, [U1(:) U2(:)],5);
c2 = reshape(c2, K,K);
h1 = subplot(1,2,1);
surf(U1,U2,c2);
xlabel('u1')
ylabel('u2')
title('Actual Copula')

alpha = 5;
M = 1000;

X = copularnd(copulaType,alpha,M);

h = .05;
c1 = empcopuladensity(X, h, K, 'betak');
h2 = subplot(1,2,2);
surf(U1,U2,c1);
xlabel('u1')
ylabel('u2')
title('Empirical Copula')

mse = mean((c1(:)-c2(:)).^2);
fprintf('MSE = %f\n', mse);

hlink = linkprop([h1,h2],{'CameraPosition','CameraUpVector'});
rotate3d on