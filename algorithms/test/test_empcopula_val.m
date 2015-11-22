% tests the empcopula_val function

clear;
clc;

error_tol = 0.02;

n = 1000;
rho = .7;
Z = mvnrnd([0 0], [1 rho; rho 1], n);
U_real = normcdf(Z);
X_real = [gaminv(U_real(:,1),2,1) tinv(U_real(:,2),5)];

K = 200;
[C,U,c] = empcopula(X_real,K);    
cc = c{end};
% reconstruct C and c from empcopula_val function, and compare the plots
% and errors between the two to verify accuracy of empcopula_val function
for ii=1:K
    for jj=1:K
        U_1 = U(ii,jj,1); U_2 = U(ii,jj,2);
        U_query = [U_1 U_2]';
        [C_u, c_u] = empcopula_val(C,cc,U_query);
        if( abs(C(ii,jj)-C_u) > error_tol )  
            fprintf('***********************************************\n');
            fprintf('Mismatch: U_1=%f U_2=%f, C(U_1,U_2)=%f C_u=%f\n', U_1,U_2,C(ii,jj),C_u);
        end
        if( abs(cc(ii,jj)-c_u) > error_tol )
            fprintf('Mismatch: U_1=%f U_2=%f, c(U_1,U_2)=%f c_u=%f\n', ...
                U_1,U_2,cc(ii,jj),c_u);
            fprintf('***********************************************\n');
        end
    end
end
