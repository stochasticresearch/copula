clear;
clc;

D = 2;
n = 1000;
% X = randn(D,n);

rho = .7;
Z = mvnrnd([0 0], [1 rho; rho 1], n);
U = normcdf(Z);
X = [gaminv(U(:,1),2,1) tinv(U(:,2),5)];

X = X';
U = U';
Z = Z';

for d=1:D
   [Z_sorted(d,:) index_alt(d,:)] = sort(X(d,:));
end 

index_alt;

for i=1:n
    for d=1:D        
        nU(d,index_alt(d,i)) =  i; 
        % values n*U(d,index_alt) in order to omit rounding errors
    end                
end
nU = nU/n;

subplot(2,2,1)
scatter(X(1,:),X(2,:))
title('Samples of F_X_Y (Generated)')
grid on
xlabel('X')
ylabel('Y')

subplot(2,2,2)
scatter(U(1,:),U(2,:))
title('Samples of U (Generated)')
grid on
xlabel('U_1')
ylabel('U_2')

subplot(2,2,3)
scatter(index_alt(1,:),index_alt(2,:))
title('Samples of Rank(F_X_Y) (Calculated)')
grid on
xlabel('X')
ylabel('Y')

subplot(2,2,4)
scatter(nU(1,:),nU(2,:))
title('Samples of U_1_2 (Calculated)')
grid on
xlabel('U_1')
ylabel('U_2')
