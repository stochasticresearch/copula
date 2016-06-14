clear;
clc;

% A test of how to do structure learning w/ copulas and optimization

% define a dag as X-->Z-->Y
% generate data from this dag

M = 1000;
Rho = [1 0.8; 0.8 1];

K = 100;
uu = linspace(0,1,K);
[U1,U2] = meshgrid(uu,uu);
c_gaussian_strong = copulapdf('Gaussian', [U1(:) U2(:)], Rho);
c_gaussian_strong = reshape(c_gaussian_strong,K,K);

U = copularnd('Gaussian', Rho, M);
U_dep = depcopularnd(c_gaussian_strong, U(:,2));
X = norminv(U(:,1),0,1);
Z = norminv(U(:,2),0,1);
Y = norminv(U_dep(:,2),0,1);

data = [X Y Z];
dataU = pobs(data);
theta0 = zeros(3,3);     % the variable to optimize
theta0(1,2) = 0.6; theta0(1,3) = 0.2; theta0(2,3) = 0.2;

f = @(x)optimfun(x,dataU);
Aeq = [1 0 0 0 0 0 0 0 0; ...
       0 0 0 0 1 0 0 0 0; ...
       0 0 0 0 0 0 0 0 1]; 
beq = [0 0 0];
A = [0 1 0 1 0 0 0 0 0; ...
     0 0 1 0 0 0 1 0 0; ...
     0 0 0 0 0 1 0 1 0];
b = [1 1 1];

lb = ones(size(theta0))*-0.99;
ub = ones(size(theta0))*0.99;
nonlcon = @nlcon;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(f,theta0,A,b,Aeq,beq,lb,ub,nonlcon,options)