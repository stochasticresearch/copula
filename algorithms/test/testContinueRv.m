% script which tests the continueRv function

M = 1000;
D = 2;
X_in = unidrnd(5,M,D)-1;
X_out = continueRv(X_in);

subplot(1,2,1);
scatter(X_in(:,1), X_in(:,2));
grid on;
title('Original Discrete RVs')
subplot(1,2,2);
scatter(X_out(:,1), X_out(:,2));
grid on;
title('Continued Discrete RVs')