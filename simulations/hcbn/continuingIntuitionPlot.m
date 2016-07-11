probs = [0.5 0.5];
a_dist = makedist('Multinomial','Probabilities',probs);
alpha = 10;
M = 1000;

U = copularnd('Frank', alpha, M);

X_hybrid = zeros(M,2);
X_hybrid(:,1) = a_dist.icdf(U(:,1));
xContinuous = unifrnd(-2,2,M,1);
isdiscrete = 0;
[fContinous,xiContinuous] = emppdf(xContinuous,isdiscrete);
FContinuous = empcdf(xContinuous,isdiscrete);
continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous,isdiscrete);

for ii=1:M
    X_hybrid(ii,2) = continuousDistInfo.icdf(U(ii,2));
end

U_in = pobs(X_hybrid);
X_continued(:,1) = continueRv(X_hybrid(:,1));
X_in = [X_continued(:,1) X_hybrid(:,2)];
U_in2 = pobs(X_in);

probs = 0.1*ones(1,10);
a_dist = makedist('Multinomial','Probabilities',probs);
alpha = 10;
M = 1000;

U = copularnd('Frank', alpha, M);

X_hybrid = zeros(M,2);
X_hybrid(:,1) = a_dist.icdf(U(:,1));
xContinuous = unifrnd(-2,2,M,1);
isdiscrete = 0;
[fContinous,xiContinuous] = emppdf(xContinuous,isdiscrete);
FContinuous = empcdf(xContinuous,isdiscrete);
continuousDistInfo = rvEmpiricalInfo(xiContinuous,fContinous,FContinuous,isdiscrete);

for ii=1:M
    X_hybrid(ii,2) = continuousDistInfo.icdf(U(ii,2));
end

U_in = pobs(X_hybrid);
X_continued(:,1) = continueRv(X_hybrid(:,1));
X_in = [X_continued(:,1) X_hybrid(:,2)];
U_in10 = pobs(X_in);

subplot(1,3,1);
scatter(U_in2(:,1), U_in2(:,2)); 
grid on; xlabel('u', 'FontSize', 20); ylabel('v', 'FontSize', 20); title('2-Outcomes Discrete RV', 'FontSize', 24);

subplot(1,3,2);
scatter(U(:,1), U(:,2)); 
grid on; xlabel('u', 'FontSize', 20); ylabel('v', 'FontSize', 20); title('True Underlying Copula', 'FontSize', 24);

subplot(1,3,3);
scatter(U_in10(:,1), U_in10(:,2)); 
grid on; xlabel('u', 'FontSize', 20); ylabel('v', 'FontSize', 20); title('10-Outcomes Discrete RV', 'FontSize', 24);
