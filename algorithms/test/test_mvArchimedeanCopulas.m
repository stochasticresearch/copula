% a test script for testing the multivariate (N>2) versions of generating
% random variates from the Frank, Gumbel, and Clayton copulas

M = 1000;
N = 3;

alpha = 5;

U_clayton = claytoncopularnd(M,N,alpha);
U_frank   = frankcopularnd(M,N,alpha);
U_gumbel  = gumbelcopularnd(M,N,alpha);

figure;
plotmatrix(U_clayton);
title('Multivariate Clayton Copula')

figure;
plotmatrix(U_frank);
title('Multivariate Frank Copula')

figure;
plotmatrix(U_gumbel);
title('Multivariate Gumbel Copula')
