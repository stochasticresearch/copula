M = 2000;
x = rand(M,1)*3-1.5;
% x = rand(M,1);

y1 = x;
y2 = 4*(x-0).^2;
y3 = 128*(x-1/3).^3-48*(x-1/3).^3-12*(x-1/3);
y4 = sin(4*pi*x);
y5 = sin(4*pi*x).*exp(-x);
y6 = sin(2*pi*(x + x.^2));
% y7=(2*binornd(1,0.5,M,1)-1) .* (sqrt(1 - (2*x - 1).^2));
% y8 = (x > 0.5);
y7 = 4*(x-.3).^2;
y8 = (x > 0.5);

figure;
subplot(4,4,1); scatter(x,y1); grid on; title('Linear'); xlabel('x'); ylabel('y');
subplot(4,4,2); scatter(x,y2); grid on; title('Parabolic'); xlabel('x'); ylabel('y');
subplot(4,4,3); scatter(x,y3); grid on; title('Cubic'); xlabel('x'); ylabel('y');
subplot(4,4,4); scatter(x,y4); grid on; title('Sine-1'); xlabel('x'); ylabel('y');
subplot(4,4,5); scatter(pobs(x)*(M+1)/M,pobs(y1)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');
subplot(4,4,6); scatter(pobs(x)*(M+1)/M,pobs(y2)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');
subplot(4,4,7); scatter(pobs(x)*(M+1)/M,pobs(y3)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');
subplot(4,4,8); scatter(pobs(x)*(M+1)/M,pobs(y4)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');
subplot(4,4,9); scatter(x,y5); grid on; title('Damped Sinusoid'); xlabel('x'); ylabel('y');
subplot(4,4,13); scatter(pobs(x)*(M+1)/M,pobs(y5)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');
subplot(4,4,10); scatter(x,y6); grid on; title('Chirp'); xlabel('x'); ylabel('y');
subplot(4,4,14); scatter(pobs(x)*(M+1)/M,pobs(y6)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');
subplot(4,4,11); scatter(x,y7); grid on; title('Parabolic Skewed'); xlabel('x'); ylabel('y');
subplot(4,4,15); scatter(pobs(x)*(M+1)/M,pobs(y7)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');
subplot(4,4,12); scatter(x,y8); grid on; title('Step'); xlabel('x'); ylabel('y');
subplot(4,4,16); scatter(pobs(x)*(M+1)/M,pobs(y8)*(M+1)/M, 'r'); grid on; xlabel('u'); ylabel('v');


%% 
figure;
x = rand(M,1)*10;
y1 = exp(x);
y2 = y1 + randn(M,1)*0.5;
y3 = 0.25*x;
y4 = y3 + randn(M,1)*0.5;

tau1 = corr(x,y1,'type','kendall');
tau2 = corr(x,y2,'type','kendall');
tau3 = corr(x,y3,'type','kendall');
tau4 = corr(x,y4,'type','kendall');

srho1 = corr(x,y1,'type','spearman');
srho2 = corr(x,y2,'type','spearman');
srho3 = corr(x,y3,'type','spearman');
srho4 = corr(x,y4,'type','spearman');


subplot(2,4,1); scatter(x,y1); grid on; 
title(sprintf('Exponential - \\tau=%0.02f \\rho_s=%0.02f',tau1,srho1)); xlabel('x'); ylabel('y');

subplot(2,4,2); scatter(x,y2); grid on; 
title(sprintf('Noisy Exponential - \\tau=%0.02f \\rho_s=%0.02f',tau1,srho1)); xlabel('x'); ylabel('y');

subplot(2,4,3); scatter(x,y3); grid on; 
title(sprintf('Linear - \\tau=%0.02f \\rho_s=%0.02f',tau1,srho1)); xlabel('x'); ylabel('y');

subplot(2,4,4); scatter(x,y4); grid on; 
title(sprintf('Noisy Linear - \\tau=%0.02f \\rho_s=%0.02f',tau1,srho1)); xlabel('x'); ylabel('y');

subplot(2,4,5); scatter(pobs(x)*(M+1)/M,pobs(y1)*(M+1)/M); grid on; xlabel('u'); ylabel('v');
subplot(2,4,6); scatter(pobs(x)*(M+1)/M,pobs(y2)*(M+1)/M); grid on; xlabel('u'); ylabel('v');
subplot(2,4,7); scatter(pobs(x)*(M+1)/M,pobs(y3)*(M+1)/M); grid on; xlabel('u'); ylabel('v');
subplot(2,4,8); scatter(pobs(x)*(M+1)/M,pobs(y4)*(M+1)/M); grid on; xlabel('u'); ylabel('v');
