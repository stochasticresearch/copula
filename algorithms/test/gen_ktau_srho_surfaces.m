% generate some plots on what Kendall's Tau and Spearman's Rho's Z-score
% surfaces look like

% Kendall's Tau
tau = 0:.05:1;
M = 50:50:500;
z_tau = zeros(length(M),length(tau));
for mm=1:length(M)
    for tt=1:length(tau)
        z_tau(mm,tt) = tau(tt)./sqrt( (2*(2*M(mm)+5))./(9*M(mm).*(M(mm))) );
    end
end
subplot(1,2,1); surf(tau,M,z_tau); xlabel('\tau'); ylabel('M'); title('\tau Z-Score Surface');

% Spearman's Rho
srho = 0:.05:1;
M = 50:50:500;
z_srho = zeros(length(M),length(srho));
for mm=1:length(M)
    for tt=1:length(srho)
        rho_val = srho(tt);
        Fr = 0.5*log( (1+rho_val)./(1-rho_val) );
        z_srho(mm,tt) = sqrt( (M(mm)-3)/1.06 ).*Fr;
    end
end
subplot(1,2,1); surf(tau,M,z_tau); xlabel('\tau'); ylabel('M'); title('\tau Z-Score Surface');
subplot(1,2,2); surf(srho,M,z_srho); xlabel('\rho_s'); ylabel('M'); title('\rho_s Z-Score Surface');
