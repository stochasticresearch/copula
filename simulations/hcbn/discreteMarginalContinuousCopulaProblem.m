clear;
clc;

u = 0:0.1:1; v = 1;

% Check the Gaussian Copula
rhoVec = -0.9:0.1:0.9;
% Check v
for rho=rhoVec
    for uu=1:length(u)
        val = copulapdf('Gaussian', [u(uu) v], rho);
        if(val~=0)
            fprintf('val=%f -- rho=%f u=%f \n', val, rho, u);
        end
    end
end

% Check u
for rho=rhoVec
    for uu=1:length(u)
        val = copulapdf('Gaussian', [v u(uu)], rho);
        if(val~=0)
            fprintf('val=%f -- rho=%f u=%f \n', val, rho, u);
        end
    end
end

% Check Frank Copula
alphaVec = 1:5;
% Check v
for alpha=alphaVec
    for uu=1:length(u)
        val = copulapdf('Frank', [u(uu) val], alpha);
        if(val~=0)
            fprintf('val=%f -- alpha=%f u=%f \n', val, alpha, u);
        end
    end
end

% Check v
for alpha=alphaVec
    for uu=1:length(u)
        val = copulapdf('Frank', [val u(uu)], alpha);
        if(val~=0)
            fprintf('v=%f -- alpha=%f u=%f \n', val, alpha, u);
        end
    end
end

% Check Gumbel Copula
alphaVec = 1:5;
% Check v
for alpha=alphaVec
    for uu=1:length(u)
        val = copulapdf('Gumbel', [u(uu) val], alpha);
        if(val~=0)
            fprintf('val=%f -- alpha=%f u=%f \n', val, alpha, u);
        end
    end
end

% Check v
for alpha=alphaVec
    for uu=1:length(u)
        val = copulapdf('Gumbel', [val u(uu)], alpha);
        if(val~=0)
            fprintf('v=%f -- alpha=%f u=%f \n', val, alpha, u);
        end
    end
end

% Check Clayton Copula
alphaVec = 1:5;
% Check v
for alpha=alphaVec
    for uu=1:length(u)
        val = copulapdf('Clayton', [u(uu) val], alpha);
        if(val~=0)
            fprintf('val=%f -- alpha=%f u=%f \n', val, alpha, u);
        end
    end
end

% Check v
for alpha=alphaVec
    for uu=1:length(u)
        val = copulapdf('Clayton', [val u(uu)], alpha);
        if(val~=0)
            fprintf('v=%f -- alpha=%f u=%f \n', val, alpha, u);
        end
    end
end