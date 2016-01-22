function [ U ] = frankcopularnd( M, N, alpha )
%FRANKCOPULARND Generates M samples from a Frank copula of dimensionality
%N, with parameter alpha
% Inputs:
%  M - the number of samples to generate
%  N - the dimensionality of the data
%  alpha - the dependency parameter of the Gumbel copula
%
% Outputs:
%  U - an M x N matrix of generated samples

if(N==2)
    U = copularnd('Frank', alpha, N);
else
    % Algorithm 1 described in both the SAS Copula Procedure, as well as the
    % paper: "High Dimensional Archimedean Copula Generation Algorithm"
    if(alpha<=0)
        error('For N>=3, alpha > 0 for the Frank Copula')
    end
    
    U = zeros(M,N);
    for ii=1:M
        p = -1.0*expm1(-1*alpha);
        if(p==1)
            % boundary case protection
            p = 1 - eps;
        end
        v = logserrnd(p, 1);
        
        % sample N independent uniform random variables
        x_i = rand(1,N);
        t = -1*log(x_i)./v;
        U(ii,:) = -1.0*log1p( exp(-t)*expm1(-1.0*alpha))/alpha;
    end
    
end % if

end % function