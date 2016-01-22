function [ U ] = gumbelcopularnd( M, N, alpha )
%GUMBELCOPULARND Generates M samples from a Gumbel copula of dimensionality
%N, with parameter alpha
% Inputs:
%  M - the number of samples to generate
%  N - the dimensionality of the data
%  alpha - the dependency parameter of the Gumbel copula
%
% Outputs:
%  U - an M x N matrix of generated samples

if(N==2)
    U = copularnd('Gumbel', alpha, N);
else
    % Algorithm 1 described in both the SAS Copula Procedure, as well as the
    % paper: "High Dimensional Archimedean Copula Generation Algorithm"
    U = zeros(M,N);
    
    for ii=1:M
        a  = 1.0/alpha;
        b  = 1;
        g  = cos(pi/(2.0*alpha)).^alpha;
        d  = 0;
        pm = 1;
        v = rstable1(1,a,b,g,d,pm);
        
        % sample N independent uniform random variables
        x_i = rand(1,N);
        t = -1*log(x_i)./v;
        
        U(ii,:) = exp(-1*(t.^1.0/alpha));
    end
    
end % if

end % function