function [ y ] = logserrnd( p, n )
%LOGSERRND generates samples from the Log-Series discrete distribution with
%parameter p
% Inputs:
%  p - shape parameter of log-series distribution
%  n - the number of random variates to generate
% Outputs:
%  y - the random variates

% We approximate the log-series distribution by converting it to a multinomial
% distribution with K = 20
K = 20;
probvec = zeros(1,K);
for kk=1:K
    probvec(kk) = - p^kk ./ (kk*log(1-p));
end

pd = makedist('Multinomial','Probabilities',probvec);

% generate random variates from this created probability distribution
y = random(pd, 1, n);

end