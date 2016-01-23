function [ y ] = logserrnd( p, n, varargin )
%LOGSERRND generates samples from the Log-Series discrete distribution with
%parameter p
% Inputs:
%  p - shape parameter of log-series distribution
%  n - the number of random variates to generate
% Optional Inputs:
%  K - the number of discrete outcomes to approximate to.  The larger the
%      better, defaults to 20
% Outputs:
%  y - the random variates

% We approximate the log-series distribution by converting it to a multinomial
% distribution with K = 20 (default)
K = 20;

nVarargs = length(varargin);
if(nVarargs==1)
    K = varargin{1};
    % validate it is an integer below 50
    if(~isinteger(K) || K>20 || K<0)
        warning('Invalid varargin{1}, defaulting to K=20');
        K = 20;
    end
end

probvec = zeros(1,K);
for kk=1:K
    probvec(kk) = - p^kk ./ (kk*log(1-p));
end

pd = makedist('Multinomial','Probabilities',probvec);

% generate random variates from this created probability distribution
y = random(pd, 1, n);

end