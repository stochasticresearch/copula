%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [ y ] = logserrnd( p, n, varargin )
%LOGSERRND generates samples from the Log-Series discrete distribution with
%parameter p - approximated as a multinomial distribution
% Inputs:
%  p - shape parameter of log-series distribution
%  n - the number of random variates to generate
% Optional Inputs:
%  K - the number of discrete outcomes to approximate to.  The larger the
%      better, defaults to 20
% Outputs:
%  y - the random variates

% We approximate the log-series distribution by converting it to a multinomial
% distribution with K = 1000 (default)
K = 1000;

nVarargs = length(varargin);
if(nVarargs==1)
    K = varargin{1};
    % validate it is an integer below 2000
    if(~isnumeric(K) || K>2000 || K<0)
        warning('Invalid varargin{1}, defaulting to K=1000');
        K = 1000;
    end
end

probvec = zeros(1,K);
for kk=1:K
    probvec(kk) = - p^kk ./ (kk*log1p(-p));
end

% normalize probvec to sum to 1
probvec = probvec/sum(probvec);

pd = makedist('Multinomial','Probabilities',probvec);

% generate random variates from this created probability distribution
y = random(pd, 1, n);

end