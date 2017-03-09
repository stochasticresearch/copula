function [ s ] = stirling2(m, n)
%STIRLING2 Computes stirling number of the seconds kind.
%   Caches numbers up to 32 dimensions.
% Acknowledgements:
%  https://github.com/mscavnicky/copula-matlab/blob/master/%2Barchim/generatorDerivative.m

persistent cache;
if isempty(cache)
   cache = NaN(32, 32);
   % Put ones on diagonal
   cache(eye(32) == 1) = 1;
   % Zeros in the first column
   cache(2:32,1) = 0;
   % Ones in the second column
   cache(2:32,2) = 1;   
end

s = cache(m+1,n+1);
if isnan(s)
   s = stirling2(m-1, n-1) + n*stirling2(m-1, n);
   cache(m+1,n+1) = s;
end

end