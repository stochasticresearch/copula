function [ y ] = log1mexp( a )
%LOG1MEXP Computes log(1-exp(-a)), a convenience function as we use it a
%lot in copula PDF computations

y = log1p(-exp(-a));

end

