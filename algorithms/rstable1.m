function [ y ] = rstable1(n, alpha, beta, gamma, delta, pm)
% Generates a value from the stable distribution.
% Algorithms copied directly from R source code of the copula package
%    - rstable1.R
%    - retstable.c

y = rstable_c(n, alpha) * gamma + delta;

end


function [ y ] = rstable_c(n, alpha)
    y = cos(pi/2.0*alpha).^(-1.0/alpha) * rstable0(alpha);
end

function [ y ] = rstable0(alpha)
    if(alpha == 1) 
        y = 1.0;
    else
        U = rand();
        while 1
            % generate non-zero exponential random variable
            W = exprnd(1);
            if(W~=0)
                break;
            end
        y = (A(pi*U,alpha)/(W.^(1.0-alpha)) ).^(1.0/alpha);
    end
end

function [ y ] = A(x, alpha)
    Ialpha = 1.0-alpha;
    y = A_3(x, alpha, Ialpha);
end

function [ y ] = A_3(x, alpha, Ialpha)
    y = ( (Ialpha* sinc(Ialpha*x/pi)).^Ialpha ) * ...
        ( ( (alpha * sinc(alpha *x/pi)).^alpha ) ./ sinc(x./pi) );
end