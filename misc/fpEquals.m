function [res] = fpEquals(a, b, tol)
%FPEQUALS - checks if two floating point numbers are equal to a optionally
%specified tolerance

if(nargin<3)
    tol = 10*eps;
end

if(abs(a-b)<tol)
    res = 1;
else
    res = 0;
end

end