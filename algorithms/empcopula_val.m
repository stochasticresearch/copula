function [ C_u, c_u ] = empcopula_val( C, c, U )
%EMPCOPULA_VAL Returns the value of C (the copula function) and 
% c (the copula density) at the specified U (vector).
% Inputs:
%  C - the empirical copula function, defined over an evenly spaced grid of
%      points
%  c - the empirical copula density, defined over an evenly spaced grid of 
%      points
%  U - A vector [U_1, U_2, ... U_d], where d is the dimensionality
%      of the copula (distribution and density functions)
%
% Outputs:
%  C_u - the value of C(U), if C is empty, C_u will be 0
%  c_u - the value of c(U), if c is empty, c_u will be 0
%
% TODO:
%  [ ] - U should have the ability to be a matrix as well


d = length(size(C));
dprime = length(size(c));

if(length(U)~=d || (d~=dprime))
        error('Error in dimensionality matching of provided arguments!');
end
if(~isequal(size(C),size(c)))
    error('The grid over which C and c are defined need to be equal!');
end
if(range(size(C)) || range(size(c)))
    error('The grid over which C and c are defined should be a square!');
end

C_u = 0;
c_u = 0;

% calculate which grid points land in the U values of interest
linearIdx = 0;
for ii=1:d
    K_ii = size(C,ii);
    idx = round(U(ii)*K_ii);
    if(idx>K_ii)
        idx = K_ii;
    elseif(idx<1)
        idx = 1;
    end
    
    if(ii==1)
        linearIdx = idx;
    else
        linearIdx = linearIdx + (idx-1)*(K_ii.^(ii-1));
    end
end

% extract value of C_u and c_u from those
if(~isempty(C_u))
    C_u = C(linearIdx);
end
if(~isempty(c_u))
    c_u = c(linearIdx);
end

end

