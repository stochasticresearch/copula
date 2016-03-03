function [ y ] = opAffine( x, A, b )
%DEPAFFINE Induces affine dependency on X w/ aX+b
% Inputs:
%  x - the variable to induce affine dependency upon - should be [m x 1]
%  A - the slope - should be [m x m]
%  b - the intercept - should be [m x 1]
% Outputs:
%  y - the output

% TODO: input argument checking

y = A*x + b;

end

