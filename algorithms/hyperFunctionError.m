function [ e ] = hyperFunctionError( f1, f2 )
%COMPUTEERRORBETWEENSURFACES Computes the error between two hyperfunctions
% as follows: sqrt( integral( (f1-f2).^2 ) )
%where the integral is computed over the entire domain of the functions.
%The two surfaces are expected to have the same domain.  I envision the 
% usage of this function to be a sort of a test of how good of an estimator 
% you have, where surf1 can be a surface (of any dimension) of the actual, 
% and surf2 can be an estimate.
% Inputs:
%  f1 - hyperfunction (of any dimension) one
%  f2 - hyperfunction (of the same dimension as f1) two
% Outputs:
%  e - the error

e = sqrt( nansum( (f1(:)-f2(:)).^2 ) ) / length(f1(:));    % normalize

end