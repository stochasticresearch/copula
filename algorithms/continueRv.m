function [ X_out ] = continueRv( X_in )
%CONTINUERV Continues a random variable according to Michel and Denuit
%(2005), and Neslehova (2007).
% Inputs:
%  X_in - Input matrix of dimension M x D, where D is the dimensionality of
%         the data, and M is the number of samples.
%
% Outputs:
%  X_out - the continued samples of the random variables

X_out = X_in + (rand(size(X_in))-1);

end

