function [ Y ] = opWeightedSum( X, W )
%CMBWEIGHTEDSUM Weights X with W and sums the columns together
% Inputs:
%  x - matrix of data to combine, should be of dimension [m x n]
%  w - a matrix of weights, should be the same dimension as x
% Outputs:
%  y - the output

% TODO: input argument checking

Y = sum(X.*W,2);  % perform element wise weighting, and sum all columns

end

