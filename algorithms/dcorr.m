function c = dcorr( X, Y )
%DCORR computes the distance correlation between two random variables
% X and Y. Rows represent the examples, and columns the variables.
% Date: 10-03-2015
% Author: Paolo Inglese (paolo.ingls@gmail.com)
% Reference: http://en.wikipedia.org/wiki/Distance_correlation 

a = pdist2(X, X); % matrix euclidean distances
b = pdist2(Y, Y);

A = a - bsxfun(@plus,mean(a),mean(a,2))+mean(a(:));
B = b - bsxfun(@plus,mean(b),mean(b,2))+mean(b(:));

dcovXY = sum(sum(A.*B)) / (size(X,1)*size(X,1));
dcovXX = sum(sum(A.*A)) / (size(X,1)*size(X,1));
dcovYY = sum(sum(B.*B)) / (size(X,1)*size(X,1));

c = sqrt(dcovXY / sqrt(dcovXX * dcovYY));

end