function A = slice_assign(A, ix, dim, B)
%SLICE_ASSIGN Assign new values to a "slice" of A
% from: https://stackoverflow.com/questions/22537326

subses = repmat({':'}, [1 ndims(A)]);
subses{dim} = ix;
A(subses{:}) = B;