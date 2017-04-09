function out = slice(A, ix, dim)
% from: https://stackoverflow.com/questions/22537326

subses = repmat({':'}, [1 ndims(A)]);
subses{dim} = ix;
out = A(subses{:});

end