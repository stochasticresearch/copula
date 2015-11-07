% A script which contains some test simulations for the Paper "Discrete
% Copula's," published in IEEE Transactions on Fuzzy Systems 2006

function [C] = irreducibleCopula(n)

P = perms(1:n);
C = zeros(size(P,1), n*n);

for zz=1:size(P,1)
    fprintf('Calculating permutation = %s\n', sprintf('%d ', P(zz,:)))
    
    ti = 1;
    for ii=1:n
        for jj=1:n
            C(zz, ti) = func(P(zz,:), ii, jj);
            ti = ti + 1;
        end
    end
end

end

function [y] = func(perm, i, j)

n = length(perm);

y = 0;
for k=1:j
    indic = 1:i;
    if( any(perm(k)==indic) )
        y = y + 1;
    end
end

y = y / n;
end