function [Z_sorted, nU, F_emp] = empVF_v3(n,D,Z)

% D Dimensionen, in any empirical VF

%Empirical distribution:

step=[1/n:1/n:1]';  % [1/n, 2/n, ... ,1]

for d=1:D
   [Z_sorted(d,:), index_alt(d,:)] = sort(Z(d,:));
end 


U=zeros(D,n);     % to speed up processing

%Find  empirical CDF  F_emp
F_emp = zeros(4,n);
for d=1:D
    F_emp(d,:) = step';
end

  % Sampling values repeatedly have different values F i/n, (i+1)/n, ...
for i=1:n
    for d=1:D        
        nU(d,index_alt(d,i)) =  i; 
        % values n*U(d,index_alt) in order to omit rounding errors
    end                
end

end