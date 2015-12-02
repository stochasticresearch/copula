clear;
clc;

D = 2;
K = 10;
x = zeros(K^D,D);

for ii=0:K^D-1
    
    value = ii;
    xIdx = D;
    while(value > 0)
        d = mod(value,K);
        
        x(ii+1,xIdx) = d;
        xIdx = xIdx - 1;
        
        value = floor(value/K);
    end
%     x(ii+1,:) = x(ii+1,:) + 1;
    x(ii+1,:) = x(ii+1,:)/(K-1);
end
