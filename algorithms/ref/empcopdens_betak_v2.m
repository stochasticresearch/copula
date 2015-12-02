function [ c, a_vec, b_vec ] = empcopdens_betak_v2( u,v,h,K )

K = K + 1;
c = zeros(K,K);
a_vec = zeros(K,K);
b_vec = zeros(K,K);
for ii=1:K
    for jj=1:K
        a = ii/K;
        b = jj/K;
        
        K1 = betapdf(u,a/h + 1,(1-a)/h + 1);
        K2 = betapdf(v,b/h + 1,(1-b)/h + 1);
        
        c(ii,jj) = sum(K1.*K2)/length(u);
        a_vec(ii,jj) = a;
        b_vec(ii,jj) = b;
    end
end

c = c(1:K-1,1:K-1);

end

