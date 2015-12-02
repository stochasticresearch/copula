function [ c ] = empcopdens_betak( u,v,h,K )

K = K + 1;
c = zeros(K,K);
for ii=1:K
    for jj=1:K
        a = ii/K;
        b = jj/K;
        c(ii,jj) = sum(betapdf(a,u/h,(1-u)/h).*betapdf(b,v/h,(1-v)/h))/length(u);
    end
end
c = c(1:K-1,1:K-1);

end

