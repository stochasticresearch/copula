function [ c ] = empcopdens_betak_3d_v2( X1,X2,X3,h,K )

K = K + 1;
c = zeros(K,K,K);

for ii=1:K
    for jj=1:K
        for kk=1:K
            uu = ii/K;
            vv = jj/K;
            ww = kk/K;
            
            K1 = betapdf(X1, uu/h + 1, (1-uu)/h + 1);
            K2 = betapdf(X2, vv/h + 1, (1-vv)/h + 1);
            K3 = betapdf(X3, ww/h + 1, (1-ww)/h + 1);

            val = sum(K1.*K2.*K3)/length(X1);
            c(ii,jj,kk) = val;
        end
    end
end

c = c(1:K-1,1:K-1,1:K-1);

end

