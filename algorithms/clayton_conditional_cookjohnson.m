function U = clayton_conditional_cookjohnson(u1,alpha)
n = length(u1);
p = rand(n,1);
if alpha < sqrt(eps)
    u2 = p;
else
    u2 = u1.*(p.^(-alpha./(1+alpha)) - 1 + u1.^alpha).^(-1./alpha);
end
U = [u1 u2];
end