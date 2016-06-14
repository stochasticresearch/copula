function [c,ceq] = nlcon(theta)

D = length(theta);
c = [];
ceq = [];
ceqidx = 1;
for ii=1:D
    for jj=1:D
        if(jj>ii)
            ceq(ceqidx) = theta(ii,jj)*theta(jj,ii);
            ceqidx = ceqidx + 1;
        end
    end
end