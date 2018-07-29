% function [ecop] = empcopula(U)
% % Sample the empirical copula based on the given points
% 
% M = size(U,1);
% ecop = zeros(M,1);
% 
% for ii=1:M
% %     ecop(ii) = sum(u(ii)>=u(1:M) & v(ii)>=v(1:M))/(M+1);
%     ecop(ii) = sum( sum( U(ii,:) >= U(1:M,:), 2 ) );
% end
% 
% end

function [ecop] = empcopula(u,v)
n = length(u);
ecop = zeros(n,1);

parfor ii=1:n
    ecop(ii) = sum(u(ii)>=u & v(ii)>=v)/(n+1);
end

end