function val = copfamprob(xx,rvEmpInfoObjs,copModel,copParams)


% compute the pseudo-observations of xx
uu = zeros(1,length(xx));
for ii=1:length(xx)
    uu(ii) = rvEmpInfoObjs{ii}.cdf(xx(ii));
end
% compute the product of the marginal densities
pi = 1;
for ii=1:length(xx)
    pi = pi * rvEmpInfoObjs{ii}.pdf(xx(ii));
end

if(strcmpi(copModel,'frank'))
    copval = frankcopulapdf(uu, copParams);
elseif(strcmpi(copModel,'gumbel'))
    copval = gumbelcopulapdf(uu, copParams);
elseif(strcmpi(copModel,'clayton'))
    copval = claytoncopulapdf(uu, copParams);
elseif(strcmpi(copModel,'gaussian'))
    copval = copulapdf('Gaussian', uu, copParams);
end
val = pi*copval*-1;     % we multiply by negative one so this turns
                        % into a maximization rather than minimization
                            
end