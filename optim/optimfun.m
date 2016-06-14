function y = optimfun(theta, data)

% WARNING - the way Rho_num and Rho_den are assigned, it makes the
% assumption that the parent variables are NOT interacting with each other
% I do this as a simplification for now, but also, another problem is, how
% do we ensure that the Rho matrix is positive-definite? how can that
% be a constraint, I don't think we can put it as a constraint ...

D = length(theta);
y = 0;
for dd=1:D
    % find which nodes this ndoe is connected to
    parentNodes = find(theta(:,dd)~=0); parentNodes = parentNodes';
    % compute the Rho matrix for the gaussian copula w/ this
    if(~isempty(parentNodes))
        Rho_num = eye(length(parentNodes)+1);
        for ii=1:length(parentNodes)
            Rho_num(1,ii+1) = theta(ii,dd);
            Rho_num(ii+1,1) = theta(ii,dd);
        end
        u_num = [data(:,dd) data(:,parentNodes)];
        Rc_num = copulapdf('Gaussian', u_num, Rho_num);
        if(length(parentNodes)>1)
            Rho_den = eye(length(parentNodes));
            u_den = data(:,parentNodes);
            Rc_den = copulapdf('Gaussian', u_den, Rho_den);
        else
            Rc_den = 1;
        end
        tmp = Rc_num./Rc_den;
        y = y + sum(tmp);
    end
end

% we are minimizing the negative log-likelihood, same as maximizing
% log-likelihood
if(y==0)
    y = 10000;
else
    y = -1*log(y);
end

end