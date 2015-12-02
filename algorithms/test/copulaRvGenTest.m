% Experiment with generating RV's from empirically generated copula's
% (First we test w/ known forms, but using copula dnesity functions instead
% of the built-in copularnd functions)

function [] = copulaRvGenTest()

    D = 2;
    M = 1000;
    u_gen = rand(M,D);

    % precompute copula density function
    K = 1000;
    uu = linspace(0,1,K);
    [U1,U2] = meshgrid(uu,uu);
    alpha = 10;
    copType = 'Gumbel';
    c2 = copulapdf(copType,[U1(:) U2(:)],alpha);
    c2 = reshape(c2,K,K);

    c1 = squeeze(sum(c2,2));

    for ii=1:M
        u = u_gen(ii,1);
        t = u_gen(ii,2);

        idx = findClosest(uu,u);
        
        C_d_num = cumsum(c2(idx,:));
        C_d_den = c1(idx);
        C_d = C_d_num./C_d_den;

        % perform numerical inverse
        v = uu(findClosest(C_d,t));

        % assign random variate
        u_gen(ii,2) = v;
    end

    subplot(1,2,1);
    scatter(u_gen(:,1),u_gen(:,2)); title('Generated')

    subplot(1,2,2);
    uu = copularnd(copType,alpha,M);
    scatter(uu(:,1),uu(:,2)); title('Built-in')

end

% function [idx] = genIdx(K,u)
%     idx = round(u*K);
%     if(idx>K)
%         idx=K;
%     elseif(idx<2)
%         idx=2;
%     end
% 
% end

function [idx] = findClosest(vec, val)
    tmp = abs(val-vec);
    [~,idx] = min(tmp);
end