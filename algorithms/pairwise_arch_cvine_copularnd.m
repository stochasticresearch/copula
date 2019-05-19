function U = pairwise_arch_cvine_copularnd(copula_type,corrvec,cases)
% depends on the MixedVineToolbox for Matlab
% https://stommac.eu/index.php/code
% also copied into stochasticresearch github repo here:
%  https://github.com/stochasticresearch/MixedVineToolbox_v1.0

% setup the vine
d = length(corrvec)+1;
vine.type = 'c-vine';
vine.families = cell(d);
vine.theta = cell(d);

for ii=2:d
    vine.families{1,ii} = copula_type;
    theta = copulaparam(copula_type, corrvec(ii-1), 'type', 'spearman');
    vine.theta{1,ii} = theta;
end
% specify the remainder of the upper triangle
for ii=2:d
    for jj=ii+1:d
        vine.families{ii,jj} = 'ind';
        vine.theta{ii,jj} = [];
    end
end

% for ii=1:d
%     for jj=1:d % we can do this due to symmetricity of copula
%         R_val = R(ii,jj);
%         if ii==jj
%             ;
%         elseif (R_val~=0)
%             vine.families{ii,jj} = copula_type;
%             % convert the R value (correlation coefficient) to clayton
%             % copula dep
%             theta = copulaparam(copula_type, R(ii,jj), 'type', 'spearman');
%             fprintf('R(ii,jj)=%0.02f, alpha=%0.02f\n', R(ii,jj),theta);
%             vine.theta{ii,jj} = theta;
%         else
%             vine.families{ii,jj} = 'ind';
%             vine.theta{ii,jj} = [];
%         end
%     end
% end

vine.families
vine.theta

% sample from the vine (this is copied from mixedvinernd.m)
w = rand(cases,d);
v = zeros(cases,d,d);
v(:,1,1) = reshape(w(:,1),[cases 1 1]);
for i = 2:d
    v(:,i,1) = reshape(w(:,i),[cases 1 1]);
    for k = (i-1):-1:1
        v(:,i,1) = copulaccdfinv(vine.families{k,i},[v(:,k,k) v(:,i,1)],vine.theta{k,i},1);
    end
    if i < d
        for j = 1:(i-1)
            v(:,i,j+1) = copulaccdf(vine.families{j,i},[v(:,j,j) v(:,i,j)],vine.theta{j,i},1);
        end
    end
end
U = reshape(v(:,:,1),size(w));

end