function [ pdf, combos ] = hist_discrete( X )
%An empirical histogram for N-dimensional discrete distributions
% TODO: currently only supports D=2 or D=3 and is very slow :(

M = size(X,1);
D = size(X,2);

uniqueVals = zeros(1,D);
for ii=1:D
    uniqueVals(ii) = length(unique(X(:,ii)));
end

if(D==1)
    combos = 1:uniqueVals(1);
elseif(D==2)
    combos = combvec(1:uniqueVals(1), 1:uniqueVals(2));
elseif(D==3)
    combos = combvec(1:uniqueVals(1), 1:uniqueVals(2), 1:uniqueVals(3));
else
    error('D>3 currently unsupported!');
end

pdf = zeros(1,size(combos,2));
for ii=1:size(combos,2)
    combo = combos(:,ii);
    % count the number of times this combo occured
    cnt = 0;
    for jj=1:M
        if(D==1)
            if(combo(1)==X(jj,1))
                cnt = cnt + 1;
            end
        elseif(D==2)
            if(combo(1)==X(jj,1) && combo(2)==X(jj,2))
                cnt = cnt + 1;
            end
        elseif(D==3)
            if(combo(1)==X(jj,1) && combo(2)==X(jj,2) && combo(3)==X(jj,3))
                cnt = cnt + 1;
            end
        end
    end
    pdf(ii) = cnt/M;
end

end

