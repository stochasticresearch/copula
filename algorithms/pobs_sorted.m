function [u,v] = pobs_sorted(x,y,scaleFlag)

% rescale back to 0-1
u = pobs(x);
v = pobs(y);

if(scaleFlag)
    M = length(x);
    u = u*(M+1)/M;
    v = v*(M+1)/M;
end

[u,I] = sort(u);
v = v(I);

end