function [u,v] = pobs_sorted_cc(x,y)

M = length(x);
data_sorted = sort(x);
[~, u] = ismember(x,data_sorted);
data_sorted = sort(y);
[~, v] = ismember(y,data_sorted);
[u,I] = sort(u); 
v = v(I);
u = u/M; v = v/M;

end