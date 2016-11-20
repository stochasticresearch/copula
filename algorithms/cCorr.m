function [metric] = cCorr( X, Y, bw, M )
% Computes the copula-correlation between realizations of random vectors X
% and Y.  The copula-correlation is loosely defined as the Lp distance
% between the copula of X and Y, and the independence copula
% Inputs:
%  X - realizations of X random vector
%  Y - realizations of Y random vector
%  bw - the bandwidth used for KDE estimation
%  M - the number of points used in numerical integration, defaults to 200
% Outputs:
%  metric - a value between 0 and 1 which is the copula correlation between
%           vectors X and Y
% See paper: Copula Correlation: An Equitable Dependence Measure and 
%            Extension of Pearson?s Correlation.
%            arXiv:1312.7214
% TODO:
%   [ ] - Extend to multivariate, currently only handles 1-D

if(~isequal(size(X),size(Y)))
    error('X and Y must be the same dimensions!');
end
if(size(X,2)>1)
    error('Currently, ccorr only supports 1-D vectors for X and Y!');
end

n = length(X);

if(nargin<3)
    bw = 0.25*n^(-1/4);
    M = 200;
elseif(nargin<4)
    M = 200;
end
maxc=ccorercpp([1:n]',[1:n]',bw,M);
minc=minfc_p(n,bw,M);
metric=(ccorercpp(X,Y,bw,M)-minc)/(maxc-minc);

if(metric>1)
    metric = 1;
elseif(metric<0)
    metric = 0;
end

end

function [s] = minfc_p(n, bw, M)

k=ceil(2*bw*(n+1));
x=floor(n/k);
xin=1:n;
yin = [];
% TODO: figure out a more efficient way to creating yin vector
for ii=1:k
    yin = [yin ii+k*(0:x)];
end
y_select = yin(yin<=n);
s = ccorercpp(xin',y_select',bw,M);

if(length(xin)~=length(y_select))
    error('lengths not equal?');
end

end

% TODO: make this a C offload function
function [s] = ccorercpp(x,y,bw,M)
n = length(x);
u = pobs(x);
v = pobs(y);
h = bw;

A = zeros(M,M);

s=0.0;

for kk=1:n
    ul=max(floor((u(kk)-h)*M)+1,1);
    uu=min(floor((u(kk)+h)*M)+1,M);
    vl=max(floor((v(kk)-h)*M)+1,1);
    vu=min(floor((v(kk)+h)*M)+1,M);
    
    for ir=ul:uu
        for ic=vl:vu
            A(ir,ic) = A(ir,ic) + 1;
        end
    end
end

for ii=1:M
    for jj=1:M
        A(ii,jj) = max(A(ii,jj)/(n*h*h*4)-1,0.0);
        s = s + A(ii,jj);
    end
end

s = s/(M^2);
end