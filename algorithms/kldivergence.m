%**************************************************************************
%* 
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>
%*
%* This program is free software: you can redistribute it and/or modify
%* it under the terms of the GNU General Public License as published by
%* the Free Software Foundation, either version 3 of the License, or
%* (at your option) any later version.
%*
%* This program is distributed in the hope that it will be useful,
%* but WITHOUT ANY WARRANTY; without even the implied warranty of
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%* GNU General Public License for more details.
%*
%* You should have received a copy of the GNU General Public License
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.
%* 
%**************************************************************************

function [ kl_val ] = kldivergence( p, q, x)
%KLDIVERGENC Computes the KL-Divergence between P and Q
% Inputs:
%  p - the "reference" UNIVARIATE distribution
%  q - the "estimating" UNIVARIATE distribution
%  x - the domain over which p and q are defined (they should be the same
%                                                 for p and q)
%      if x is not provided, we approximate integration w/ a sum.  If x is
%      provided, then we use trapz to perform numerical integration

p = p(:);
q = q(:);

% tol = 1e-5;
% p(p<=tol)=tol;       % TODO: make this an argument?  provides log(0) protection
% q(q<=tol)=tol;       % TODO: make this an argument?  provides /0 protection

y = p.*log(p./q);
y(p==0) = 0;        % this is the usual convention when implementing KL-Div
if(nargin>2)
    kl_val = trapz(x,y);
else
    kl_val = sum(y);
end

end