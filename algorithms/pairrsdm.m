function [R] = pairrsdm( X )
%PAIRRSDM - computes pairwise dependency metrics of a given vector
% Inputs:
%  X - a matrix of observations from which pairwise RSDM metrics are
%  computed.
%  varargin{1} - minscanincr - the minimum scanning increment.  Large
%                              values will filter out high frequency
%                              dependencies, small values decrease the
%                              statistical power of the dependency metric
%  varargin{2} - diffthresh  - the threshold at which a change in
%                              concordance amount is detected.  Larger
%                              values are more robust to noise, but tend to
%                              miss high frequency changes.
%  varargin{3} - alpha       - the value used to determine significance
%                              level of a box's concordance level
% Outputs:
%  R - a matrix which contains pairwise correlations of all columns in X
%**************************************************************************
%*                                                                        *
%* Copyright (C) 2016  Kiran Karra <kiran.karra@gmail.com>                *
%*                                                                        *
%* This program is free software: you can redistribute it and/or modify   *
%* it under the terms of the GNU General Public License as published by   *
%* the Free Software Foundation, either version 3 of the License, or      *
%* (at your option) any later version.                                    *
%*                                                                        *
%* This program is distributed in the hope that it will be useful,        *
%* but WITHOUT ANY WARRANTY; without even the implied warranty of         *
%* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
%* GNU General Public License for more details.                           *
%*                                                                        *
%* You should have received a copy of the GNU General Public License      *
%* along with this program.  If not, see <http://www.gnu.org/licenses/>.  *
%*                                                                        *
%**************************************************************************

% % minscanincr_in = 0;
% % diffthresh_in = 0;
% % alpha_in = 0;
% % 
% % % overwrite defaults with user-inputted values
% % nVarargin = length(varargin);
% % switch nVarargin
% %     case 0
% %     case 1
% %         minscanincr = varargin{1};
% %         
% %         minscanincr_in = 1;
% %     case 2
% %         minscanincr = varargin{1};
% %         diffthresh = varargin{2};
% %         
% %         minscanincr_in = 1;
% %         diffthresh_in = 1;
% %     otherwise
% %         % the >=3 case
% %         minscanincr = varargin{1};
% %         diffthresh = varargin{2};
% %         alpha = varargin{3};
% %         
% %         minscanincr_in = 1;
% %         diffthresh_in = 1;
% %         alpha_in = 1;
% % end

n = size(X,2);      % the dimensionality of the dataset
R = zeros(n,n);

for ii=1:n
    x = X(:,ii);
    parfor jj=ii+1:n
        y = X(:,jj);
        
% % %         if(minscanincr_in && diffthresh_in && alpha_in)
% % %             rsdmVal = rsdm(x,y,minscanincr,diffthresh,alpha);
% % %         elseif(minscanincr_in && diffthresh_in)
% % %             rsdmVal = rsdm(x,y,minscanincr,diffthresh);
% % %         elseif(minscanincr_in)
% % %             rsdmVal = rsdm(x,y,minscanincr);
% % %         else
% % %             rsdmVal = rsdm(x,y);
% % %         end
        rsdmVal = rsdm(x,y);
        R(ii,jj) = rsdmVal;
    end
end

% assign the lower triangle by copying the upper-triangle of the R matrix
% make matrix symmetric.  The below works because R is initialized to
% zeros.
R=R+R';
R(1:n+1:n*n) = 1;   % set the diagonal to 1

end