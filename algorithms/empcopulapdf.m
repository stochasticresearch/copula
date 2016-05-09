function [ c ] = empcopulapdf( U, h, K, method )
%EMPCOPULAPDF Calculates the empirical copula density over the grid
%points with K points in each dimension, using the specified method
% 
% Inputs
%  U - the pseudosamples of the copula (i.e.) F(X) = [F(X_1) ... F(X_D)],
%      should be a M x D input, where M is the number of samples and D is
%      the dimensionality
%  h - the kernel bandwidth
%  K - the spacing of the gridpoints over which to calculate the empirical
%      copula density
%  method - options are:
%           betak - Use Beta Kernels to estimate copula density, as
%                   specified by Copula Density Estimation Book - by Arthur
%                   Charpentier
%
% Outputs
%  c - the copula density
%
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
D = size(U,2);
sz = ones(1,D)*K;
ndgridInput = cell(1,D);
for dd=1:D
    ndgridInput{dd} = linspace(0,1,K);
end
ndgridOutput = cell(1,numel(ndgridInput));
[ndgridOutput{:}] = ndgrid(ndgridInput{:});

% the loop below to generate the gridPoints emulates the behavior of teh
% following: [U1,U2,U3] = ndgrid(u); copulapdf('..', [U1(:) U2(:) U3(:)], ..]
gridPoints = zeros(numel(ndgridOutput{1}), D);
for dd=1:D
    gridPoints(:,dd) = reshape(ndgridOutput{dd},numel(ndgridOutput{dd}),1);
end

% gridPoints = reshape(cell2mat(ndgridOutput),K^D,D);  %% WARNING
if(strcmpi(method, 'betak'))
    % Beta-Kernel method w/ C implementation for speed improvement
    c = empcopulapdf_c(U, K, h, gridPoints);
    c = reshape(c,sz);
elseif(strcmpi(method, 'betak-matlab'))
    M = size(U,1);
    c = zeros(1,K^D);       % we create it as a vector, then reshape at the
                            % end to the proper dimensionality
    for ii=1:K^D
        % make the beta pdf param's vector
        Kernel_vec = ones(M,1);
        % NOTE: According to different Matlab threads, optimizing the loop
        % below to use cellfun would yield similar or worse performance,
        % because apparently cellfun is slower than a for loop when not
        % using it with built-in functions?
        % See:
        %  https://www.mathworks.com/matlabcentral/newsreader/view_thread/301894
        %  https://www.mathworks.com/matlabcentral/newsreader/view_thread/253815
        %  https://www.mathworks.com/matlabcentral/newsreader/view_thread/253596
        %  https://www.mathworks.com/matlabcentral/newsreader/view_thread/251700
        %  http://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
        for jj=1:D 
            gridPoint = gridPoints(ii,jj);
            Kernel_vec = Kernel_vec .* betapdf(U(:,jj), gridPoint/h + 1, (1-gridPoint)/h + 1);            
        end
        c(ii) = sum(Kernel_vec)/M;
    end                            
    c = reshape(c,sz);
elseif(strcmpi(method, 'betak-matlab-old'))
    M = size(U,1);
    c = zeros(1,K^D);       % we create it as a vector, then reshape at the
                            % end to the proper dimensionality
    
    for ii=0:K^D-1
        % make the beta pdf param's vector
        gridPoints = zeros(1,D);
        
        value = ii;
        xIdx = D;
        while(value > 0)
            d = mod(value,K);

            gridPoints(xIdx) = d;
            xIdx = xIdx - 1;

            value = floor(value/K);
        end
        gridPoints = gridPoints/(K-1);      % this calculates 0 ... 1 (spacing will be 1/K)
                                            % to maintain K points
        gridPoints = fliplr(gridPoints);    % flip this to calculate column wise so that reshape
                                            % will restore c to expected dimensions
        Kernel_vec = ones(M,1);
        for jj=1:D 
            gridPoint = gridPoints(jj);
            Kernel_vec = Kernel_vec .* betapdf(U(:,jj), gridPoint/h + 1, (1-gridPoint)/h + 1);            
        end

        c(ii+1) = sum(Kernel_vec)/M;
    end
                            
    c = reshape(c,sz);
else
    error('Only beta kernel method supported currently!');
end

end

