function [optimalModel, modelParams] = copulamodelselect(U, varargin)
%COPULAMNSIG - uses a model selection algorithm to determine the optimal
%copula to fit the data U based on the algorithm in the paper described
%below.
% Inputs:
%   U - an [M x D] matrix of pseudo observations from which the multinomial
%       signature should be calculated
%   varargin{1} - K, the binning size of the multinomial signature
%                 calculation
%   varargin{2} - a cell array of all copula models to test against, if not
%                 provided, it will be Gaussian, Frank, Gumbel, Clayton
% Outputs:
%   model - the recommended copula model
%   param - the parameter of the copula model
%
% The multinomial signature of a copula is defined in the paper:
%  HELM: Highly Efficient Learning of Mixed copula networks
%  Yaniv Tenzer, Gal Elidan
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
%* along with this program.  If not, see <http://www.gnu.org/licenses/>
%* 
%**************************************************************************

D = size(U,2);

K = 4+1;        % we do a +1 b/c of the way linspace works
copulaFamilyOptions = {'Gaussian', 'Frank', 'Gumbel', 'Clayton'};

nVarargin = length(varargin);
if(nVarargin>0)
    if(isnumeric(varargin{1}))
        K = floor(varargin{1})+1;
    end
end
if(nVarargin>1)
    if(iscell(varargin{2}))
        % TODO: check to make sure each element in the cell array is a
        % valid copula type
        copulaFamilyOptions = varargin{2};
    end
end

minDiv = Inf;
for copulaFamily=copulaFamilyOptions
    famDiv = computeMultivariateDivergence(U, K, copulaFamily);
    if(famDiv < minDiv)
        minDiv = famDiv;
        optimalModel = copulaFamily;
    end
end

% compute the parameters of this optimal model. If it is Gaussian, we will
% get a corr matrix out, so there is no problem.  If it is Archimedean, we
% average all the pairwise dependency parameters as the optimal model
% parameter. This is due to the constraint that the pairwise dependency
% parameters for an Archimedean copula must be the same!
if(strcmpi(optimalModel, 'Gaussian'))
    modelParams = copulafit(optimalModel, U);
else
    alphaVec = zeros(1,nchoosek(D,2));
    alphaVecIdx = 1;
    for jj=1:D
        for ii=1:D
            if(ii<jj)
                U_splice = [U(:,ii) U(:,jj)];
                alphaHat = copulafit(optimalModel, U_splice);
                alphaVec(alphaVecIdx) = alphaHat;
                alphaVecIdx = alphaVecIdx + 1;
            end
        end
    end
    modelParams = mean(alphaVec);
end

end

function [div] = computeMultivariateDivergence(U, K, refCopula)
% computes the divergence between an empirical dataset and its fit to the
% refCopula specified, using the algorithm described in the paper
% referenced above

D = size(U,2);
div = 0;
for jj=1:D
    for ii=1:D
        if(ii<jj)
            U_splice = [U(:,ii) U(:,jj)];
            div = div + computeBivariateDivergence(U_splice, K, refCopula);
        end
    end
end

end

function [div] = computeBivariateDivergence(U, K, refCopula)
% computes the score of the bivariate signature of an empirical dataset
% calculated via the bivariateSignature function against a reference copula
% that has a dependency parameter which is functionally related to the
% empirical spearman's rho of the dataset

% compute the multinomial signature of this reference dataset
empirical_signature = bivariateSignature(U, K);
tau_hat = corr(U, 'type', 'Kendall');
alpha = copulaparam( refCopula, tau_hat(1,2));

% compute the reference signature for this copula
ref_signature_M = 2000;
U_ref = copularnd(refCopula, alpha, ref_signature_M);
ref_signature = bivariateSignature(U_ref, K);

% compare the divergence between the empirical signature and the reference
% signature
div = kldivergence(ref_signature, empirical_signature);

if(isinf(div) || isnan(div))
    error('Error in calculating KL_DIVERGENCE!');
end

end

function [sig] = bivariateSignature(U, K)

ADD_PROB = 0.001;

% Computes the bivariate multinomial signature

edges = ndgrid(linspace(0,1,K));
multinomialBinnedIdxs = discretize(U, edges);

% count how many times each combination occurs
combos = combvec(1:K-1,1:K-1); combos = combos';
numCombos = size(combos,1);
sig = zeros(1,numCombos);
M = size(multinomialBinnedIdxs,1);
for comboIdx=1:numCombos
    combo = combos(comboIdx,:);
    % count number of times this combination occurs
    count = 0;
    for mm=1:M
        if(isequal(multinomialBinnedIdxs(mm,:),combo))
            count = count + 1;
        end
    end
    sig(comboIdx) = count/M;
end

% ensure there is no 0 probability
zeroIdxs = find(sig==0);
nonZeroIdxs = find(sig~=0);

addProb = length(zeroIdxs)*ADD_PROB;
if(addProb>0)
    subtractProb = length(nonZeroIdxs)/addProb;
    sig(zeroIdxs) = ADD_PROB;
    sig(nonZeroIdxs) = sig(nonZeroIdxs) - subtractProb;
end

end