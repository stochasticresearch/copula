function [optimOut,fval] = maxcopfamprob(copModel,copParams,x,rvEmpInfoObjs,missingIdxs,intIdxs)
%OPTIMFAMPROB Maximizes a copula family's probability, given the copula
%model, the data, and the missing data indices
% Inputs:
%  copModel - the copula model which describes the dependency structure for
%             this particular copula family.  Valid options are:
%             1.) frank
%             2.) gumbel
%             3.) clayton
%             4.) gaussian
%             6.) empirical (not yet supported)
%  copParams - the parameters associated with this copula model.  This
%              should be a list of parameters which describe the copula,
%              for example, a Rho matrix for a Gaussin copula, etc.. For
%              Archimedean copulas -- please be careful! The parameter
%              should be according to Mathwork's convention 1/alpha = theta
%              where theta is the dependence parameter of the copula as
%              defined by virtually everyone else in the world ... 
%  x - a vector which has the observed and unobserved nodes.  The order in
%      which these elements appear in each row should match the dimensions 
%      of the copula.  Probability is maximized for each row.  For missing
%      values, the value provided in x for the missing index will be
%      replaced with the expected value of the empirical probability
%      distribution representing that index's marginal distribution.  
%  rvEmpInfoObjs - a cell array of rvEmpiricalInfo objects.  The
%             order must match the order of arguments for x.
%  missingIdxs - indices of x which are the missing data.
%  intIdxs - indices of x which are integer bound (for hybrid or discrete
%        networks)
%
% Outputs:
%  
%
% Acknowledgements:
% This code modeled after the code in the CoImp R package.  See
% specifically: https://github.com/cran/CoImp/blob/master/R/CoImp.R
%
%**************************************************************************
%* 
%* Copyright (C) 2017  Kiran Karra <kiran.karra@gmail.com>
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
D = length(x);
xx = x;
% replace the missing idx's with their mean
for ii=missingIdxs
    xx(ii) = rvEmpInfoObjs{ii}.mean;
end
remainingIdxs = setdiff(1:D,missingIdxs);

% put the variables to be optimized into their own vector
optimArgs = xx(missingIdxs);
% put the remaining parameters needed by the optimization routine into
% their own vector
params = xx(remainingIdxs);
f = @(zz)optimf(zz,params,missingIdxs,remainingIdxs,...
                rvEmpInfoObjs,copModel,copParams);

lb = zeros(1,length(missingIdxs)); ub = zeros(1,length(missingIdxs));
for ii=1:length(missingIdxs)
    jj = missingIdxs(ii);
    lb(ii) = min(rvEmpInfoObjs{jj}.domain);
    ub(ii) = max(rvEmpInfoObjs{jj}.domain);
end
            
% determine if any of the missing indices are integer, if so, we use the
% mixed integer optimization provided by "ga", otherwise, we use
% "fminsearch"
mixedOpt = any(ismember(missingIdxs,intIdxs));
if(mixedOpt)
    % determine which of the indices we should constrain to integer values
    IntCon = zeros(1,length(missingIdxs));
    for ii=1:length(missingIdxs)
        if(any(ismember(missingIdxs(ii),intIdxs)))
            IntCon(ii) = ii;
        end
    end
    IntCon(IntCon==0) = [];
    [optimOut,fval] = ga(f,length(missingIdxs),[],[],[],[],...
                  lb,ub,[],IntCon);
else
    if(length(missingIdxs)==1)
        [optimOut,fval] = fminbnd(f,lb,ub);
    else
        [optimOut,fval] = fmincon(f,optimArgs,[],[],[],[],lb,ub);
    end
    
end

end

% create the optimization function.
function val = optimf(optimArgs,params,missingIdxs,remainingIdxs,...
                      rvEmpInfoObjs,copModel,copParams)
    % join optimArgs & params according to information in missingIdxs and
    % remainingIdxs
    idxsJoined = [missingIdxs remainingIdxs];
    [~,I] = sort(idxsJoined);
    xxJoined = [optimArgs params];
    xxJoined = xxJoined(I);
    
    val = copfamprob(xxJoined,rvEmpInfoObjs,copModel,copParams)*-1;
end