classdef cbn < handle
    %CBN Definition of a Copula Bayesian Network
    %
    %**********************************************************************
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
    %**********************************************************************
    
    properties
        dag;        % the Directed Acyclic Graph structure.  The format of
                    % the DAG is an adjancency matrix.  Element (i,j)
                    % represents a connection from the ith random variable
                    % to the jth random variable
                    % TODO: consider changing this to a sparse matrix
                    % representation
        nodeNames;  % a cell array of the list of node names  
        nodeVals;   % an array of the topological order of the nodes w.r.t
                    % the DAG
        D;          % the number of nodes in the graph
        empInfo;    % a cell array of the empirical information for each
                    % node.  empInfo{i} is the empirical information for
                    % the node nodeNames{nodeVals{i}}.  Empirical
                    % information consits of 3 things, 1.) the domain of
                    % the random variable, the empirical marginal density 
                    % function for node i, and the empirical marginal 
                    % distribution function for node i
        copulaFamilies; % a cell array of the copula of the family associated 
                        % with nodeNames{nodeVals{i}}
        X;          % an M x D matrix of the observable data.  This can be
                    % considered the "training" data.
        disableModelSelectFlag;  % Forces only Gaussian copula's to be used
                                 % by disabling copula model selection
                    
        LOG_CUTOFF;
        PSEUDO_OBS_CALC_METHOD;     % should be either RANK or ECDF
    end
    
    methods
        function obj = cbn(bntPath, X, nodeNames, varargin)
            % CBN - Constructs a CBN object
            %  Inputs:
            %   X - a M x D matrix of the the observable data which the
            %       HCBN will model
            %   nodes - a cell array of names which represent the columns,
            %           where the ith 
            %  Optional Inputs:
            %   dag - an adjacency matrix representing the DAG structure.
            %         If this is not specified, it is assumed that the user
            %         will learn the structure through an available
            %         learning algorithm
            %
            %  TODO
            %   [ ] - 
            %
            
            obj.LOG_CUTOFF = 10^-5;
            obj.PSEUDO_OBS_CALC_METHOD = 'RANK';    % can be either RANK or ECDF
            
            % add BNT to the path
            addpath(genpath(bntPath));
            
            obj.D = size(X,2);      
            
            obj.nodeNames = nodeNames;
            obj.nodeVals = 1:obj.D;
            obj.copulaFamilies = cell(1,obj.D);
            obj.empInfo = cell(1,obj.D);
            
            obj.dag = zeros(obj.D,obj.D);
            
            obj.X = X;
            
            obj.calcEmpInfo();
            obj.disableModelSelectFlag = 0;
            
            nVarargs = length(varargin);
            if(nVarargs>0)
                candidateDag = varargin{1};
                if(~acyclic(candidateDag))
                    error('Specified DAG is not acyclic!\n');
                end
                if(nVarargs>4)
                    opt1 = varargin{2};
                    if(strcmpi(opt1,'copulaFamilyInput'))
                        obj.setDag(candidateDag, 0);    % do NOT compute estimated copula families
                                                        % b/c they were provided as debugging
                                                        % input
                        copulaFamiliesInput = varargin{3};
                        obj.inflateCopulaInformation(copulaFamiliesInput);
                    elseif(strcmpi(opt1,'empInfoInput'))
                        % overwrite the empirical info w/ the provided
                        % debugging empirical info
                        empInfoInput = varargin{3};
                        obj.inflateEmpInfo(empInfoInput);
                    else
                        error('Unrecognized var-input!');
                    end
                    opt2 = varargin{4};
                    if(strcmpi(opt2,'copulaFamilyInput'))
                        obj.setDag(candidateDag, 0);    % do NOT compute estimated copula families
                                                        % b/c they were provided as debugging
                                                        % input
                        copulaFamiliesInput = varargin{5};
                        obj.inflateCopulaInformation(copulaFamiliesInput);
                    elseif(strcmpi(opt2,'empInfoInput'))
                        % overwrite the empirical info w/ the provided
                        % debugging empirical info
                        empInfoInput = varargin{5};
                        obj.inflateEmpInfo(empInfoInput);
                    else
                        error('Unrecognized var-input');
                    end
                elseif(nVarargs>2)
                    opt1 = varargin{2};
                    if(strcmpi(opt1,'copulaFamilyInput'))
                        obj.setDag(candidateDag, 0);    % do NOT compute estimated copula families
                                                        % b/c they were provided as debugging
                                                        % input
                        copulaFamiliesInput = varargin{3};
                        obj.inflateCopulaInformation(copulaFamiliesInput);
                    elseif(strcmpi(opt1,'empInfoInput'))
                        % overwrite the empirical info w/ the provided
                        % debugging empirical info
                        empInfoInput = varargin{3};
                        obj.inflateEmpInfo(empInfoInput);
                        obj.setDag(candidateDag);
                    else
                        error('Unrecognized var-input!');
                    end
                else
                    obj.setDag(candidateDag);       % estimate the copula families
                end
            end
            if(nVarargs>1)
                obj.disableModelSelectFlag = varargin{2};
                obj.setDag(candidateDag);
            end
        end
        
        function [] = inflateCopulaInformation(obj, copulaFamiliesInput)
            for dd=1:obj.D
                if(isempty(copulaFamiliesInput{dd}))
                    obj.copulaFamilies{dd} = [];
                else
                    node = obj.nodeNames{dd};
                    nodeIdx = obj.nodeVals(dd);
                    [parentIdxs, parentNames] = obj.getParents(nodeIdx);
                    C_all = copulaFamiliesInput{dd}{1};
                    C_all_params = copulaFamiliesInput{dd}{2};
                    C_parents = copulaFamiliesInput{dd}{3};
                    C_parents_params = copulaFamiliesInput{dd}{4};
                    
                    copFam = cbnfamily(node, nodeIdx, parentNames, parentIdxs, ...
                        C_all, C_all_params, C_parents, C_parents_params);
                    obj.copulaFamilies{nodeIdx} = copFam;
                end
            end
        end
        
        function [] = inflateEmpInfo(obj, empInfoInput)
            for dd=1:obj.D
                obj.empInfo{dd} = empInfoInput{dd};
            end
        end
        
        function [] = calcEmpInfo(obj)
            %CALCEMPINFO - calculates the empirical distribution function
            %              and the empirical density function via kernel
            %              based methods
            for ii=1:obj.D
                % check if this is a discrete or continuous node for
                % density estimation, we handle these separately
                isdiscrete = 0;
                [F,x] = empcdf(obj.X(:,ii), isdiscrete);
                f = emppdf(obj.X(:,ii), isdiscrete);
                empInfoObj = rvEmpiricalInfo(x, f, F, isdiscrete);
                obj.empInfo{ii} = empInfoObj;
            end
        end
        
        function [parentIdxs, parentNames] = getParents(obj, node)
            %GETPARENTS - returns the indices and the names of a nodes
            %                 parents
            %
            % Inputs:
            %  node - the node index or name of the node for which the
            %         parents are desired
            %
            % Outputs:
            %  parentIdxs - a vector of the indices of all the parents of
            %               this node
            %  parentNames - a cell array of the names of the parents of
            %               this node
            
            if(ischar(node))
                % node name was provided
                nodeName = node;
                for ii=1:obj.D
                    if(isequal(nodeName, obj.nodeNames{ii}))
                        nodeIdx = ii;
                        break;
                    end
                end
            else
                % assume node index was provided
                nodeIdx = node;
            end
            
            % find the node's parents
            parentIdxs = find(obj.dag(:,nodeIdx))';
            parentNames = cell(1,length(parentIdxs));
            for jj=1:length(parentIdxs)
                parentNames{jj} = obj.nodeNames{parentIdxs(jj)};
            end
        end
        
        function [] = setDag(obj, candidateDag, varargin)
            % if you set varargin=1, then the copula families will not be
            % estimated automatically
            obj.dag = candidateDag;
            nVarargin = length(varargin);
            if(nVarargin==0 || varargin{1}~=0)
                obj.estFamilyCopula();
            end
        end
        
        function [] = estFamilyCopula(obj)
            %ESTFAMILYCOPULA - estimates the copulas for each family (child
            %                  and parents), based on the current definition
            %                  of the dag, which is stored in obj.dag
            
            % for each node, estimate the copula of that node and it's
            % parents
            for ii=1:obj.D
                node = obj.nodeNames{ii};
                nodeIdx = obj.nodeVals(ii);
                [parentIdxs, parentNames] = obj.getParents(nodeIdx);
                
                if(isempty(parentIdxs))
                    obj.copulaFamilies{nodeIdx} = [];
                else
                    % grab the appropriate values 
                    X_in = zeros(size(obj.X,1), 1+length(parentNames));
                    X_in(:,1) = obj.X(:,nodeIdx);
                    kk = 2:2+length(parentIdxs)-1;
                    X_in(:,kk) = obj.X(:,parentIdxs);
                    
                    if(strcmpi(obj.PSEUDO_OBS_CALC_METHOD, 'ecdf'))
                        U_in = pobs(X_in, 'ecdf', 100);
                    else
                        U_in = pobs(X_in);
                    end
                    U_parents = U_in(:,2:end);
                    
                    % estimate the best copula for the family
                    if(obj.disableModelSelectFlag || (size(U_in,2)>2))
                        % force the Gaussian Copula here, b/c it is the
                        % only parametric model we have for situations
                        % where D>2
                        C_all = 'Gaussian';
                        C_all_params = copulafit('Gaussian', U_in);
                    else
                        [C_all, C_all_params] = copulamodelselect(U_in);
                    end
                    
                    
                    % estimate the best copula for the parents
                    if(length(parentNames)==1)
                        C_parents = [];
                        C_parents_params = [];
                    else
                        % in order to maintain copula compatibility, we use
                        % the same copula for the parents.  We use
                        % copulafit to find the parameters
                        C_parents = C_all;
                        C_parents_params = copulafit(C_all, U_parents);
                    end
                    copFam = cbnfamily(node, nodeIdx, parentNames, parentIdxs, ...
                        C_all, C_all_params, C_parents, C_parents_params);
                    obj.copulaFamilies{nodeIdx} = copFam;
                end
                
            end
        end
        
        function [ll_val] = dataLogLikelihood(obj, X)
            M = size(X,1);
            if(size(X,2)~=obj.D)
                error('Input data for LL calculation must be the same dimensions as the BN!');
            end      
            ll_val = 0;
            for mm=1:M
                for dd=1:obj.D
                    f_Xi = obj.empInfo{dd}.pdf(X(mm,dd));
                    R_ci = obj.computeCopulaRatio(dd, X(mm,:));
                    
                    if(f_Xi < obj.LOG_CUTOFF)
                        f_Xi = obj.LOG_CUTOFF;
                    end
                    if(R_ci < obj.LOG_CUTOFF)
                        R_ci = obj.LOG_CUTOFF;
                    end
                    
                    ll_val = ll_val + log(f_Xi) + log(R_ci);
                    if(~isreal(ll_val))
                        error('LL Value imaginary!');
                    end
                end
            end
        end
        
        function [ll_val, R_ci_mat] = copulaLogLikelihood(obj, X, pobsFlag)
            M = size(X,1);
            if(size(X,2)~=obj.D)
                error('Input data for LL calculation must be the same dimensions as the BN!');
            end      
            ll_val = 0;
            R_ci_mat = zeros(M,obj.D);
            for mm=1:M
                for dd=1:obj.D
                    R_ci = obj.computeCopulaRatio(dd, X(mm,:), pobsFlag);
                    R_ci_mat(mm,dd) = R_ci;
                    
                    if(isinf(R_ci) || isnan(R_ci))
                        error('R_ci is inf/nan!');
                    end
                    
                    if(R_ci < obj.LOG_CUTOFF)
                        R_ci = obj.LOG_CUTOFF;
                    end
                    
                    ll_val = ll_val + log(R_ci);
                    if(~isreal(ll_val))
                        error('LL Value imaginary!');
                    end
                end
            end
        end
        
        function [rcVal] = computeCopulaRatio(obj, nodeIdx, x, pobsFlag)
            %COMPUTECOPULARATIO - computes the copula ratio as defined by
            %Elidan in Copula Bayesian Networks
            % Inputs:
            %  node - the node for which the copula ratio is to be computed
            %  x - a [1 x D] vector of observations.  This should be all
            %      the observations in the original node order.  This
            %      function will extract the appropriate values from the
            %      vector based on the DAG structure and convert to
            %      pseudo-observations
            %
            % Outputs:
            %  rcVal - the copula ratio value
            
            copFam = obj.copulaFamilies{nodeIdx};
            if(isempty(copFam))
                rcVal = 1;
            else
                idxs_all = [nodeIdx copFam.parentNodeIdxs];
                x_all = x(idxs_all);
                if(pobsFlag)
                    u_all = x_all;
                else
                    u_all = zeros(1,length(x_all));
                    % convert to pseudo-observations via ECDF
                    for ii=1:length(x_all)
                        u_all(ii) = obj.empInfo{idxs_all(ii)}.cdf(x_all(ii));
                    end
                end
                
                u_parents = u_all(2:end);
                if(strcmpi(copFam.C, 'Gaussian'))
                    rcNumVal = copulapdf('Gaussian', u_all, copFam.C_param);
                    if(length(u_parents)>1)
                        rcDenVal = copulapdf('Gaussian', u_parents, copFam.C_parents_param);
                    else
                        rcDenVal = 1;
                    end
                elseif(strcmpi(copFam.C, 'Frank'))
%                     rcNumVal = frankcopulapdf(u_all, copFam.C_param);
                    rcNumVal = copulapdf('Frank', u_all, copFam.C_param);
                    if(length(u_parents)>1)
%                         rcDenVal = frankcopulapdf(u_parents, copFam.C_parents_param);
                        rcDenVal = copulapdf('Frank', u_parents, copFam.C_parents_param);
                    else
                        rcDenVal = 1;
                    end
                elseif(strcmpi(copFam.C, 'Gumbel'))
%                     rcNumVal = gumbelcopulapdf(u_all, copFam.C_param);
                    rcNumVal = copulapdf('Gumbel', u_all, copFam.C_param);
                    if(length(u_parents)>1)
%                         rcDenVal = gumbelcopulapdf(u_parents, copFam.C_parents_param);
                        rcDenVal = copulapdf('Gumbel', u_parents, copFam.C_parents_param);
                    else
                        rcDenVal = 1;
                    end
                elseif(strcmpi(copFam.C, 'Clayton'))
%                     rcNumVal = claytoncopulapdf(u_all, copFam.C_param);
                    rcNumVal = copulapdf('Clayton', u_all, copFam.C_param);
                    if(length(u_parents)>1)
%                         rcDenVal = claytoncopulapdf(u_parents, copFam.C_parents_param);
                        rcDenVal = copulapdf('Clayton', u_parents, copFam.C_parents_param);
                    else
                        rcDenVal = 1;
                    end
                end
                if(~isreal(rcNumVal))
                    if(imag(rcNumVal)<0.001)
                        rcNumVal = real(rcNumVal);
                    else
                        error('rcNumVal has a significant imaginary component!');
                    end
                end
                if(~isreal(rcDenVal))
                    if(imag(rcDenVal)<0.001)
                        rcDenVal = real(rcDenVal);
                    else
                        error('rcNumVal has a significant imaginary component!');
                    end
                end 
                rcVal = rcNumVal/rcDenVal;
            end
        end
        
        function [] = learnStruct_hc(obj, seeddag)
            %LEARNSTRUCT_HC learn the structure of the HCBN network using
            %               the hill climbing algorithm as described in 
            %               Koller and Friedman (2009).  The code for this 
            %               is based off the structure learning toolbox 
            %               within BNT.  See BNT/SLP/learn_struct_hc.m.
            %               This will use the X dataset that was passed to 
            %               the constructor to learn the structure.
            %
            % Inputs:
            %  seeddag - a DAG which can be used as a reference DAG as a
            %            starting point for the search process.  This is
            %            optional, and can be an empty array if no seed is
            %            desired.
            
            % ensure that the seeddag is acyclic
            if(~obj.acyclic(seeddag))
                obj.setDag(zeros(obj.D,obj.D));
            else
                obj.setDag(seeddag);
            end
            
            % get the baseline score
            bestScore = obj.dataLogLikelihood(obj.X);
            done = 0;
            while ~done
                % make dag's which are addition, reversal, and subtraction
                % of edges
                [candidateDags,~,~] = mk_nbrs_of_dag(obj.dag);
                
                % score all the dags
                scores = -Inf*ones(1,length(candidateDags));
                for ii=1:length(candidateDags)
                    obj.setDag(candidateDags{ii});
                    scores(ii) = obj.hcbnLogLikelihood(obj.X);
                end
                
                % find the maximum scoring DAG, and see if it is better
                % than the current best
                maxScore = max(scores);
                % see if multiple candidate dag's had the same maximum
                % score, and if so, choose randomely among those dag's
                new = find(scores == maxScore );
                % update best candidate dag as new dag and continue search
                if ~isempty(new) && (maxScore > bestScore)
                    p = sample_discrete(normalise(ones(1, length(new))));
                    bestScore = maxScore;
                    obj.setDag(candidateDags{new(p)});
                else
                    done = 1;
                end 
            end
            
            % TODO: topo-sort the DAG, and sort the names cell array to
            % match the topologically sorted DAG
            
            % TODO: print out DAG structure
        end
        
    end
end