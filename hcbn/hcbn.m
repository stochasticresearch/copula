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

classdef hcbn < handle
    %HCBN Definition of a Hybrid Copula Bayesian Network
    
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
                    % the random variable, the empirical density function,
                    % and the empirical distribution function
        copulaFamilies; % a cell array of the copula of the family associated 
                        % with nodeNames{nodeVals{i}}
        X;          % an N x D matrix of the observable data.  This can be
                    % considered the "training" data.
        X_xform;    % a transformed N x D matrix of the observable data,
                    % where we "continue" the discrete random variables, as
                    % described in Neslehova's paper, and Denuit and
                    % Lambert's paper
        K;          % the number of points over the unit-hypercube for
                    % which to calculate the empirical copula.  The higher
                    % the value of K, the greater the accuracy, but also
                    % the higher the memory and processing time
                    % requirements
        discreteNodes;  % the names of the nodes which should be modeled as
                        % discrete random variables
        discNodeIdxs;   % the indices of of the nodes of the discrete
                        % random variables
        
        h;          % the beta-kernel density estimation weighting
                        
        LOG_CUTOFF;
                                
        DEBUG_MODE; % if turned on, will print out extra stuff to screen to
                    % monitor what is going on and store extra debugging
                    % variables (slower)
        SIM_NUM;    % used for debugging
        TYPE_NAME;  % used for debugging
    end
    
    methods
        function obj = hcbn(bntPath, X, nodes, discreteNodes, varargin)
            % HCBN - Constructs a HCBN object
            %  Inputs:
            %   X - a N x D matrix of the the observable data which the
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
            
            obj.DEBUG_MODE = 0;
            obj.SIM_NUM = 0;
            obj.TYPE_NAME = '';
            obj.LOG_CUTOFF = 10^-5;
            
            % add BNT to the path
            addpath(genpath(bntPath));
            
            obj.D = size(X,2);      
            
            obj.nodeNames = nodes;
            obj.nodeVals = 1:obj.D;
            obj.copulaFamilies = cell(1,obj.D);
            obj.empInfo = cell(1,obj.D);
            
            obj.dag = zeros(obj.D,obj.D);
            
            obj.K = 25;    % hard coded, this value seems to be a
                            % reasonable tradeoff between accuracy and
                            % memory/computational requirements
            
            obj.X = X;
            obj.X_xform = X;

            obj.h = 0.05;
            
            obj.discreteNodes = discreteNodes;
            obj.discNodeIdxs = zeros(1,length(obj.discreteNodes));
            for ii=1:length(obj.discNodeIdxs)
                nodeName = obj.discreteNodes{ii};
                for jj=1:length(obj.nodeNames)
                    if(isequal(nodeName,obj.nodeNames{jj}))
                        obj.discNodeIdxs(ii) = jj;
                        break;
                    end
                end
            end
            
            % continue the discrete random variables.  X* = X + (U - 1)
            for idx=obj.discNodeIdxs
                obj.X_xform(:,idx) = continueRv(obj.X(:,idx));
            end
            
            obj.calcEmpInfo();
            
            nVarargs = length(varargin);
            if(nVarargs>0)
                candidateDag = varargin{1};
                if(~acyclic(candidateDag))
                    error('Specified DAG is not acyclic!\n');
                end
                obj.setDag(candidateDag);
            end
        end
        
        function [] = setSimNum(obj, simnum)
            obj.SIM_NUM = simnum;
        end
        
        function [] = setTypeName(obj, typeName)
            obj.TYPE_NAME = typeName;
        end
        
        function [] = calcEmpInfo(obj)
            %CALCEMPINFO - calculates the empirical distribution function
            %              and the empirical density function via kernel
            %              based methods
            for ii=1:obj.D
                % check if this is a discrete or continuous node for
                % density estimation, we handle these separately
                isdiscrete = 0;
                if(any(obj.discNodeIdxs==ii))
                    isdiscrete = 1;
                end
                [F,x] = empcdf(obj.X(:,ii), isdiscrete);
                f = emppdf(obj.X(:,ii), isdiscrete);
                empInfoObj = rvEmpiricalInfo(x, f, F);
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
            % if you set varargin=1, then 
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
                
                if(obj.DEBUG_MODE)
                    fprintf('Estimating Copula for Node=%s:%d <-- ', node, nodeIdx);
                    for jj=1:length(parentNames)
                        fprintf('%d:%s ', parentIdxs(jj), parentNames{jj});
                    end
                    fprintf('\n');
                end
                
                if(isempty(parentIdxs))
                    % no parents situation
                    obj.copulaFamilies{nodeIdx} = [];
                else
                    % grab the appropriate values 
                    X_in = zeros(size(obj.X_xform,1), 1+length(parentNames));
                    X_in(:,1) = obj.X_xform(:,nodeIdx);
                    kk = 2:2+length(parentIdxs)-1;
                    X_in(:,kk) = obj.X_xform(:,parentIdxs);
                    % check to see if any of the indices that we extracted
                    % from were discrete, if so, we calculate the empirical
                    % copula, if not we fit to a copula model (for now only
                    % Gaussian copula)
                    allIdxs = [nodeIdx parentIdxs];
                    if(sum(ismember(allIdxs, obj.discNodeIdxs)))
                        U_in = pseudoobs(X_in);     % WARNING: before, it used to be pseudoobs(X_in, 'ecdf', 100);
                        
                        M = size(X_in,1);
                        C = [];         % TODO: estimate this, but not necessary for hcbn
                        c = empcopulapdf(U_in, obj.h, obj.K, 'betak');
                        type = 'empirical';
                        Rho = [];
                        Rho_parents = [];
                        if(length(parentIdxs)==1)
                            % the density will be used directly, so there
                            % is no need to calculate these copula's (they
                            % don't actually make sense because there is
                            % only one dimension, no concept of "joint"
                            % distribution)
                            C_parents = [];
                            c_parents = [];
                        else
                            U_in_parents = U_in(:,2:end);
                            C_parents = []; % TODO: estimate this, but not necessary for hcbn
                            c_parents = empcopulapdf(U_in_parents, obj.h, obj.K, 'betak');
                        end
                    else
                        % all continuous marginals, lets fit to a copula
                        % model (for now, we only do Gaussian)
                        % TODO: use other models as well
                        type = 'model';
                        C = []; c = [];
                        C_parents = []; c_parents = [];
                        u = zeros(size(X_in));
                        M = size(X_in,1);
                        for m=1:M
                            uIdx = 1;
                            for jj=allIdxs
                                u(m,uIdx) = obj.empInfo{jj}.queryDistribution(X_in(m,uIdx));
                                uIdx = uIdx + 1;
                            end
                        end
                        Rho = copulafit('Gaussian', u);
                        if(length(parentIdxs)==1)
                            Rho_parents = [];
                        else
                            u_parents = u(:,2:end);
                            Rho_parents = copulafit('Gaussian', u_parents);
                        end
                    end
                    copFam = copulafamily(node, nodeIdx, parentNames, parentIdxs, type, ...
                            C, c, Rho, C_parents, c_parents, Rho_parents);
                    obj.copulaFamilies{nodeIdx} = copFam;
                end
            end
        end
        
        function [X] = genFamilySamples(obj, node, M)
            %GENFAMILYSAMPLES - generates samples from a specified family
            % Inputs:
            %  node - the child node, whose family is sampled
            %  M - the number of samples to generate
            % Outputs:
            %  X - the samples from the family associated with the input
            %      node.  The columns of X are [node, parent1, parent2 ...]
            
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
            
            % generate U from the copula that was "learned"
            copFam = obj.copulaFamilies{nodeIdx};
            if(strcmpi(copFam.type, 'empirical'))
                U = empcopularnd(copFam.c, M);
            else
                % TODO: when other copula models are introduced into the
                % continuous case, we need to change copularnd here to
                % incorporate the other copula types
                U = copularnd('Gaussian', copFam.Rho, M);
            end
            D_family = size(U,2);
            X = zeros(size(U));
            % invert each U appropriately to generate X
            allIdxs = [nodeIdx copFam.parentNodeIdxs];
            for ii=1:M
                for dd=1:D_family
                    X(ii,dd) = obj.empInfo{allIdxs(dd)}.invDistribution(U(ii,dd));
                end
            end
        end
        
        function [ ll_val, Rc_val_mat, Rc_num_mat, Rc_den_mat ] = dataLogLikelihood(obj, X)
            %DATALOGLIKELIHOOD - calculates the log-likelihood of the HCBN
            %                    model to the provided data
            % 
            % Inputs
            %  X - the data for which to calculate the log-likelihood for
			localFunctionDebugMode = 0;
			if(nargout>1)
				localFunctionDebugMode = 1;
				Rc_val_mat = zeros( obj.D, size(X,1) );
				Rc_num_mat = zeros( obj.D, size(X,1) );
				Rc_den_mat = zeros( obj.D, size(X,1) );
			end
			
			ll_val = 0;
			
			for ii=1:obj.D
				if(localFunctionDebugMode)
					[tmp, Rc_val_vec, Rc_num_vec, Rc_den_vec] = obj.copulall(ii,X);
					Rc_val_mat(ii,:) = Rc_val_vec;
					Rc_num_mat(ii,:) = Rc_num_vec;
					Rc_den_mat(ii,:) = Rc_den_vec;
				else
					tmp = obj.copulall(ii,X);
				end
				ll_val = ll_val + tmp;                
			end
        end
        
        function [ ll_val, Rc_val_vec, Rc_num_vec, Rc_den_vec ] = copulall(obj, nodeIdx, X )
            %COPULALL Computes the log-likelihood of a copula density matching the
            %         given data
            %
            % Inputs:
            %  nodeIdx - the node for which to calculate the copula
            %            likelihood
            %  X - the data against which to calculate the likelihood. This
            %      should be the full data (i.e. all the columns for all
            %      the random variables).  The correct columns will be
            %      extracted by the copulall function
            %
            % Outputs:
            %  ll_val - the log likelihood, defined as follows:
            %           ll_val = sum(1,M, log(c(u_1[m], ..., u_n[m])/c(1,u_2[m], ..., u_n[m])) +
            %                    sum(1,M, log(f(x_1[m]))
            
            % get the parents associated w/ this node
            parentIdxs = obj.getParents(nodeIdx);
            % grab the appropriate values 
            X_in = zeros(size(X,1), 1+length(parentIdxs));
            X_in(:,1) = X(:,nodeIdx);

            kk = 2;
            for jj=parentIdxs
                X_in(:,kk) = X(:,jj);
                kk = kk + 1;
            end
            allIdxs = [nodeIdx parentIdxs];
            
            % compute the copularatio for each data point
            M = size(X_in,1);
            ll_val = 0;
            u = zeros(1, size(X_in,2));
            %%%%%%%%%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%
            if(nargout>1)
                Rc_val_vec = zeros(1,M);
                Rc_num_vec = zeros(1,M);
                Rc_den_vec = zeros(1,M);
            end
            %%%%%%%%%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%
            
            for m=1:M
                % compute empirical distribution log likelihood (for just
                % the child node)
                f_xi = obj.empInfo{nodeIdx}.queryDensity(X_in(m,1));
                if(f_xi<=obj.LOG_CUTOFF)     % do this for numerical reasons
                    f_xi = obj.LOG_CUTOFF;
                end
                ll_val = ll_val + log(f_xi);
                
                % generate u, avoid values of u being 0 or 1
                uIdx = 1;
                for jj=allIdxs
                    uu = obj.empInfo{jj}.queryDistribution(X_in(m,uIdx));
                    if(abs(uu-1)<=.001)
                        uu = 0.999;
                    elseif(abs(uu-.001)<=0.001)
                        uu = 0.001;
                    end
                    u(uIdx) = uu;
                    uIdx = uIdx + 1;
                end
                
                % compute copula ratio value
                [Rc, num_val, den_val] = obj.copulaRatioVal(nodeIdx, u);
                if(Rc<=obj.LOG_CUTOFF)
                    Rc = obj.LOG_CUTOFF;
                end 
                %%%%%%%%%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%
                if(nargout>1)
                    Rc_val_vec(m) = Rc;
                    Rc_num_vec(m) = num_val;
                    Rc_den_vec(m) = den_val;
                end
                ll_val_prev = ll_val;
                %%%%%%%%%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%
                
                ll_val = ll_val + log(Rc);
                if(isnan(ll_val) || isinf(ll_val))
                    warning('Wrote DEBUG INFORMATION!');
                    1;      % PUT DEBUG ON AT THIS LINE FOR INTERACTIVE DEBUGGING
                    fid = fopen(sprintf('/home/kiran/ownCloud/PhD/sim_results/copulall_%d.txt', obj.SIM_NUM), 'a');
                    fprintf(fid, 'sample=%d Rc_val=%f Rc_num=%f Rc_den=%f nodeIdx=%d ll_val_prev=%f u=%s \n', ...
                        m, Rc, num_val, den_val, nodeIdx, ll_val_prev, sprintf('%f,', u) );
                    fclose(fid);
                    % save all the workspace variables for debugging
                    save(sprintf('/home/kiran/ownCloud/PhD/sim_results/%s_%d.mat', ...
                        obj.TYPE_NAME, obj.SIM_NUM));
                end
            end
        end
        
        function [Rc_val, c_val_numerator, c_val_denominator] = copulaRatioVal(obj, nodeIdx, u)
            %COPULARATIOVAL - calculates the copula ratio for a node at a
            %                 location in the unit hypercube
            % Inputs:
            %  nodeIdx - the node for which the copula ratio is to be
            %            calculated.  This is the node index.
            %  u - a vector of a point in the unit-hypercube where
            %      the copula ratio will be calculated
            
            % find the associated copula family
            copFam = obj.copulaFamilies{nodeIdx};
            
            if(isempty(copFam))
                % if copFam is empty, this means this node has no parents
                % and thus by defintion the copula ratio here is defined to
                % be 1
                c_val_numerator = 1;
                c_val_denominator = 1;
            elseif(strcmp(copFam.type, 'empirical'))
                % get the copula density for this family
                c = copFam.c;
                c_parents = copFam.c_parents;

                % query it for the specified point with empcopula_val
                c_val_numerator = empcopula_val(c,u);
                u_denom = u(2:end);
                if(length(u_denom)==1)
                    c_val_denominator = 1;
                else
                    c_val_denominator = empcopula_val(c_parents,u_denom);
                end
            else
                % assume we fit the Gaussian model to the data
                c_val_numerator = copulapdf('Gaussian', u, copFam.Rho);
                
                u_denom = u(2:end);
                if(length(u_denom)==1)
                    c_val_denominator = 1;
                else
                    c_val_denominator = copulapdf('Gaussian', u_denom, copFam.Rho_parents);
                end
            end
            
            Rc_val = c_val_numerator/c_val_denominator;
        end
        
    end
end
