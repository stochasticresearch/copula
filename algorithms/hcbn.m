classdef hcbn
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
    end
    
    methods
        function obj = hcbn(bntPath, D, X, nodes, discreteNodes, varargin)
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
            
            % add BNT to the path
            addpath(genpath(bntPath));
            
            obj.D = D;      
            
            obj.nodeNames = nodes;
            obj.nodeVals = 1:obj.D;     % TODO: need to topologically order the nodeVals
            obj.copulaFamilies = cell(1,obj.D);
            obj.empInfo = cell(1,obj.D);
            
            obj.dag = zeros(D,D);
            
            obj.K = 200;    % hard coded to 200, which seems to be a
                            % reasonable tradeoff between accuracy and
                            % memory/computational requirements
            
            nVarargs = length(varargin);
            if(nVarargs>0)
                obj.dag = varargin{1};
                if(~obj.acyclicCheck())
                    error('Specified DAG is not acyclic!\n');
                end
                % since the DAG is specified, we can estimate the copulas
                obj.estFamilyCopula();
            end
            
            % ensure X is the correct dimensionality
            if(size(X,2)~=D)
                error('Specified data matrix dimensionality does not match D!\n');
            end
            obj.X = X;
            obj.X_xform = X;
            
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
                obj.X_xform(:,idx) = obj.X(:,idx) + (rand(size(X,1),1)-1);
            end
            
            obj.calcEmpInfo();
        end
        
        function [] = calcEmpInfo(obj)
            %CALCEMPINFO - calculates the empirical distribution function
            %              and the empirical density function via kernel
            %              based methods
            for ii=1:obj.D
                % estimate the ECDF
                M = size(obj.X,1);
                [F,x] = ecdf(obj.X(:,ii));
                F = F(2:end);
                x = x(2:end);
                
                % check if this is a discrete or continuous node for
                % density estimation
                if(any(obj.discNodeIdxs==ii))
                    % means this is a discrete node, we handle separately
                    f = zeros(1,length(x));
                    idx = 1;
                    for ii=x
                        f(idx) = sum(obj.X(:,ii)==ii)/M;
                        idx = idx + 1;
                    end
                else
                    f = ksdensity(obj.X(:,ii),x);
                end
                obj.empInfo{ii} = rvEmpiricalInfo(x, f, F);
            end
        end
        
        function [res] = acyclicCheck(obj)
            %ACYCLICCHECK - checks to see if the specified graph is
            %               acyclic.
            %
            % Output:
            %  res - 1 if it is acyclic, 0 if it is not
            %
            % NOTE - this requires code from BNT, so make sure BNT is in
            %        the path
            res = acyclic(obj.dag);
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
            parentIdxs = find(obj.dag(nodeIdx,:));
            parentNames = cell(1,length(parentIdxs));
            for jj=1:parentIdxs
                parentNames{jj} = obj.nodeNames{jj};
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
                nodeIdx = obj.nodeVals{ii};
                
                [parentIdxs, parentNames] = obj.getParents(node);
                
                fprintf('Estimating Copula for Node=%s <-- ', node);
                for jj=1:length(parentNames)
                    fprintf(' %s ', parentNames{jj});
                end
                fprintf('\n');
                
                if(isempty(parentIdxs))
                    % no parents situation
                    obj.copulaFamilies{nodeIdx} = [];
                else
                    % grab the appropriate values 
                    X_in = zeros(size(obj.X_xform,1), 1+length(parentNames));
                    X_in(:,1) = obj.X_xform(:,nodeIdx);

                    kk = 2;
                    for jj=parentIdxs
                        X_in(:,kk) = obj.X_xform(:,jj);
                        kk = kk + 1;
                    end
                    
                    % check to see if any of the indices that we extracted
                    % from were discrete, if so, we calculate the empirical
                    % copula, if not we fit to a copula model (for now only
                    % Gaussian copula)
                    allIdxs = [nodeIdx parentIdxs];
                    if(sum(ismember(allIdxs, obj.discNodeIdxs)))
                        [ C, U, c ] = empcopula(X_in, obj.K);
                        type = 'empirical';
                        Rho = [];
                    else
                        % all continuous marginals, lets fit to a copula
                        % model (for now, we only do Gaussian)
                        type = 'model';
                        C = []; U = []; c = [];
                        u = zeros(size(X_in));
                        M = size(X_in,1);
                        for m=1:M
                            uIdx = 1;
                            for jj=allIdxs
                                u(uIdx) = obj.empInfo{jj}.queryDistribution(X_in(m,uIdx));
                                uIdx = uIdx + 1;
                            end
                        end
                        Rho = copulafit('Gaussian', u);
                    end
                    copFam = copulafamily(node, nodeIdx, parentNames, parentIdxs, ...
                            type, C, U, c, Rho);
                    obj.copulaFamilies{nodeIdx} = copFam;
                end
            end
        end
        
        function [ ll_val ] = hcbnLogLikelihood(obj, X)
            %HCBNLOGLIKELIHOOD - calculates the log-likelihood of the HCBN
            %                    model to the provided data
            % 
            % Inputs
            %  X - the data for which to calculate the log-likelihood for
            ll_val = 0;
            for ii=1:obj.D
                ll_val = ll_val + obj.copulall(ii,X);
            end
        end
        
        function [ ll_val ] = copulall(obj, nodeIdx, X )
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
            [parentIdxs, parentNodes] = obj.getParents(nodeIdx);
            % grab the appropriate values 
            X_in = zeros(size(X,1), 1+length(parentNames));
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
            for m=1:M
                ll_val = 0;
                % compute empirical distribution log likelihood
                colIdx = 1;
                for jj=allIdxs
                    ll_val = ll_val + log(obj.empInfo{jj}.queryDensity(X_in(m,colIdx)));
                    colIdx = colIdx + 1;
                end
                
                % generate u
                uIdx = 1;
                for jj=allIdxs
                    u(uIdx) = obj.empInfo{jj}.queryDistribution(X_in(m,uIdx));
                    uIdx = uIdx + 1;
                end
                
                % compute copula ratio value
                ll_val = ll_val + log(obj.copularatio_val(nodeIdx, u));
            end
        end
        
        function [c_val] = copularatio_val(obj, nodeIdx, u)
            %COPULARATIO_VAL - calculates the copula ratio for a node at a
            %                  location in the unit hypercube
            % Inputs:
            %  nodeIdx - the node for which the copula ratio is to be
            %            calculated.  This is the node index.
            %  u - a vector of a point in the unit-hypercube where
            %      the copula ratio will be calculated
            
            fprintf('Calculating Copula Ratio for Node %s\n', ...
                obj.nodeNames{nodeIdx});
            
            % find the associated copula family
            copFam = obj.copulaFamilies{nodeIdx};
            
            if(isempty(copFam))
                % if copFam is empty, this means this node has no parents
                % and thus by defintion the copula ratio here is defined to
                % be 1
                c_val = 1;
            elseif(strcmp(copFam.type, 'empirical'))
                % get the copula density for this family
                C = copFam.C; c = copFam.c{end};

                % query it for the specified point with empcopula_val
                [~, c_val_numerator] = empcopula_val(C,c,u);
                u_denom = u; u_denom(1) = 1;
                [~, c_val_denominator] = empcopula_val(C,c,u_denom);
                c_val = c_val_numerator/c_val_denominator;
            else
                % assume we fit the Gaussian model to the data
                c_val_numerator = copulapdf('Gaussian', u, copFam.Rho);
                u_denom = u; u_denom(1) = 1;
                c_val_denominator = copulapdf('Gaussian', u_denom, copFam.Rho);
                c_val = c_val_numerator/c_val_denominator;
            end
        end
    end
    
end

