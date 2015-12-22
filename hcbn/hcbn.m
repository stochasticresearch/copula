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
                        
        dagLL_scalingFactor;    % a scaling factor for the log-likelihood
                                % of the DAG with the training data.  This
                                % scaling factor will make the likelihood
                                % 1, and when we compare this against
                                % "test" datasets, we will see how well the
                                % model holds against the other datasets
        
                                
                                
        DEBUG_MODE; % if turned on, will print out extra stuff to screen to
                    % monitor what is going on
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
            
            obj.DEBUG_MODE = 1;
            
            % add BNT to the path
            addpath(genpath(bntPath));
            
            obj.D = size(X,2);      
            
            obj.nodeNames = nodes;
            obj.nodeVals = 1:obj.D;
            obj.copulaFamilies = cell(1,obj.D);
            obj.empInfo = cell(1,obj.D);
            
            obj.dag = zeros(obj.D,obj.D);
            
            obj.K = 25;     % hard coded to 25, which seems to be a
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
                    for jj=1:length(x)
                        f(idx) = sum(obj.X(:,ii)==x(jj))/M;
                        idx = idx + 1;
                    end
                else
                    f = ksdensity(obj.X(:,ii),x);
                end
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
        
        function [] = setDag(obj, candidateDag)
            obj.dag = candidateDag;
            obj.estFamilyCopula();
            llTrain = obj.hcbnLogLikelihood(obj.X);
            obj.dagLL_scalingFactor = llTrain;
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
                    ecdfNumPts = 100;
                    X_in = zeros(size(obj.X_xform,1), 1+length(parentNames));
                    FX_in = zeros(ecdfNumPts, 1+length(parentNames));
                    X_in(:,1) = obj.X_xform(:,nodeIdx);
                    
                    kk = 2:2+length(parentIdxs)-1;
                    X_in(:,kk) = obj.X_xform(:,parentIdxs);
                    
                    % generate pseudo-observations
                    for jj=1:size(X_in,2)
                        FX_in(:,jj) = ksdensity(X_in(:,jj),...
                            linspace(min(X_in(:,jj)),max(X_in(:,jj)),ecdfNumPts),...
                            'function','cdf')';
                    end
                    
                    % check to see if any of the indices that we extracted
                    % from were discrete, if so, we calculate the empirical
                    % copula, if not we fit to a copula model (for now only
                    % Gaussian copula)
                    allIdxs = [nodeIdx parentIdxs];
                    if(sum(ismember(allIdxs, obj.discNodeIdxs)))
                        M = size(X_in,1);
                        C = [];         % TODO: estimate this, but not necessary for hcbn
                        c = empcopuladensity(FX_in, obj.h, obj.K, 'betak');
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
                            FX_in_parents = FX_in(:,2:end);
                            C_parents = []; % TODO: estimate this, but not necessary for hcbn
                            c_parents = empcopuladensity(FX_in_parents, obj.h, obj.K, 'betak');
                        end
                    else
                        % all continuous marginals, lets fit to a copula
                        % model (for now, we only do Gaussian)
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
        
        function [ ll_val ] = hcbnLogLikelihood(obj, X)
            %HCBNLOGLIKELIHOOD - calculates the log-likelihood of the HCBN
            %                    model to the provided data
            % 
            % Inputs
            %  X - the data for which to calculate the log-likelihood for
            if(isempty(obj.dag))
                ll_val = -Inf;
            else
                ll_val = 0;
                for ii=1:obj.D
                    ll_val = ll_val + obj.copulall(ii,X);
                end
            end
        end
        
        function [ ll_val, Rc_val_vec ] = copulall(obj, nodeIdx, X )
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
            end
            %%%%%%%%%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%
            
            for m=1:M
                % compute empirical distribution log likelihood (for just
                % the child node)
                ll_val = ll_val + log(obj.empInfo{nodeIdx}.queryDensity(X_in(m,1)));
                
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
                Rc = obj.copulaRatioVal(nodeIdx, u);
                
                %%%%%%%%%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%
                if(nargout>1)
                    Rc_val_vec(m) = Rc; % for DEBUG purposes
                end
                %%%%%%%%%%%%%%%%%%%% DEBUGGING %%%%%%%%%%%%%%%%%%%%%%%
                
                ll_val = ll_val + log(Rc);
            end
        end
        
        function [Rc_val] = copulaRatioVal(obj, nodeIdx, u)
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
            bestScore = obj.hcbnLogLikelihood(obj.X);
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

