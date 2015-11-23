classdef hcbn
    %HCBN Definition of a Hybrid Copula Bayesian Network
    
    properties
        dag;        % the Directed Acyclic Graph structure.  The format of
                    % the DAG is an adjancency matrix.  Element (i,j)
                    % represents a connection from the ith random variable
                    % to the jth random variable
        nodeNames;  % a cell array of the list of node names  
        nodeVals;   % an array of the topological order of the nodes w.r.t
                    % the DAG
        D;          % the number of nodes in the graph
        empDists;   % a cell array of the empirical distributions for each
                    % node.  empDists{i} is the empirical distribution of
                    % the node nodeNames{nodeVals{i}}
        copulaFamilies;  % a cell array of the copula of the family associated 
                         % with nodeNames{nodeVals{i}}
        X;          % an N x D matrix of the observable data
        K;          % the number of points over the unit-hypercube for
                    % which to calculate the empirical copula.  The higher
                    % the value of K, the greater the accuracy, but also
                    % the higher the memory and processing time
                    % requirements
    end
    
    methods
        function obj = hcbn(D, X, nodes, varargin)
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
            obj.D = D;      
            
            obj.nodeNames = nodes;
            obj.empDists = cell(1,obj.D);
            obj.nodeVals = zeros(1,obj.D);
            obj.copulaFamilies = cell(1,obj.D);
            
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
            
        end
        
        function [] = calcMarginalDist(obj)
            %CALCMARGINALDIST - calculates the empirical distributions of
            %                   univariate marginals in the dataset
            % Outputs:
            %  N/A - the object is updated
        end
        
        function [res] = acyclicCheck(obj)
            %ACYCLICCHECK - checks to see if the specified graph is
            %               acyclic.
            %
            % Output:
            %  res - 1 if it is acyclic, 0 if it is not
            
            % TODO: actually check the adjacency matrix
            res = 1;
        end
        
        function [] = estFamilyCopula(obj)
            %ESTFAMILYCOPULA - estimates the copulas for each family (child
            %                  and parents), based on the current definition
            %                  of the dag, which is stored in obj.dag
            
            % for each node, estimate the copula of that node and it's
            % parents
            for ii=1:D
                node = obj.nodeNames{ii};
                nodeIdx = obj.nodeVals{ii};
                % find the node's parents
                parentIdxs = find(obj.dag(nodeIdx,:));
                parentNames = cell(1,length(parentIdxs));
                for jj=1:parentIdxs
                    parentNames{jj} = obj.nodeNames{jj};
                end
                
                fprintf('Estimating Copula for Node=%s <-- ', node);
                for jj=1:length(parentNames)
                    fprintf(' %s ', parentNames{jj});
                end
                fprintf('\n');
                
                % grab the appropriate values 
                X_in = zeros(size(obj.X,1), 1+length(parentNames));
                X_in(:,1) = X(:,nodeIdx);
                kk = 2;
                for jj=parentIdxs
                    X_in(:,kk) = X(:,jj);
                    kk = kk + 1;
                end
                
                [ C, U, c ] = empcopula(X_in, obj.K);
                copFam = copulafamily(node, nodeIdx, parentNames, parentIdxs, C, U, c);
                obj.copulaFamilies{nodeIdx} = copFam;
            end
            
        end
        
        function [c_val] = copularatio_val(obj, nodeIdx, u)
            %COPULARATIO_VAL - calculates the copula ratio for a node at a
            %                  location in the unit hypercube
            % Inputs:
            %  nodeIdx - the node for which the copula ratio is to be
            %            calculated.  This is the node index.
            %  u - a column vector of a point in the unit-hypercube where
            %      the copula ratio will be calculated
            
            fprintf('Calculating Copula Ratio for Node %s\n', ...
                obj.nodeNames{nodeIdx});
            
            % find the associated copula family
            copFam = obj.copulaFamilies{nodeIdx};
            
            % get the copula density for this family
            C = copFam.C; c = copFam.c{end};
            
            % query it for the specified point with empcopula_val
            [~, c_val_numerator] = empcopula_val(C,c,u);
            u_denom = u; u_denom(1) = 1;
            [~, c_val_denominator] = empcopula_val(C,c,u_denom);
            c_val = c_val_numerator/c_val_denominator;
        end
    end
    
end

