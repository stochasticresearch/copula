classdef copulafamily
    %COPULAFAMILY - has all information related to a copula family
    
    properties
        nodeName;
        nodeIdx;
        parentNodeNames;
        parentNodeIdxs;
        
        type;
        
        C;
        U;
        c;
        
%         D;              % the total dimensionality of the copula family
        
        C_parents;      % the copula of the parents (which is the
                        % marginalized coupla of the entire family, with
                        % the first argument = 1)
        c_parents;      % the copula density of C_parents
        U_parents;      % the grid points in the unit hypercube over which
                        % the parents copula function and copula density
                        % are defined
        
        Rho;
        Rho_parents;    
    end
    
    methods
        function obj = copulafamily(nn, ni, pnn, pni, t, CC, UU, cc, RR, ...
                CC_parents, UU_parents, cc_parents, RR_parents)
            obj.nodeName = nn;
            obj.nodeIdx = ni;
            obj.parentNodeNames = pnn;
            obj.parentNodeIdxs = pni;
            
            obj.type = t;
            
            obj.C = CC;
            obj.U = UU;
            obj.c = cc;
            
            obj.Rho = RR;
            
%             obj.D = 1+length(obj.parentNodeIdxs);
            
            % default definitions, if the type is 'model', then we can
            % evaluate these values using built in matlab functions, so
            % these don't need to be updated
            obj.C_parents = CC_parents;
            obj.c_parents = cc_parents;
            obj.U_parents = UU_parents;
            obj.Rho_parents = RR_parents;
            
%             if(strcmp(obj.type,'empirical'))
%                 % marginalize C
%                 sz = size(obj.C);
%                 if(length(sz)==2)
%                     new_sz = [1 sz(end)];
%                 else
%                     new_sz = sz(2:end);
%                 end;
%                 % we will reshape these arrays 
%                 obj.C_parents = zeros(new_sz);
%                 
%                 % compute DF of marginalized C_parents
%                 u_parents = zeros(1, length(obj.D-1));  %[ jj kk ll mm ...], jj = dim2, kk=dim3 ...
%                 for ii=1:length(obj.C_parents)
%                     % generate u_query
%                     remVal = ii - 1;
%                     divVal = prod(new_sz(1:end-1));
%                     for jj=1:length(new_sz)
%                         u_parents(jj) = floor( remVal/divVal );
%                         remVal = remVal - u_parents(jj)*divVal;
%                         divVal = divVal / new_sz(1);
%                     end
%                     u_parents = u_parents + 1;
%                     u_parents = u_parents / new_sz(end);
%                     u_query = [1 u_parents];
%                     
%                     % get the value of the copula function at this u
%                     [C_u] = empcopula_val(obj.C, obj.c{end}, u_query);
%                     
%                     % store into C_parents (using linear index format)
%                     obj.C_parents(ii) = C_u;
%                 end
%                 
%                 % differentiate to get c_parents
%                 idxShiftArr = 1:obj.D-1;
%                 shiftAmt = 0;
%                 for ii=1:obj.D-1
%                     obj.c_parents = gradient(permute(obj.C_parents, circshift(idxShiftArr',shiftAmt)'));
%                     shiftAmt = shiftAmt + 1;
%                 end                
%             end
        end
    end
    
end