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
        
        Rho;
    end
    
    methods
        function obj = copulafamily(nn, ni, pnn, pni, t, CC, UU, cc, RR)
            obj.nodeName = nn;
            obj.nodeIdx = ni;
            obj.parentNodeNames = pnn;
            obj.parentNodeIdxs = pni;
            
            obj.type = t;
            
            obj.C = CC;
            obj.U = UU;
            obj.c = cc;
            
            obj.Rho = RR;
        end
    end
    
end