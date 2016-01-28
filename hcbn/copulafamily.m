%******************************************************************************
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

classdef copulafamily
    %COPULAFAMILY - has all information related to a copula family
    
    properties
        nodeName;
        nodeIdx;
        parentNodeNames;
        parentNodeIdxs;
        
        type;
        
        C;              % the copula of the node and it's parents
        c;              % the copula density of the node and it's parents
        
        C_parents;      % the copula of the parents (which is the
                        % marginalized coupla of the entire family, with
                        % the first argument = 1)
        c_parents;      % the copula density of the parent's of the node which represents this family
        
        Rho;
        Rho_parents;    
    end
    
    methods
        function obj = copulafamily(nn, ni, pnn, pni, t, CC, cc, RR, ...
                CC_parents, cc_parents, RR_parents)
            obj.nodeName = nn;
            obj.nodeIdx = ni;
            obj.parentNodeNames = pnn;
            obj.parentNodeIdxs = pni;
            
            obj.type = t;
            
            obj.C = CC;
            obj.c = cc;
            
            obj.Rho = RR;
            
            obj.C_parents = CC_parents;
            obj.c_parents = cc_parents;
            obj.Rho_parents = RR_parents;
        end
    end
    
end