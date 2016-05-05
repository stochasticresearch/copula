classdef cbnfamily
    %CBNFAMILY - has all information related to a copula family in the
    %             context of the CBN construction
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
        nodeName;
        nodeIdx;
        parentNodeNames;
        parentNodeIdxs;
        
        C;              % the copula of the node and it's parents
        C_param;        % the parameters of the above copula
                                
        C_parents;      % the copula of the parents (which is the
                        % marginalized coupla of the entire family, with
                        % the first argument = 1)
        C_parents_param;% the parameters of the copula model for the parents
        
    end
    
    methods
        function obj = cbnfamily(nn, ni, pnn, pni, CC, CC_param, CC_parents, CC_parents_param)
            obj.nodeName = nn;
            obj.nodeIdx = ni;
            obj.parentNodeNames = pnn;
            obj.parentNodeIdxs = pni;
            
            obj.C = CC;
            obj.C_param = CC_param;
            
            obj.C_parents = CC_parents;
            obj.C_parents_param = CC_parents_param;
            
        end
    end
    
end