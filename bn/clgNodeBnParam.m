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

classdef clgNodeBnParam < handle
    properties
        node;
        parentCombination;
        parentCombinationProbability;
        Mean;
        Covariance;
    end
    
    methods
        function obj = clgNodeBnParam(node, combo, comboProb, Mean, Covariance)
            obj.node = node;
            obj.parentCombination = combo;
            obj.parentCombinationProbability = comboProb;
            obj.Mean = Mean;
            obj.Covariance = Covariance;
        end
    end
    
end