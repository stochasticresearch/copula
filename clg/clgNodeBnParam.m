classdef clgNodeBnParam < handle
    properties
        node;
        combo;
        Mean;
        Covariance;
    end
    
    methods
        function obj = clgNodeBnParam(node, combo, Mean, Covariance)
            obj.node = node;
            obj.combo = combo;
            obj.Mean = Mean;
            obj.Covariance = Covariance;
        end
    end
    
end