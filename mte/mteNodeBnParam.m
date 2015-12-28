classdef mteNodeBnParam < handle
    properties
        node;
        combo;
        mte_info;
    end
    
    methods
        function obj = mteNodeBnParam(node, combo, mte_info)
            obj.node = node;
            obj.combo = combo;
            obj.mte_info = mte_info;
        end
    end
    
end