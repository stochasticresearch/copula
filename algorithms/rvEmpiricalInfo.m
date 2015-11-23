classdef rvEmpiricalInfo
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain;
        density;
        distribution;
    end
    
    methods
        function obj = rvEmpiricalInfo(domain, density, distribution)
            obj.domain = domain;
            obj.density = density;
            obj.distribution = distribution;
        end
        
        function [idx, domainVal] = findClosestDomainVal(q)
            % find the closest value in the domain
            
            % TODO: currently, in discrete RV's it assumes that the query
            %       will be in the domain of the empirical density already
            %       created.  Do we want some error protection here?
            
            if(q > obj.domain(end))
                domainVal = obj.domain(end);
                idx = length(obj.domain);
            elseif(q < obj.domain(1))
                domainVal = obj.domain(1);
                idx = 1;
            else
                tmp = abs(obj.domain-q);
                [~, idx] = min(tmp); %index of closest value
                domainVal = f(idx); %closest value
            end
        end
        
        function [val] = queryDensity(obj, q)
            idx = obj.findClosestDomainVal(q);
            val = obj.density(idx);
        end
        
        function [val] = queryDistribution(obj, q)
            idx = obj.findClosestDomainVal(q);
            val = obj.distribution(idx);
        end
    end
    
end

