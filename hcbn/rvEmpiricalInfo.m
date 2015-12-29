classdef rvEmpiricalInfo
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        domain;
        density;
        distribution;
    end
    
    methods(Static)
        function [v_out] = adjustVal(v_in)
            if(v_in==0)
                v_out = v_in + eps;
            elseif(v_in==1)
                v_out = v_in - eps;
            else
                v_out = v_in;
            end
        end
    end
    
    methods
        function obj = rvEmpiricalInfo(domain, density, distribution)
            obj.domain = domain;
            obj.density = density;
            obj.distribution = distribution;
        end
        
        function [idx, domainVal] = findClosestDomainVal(obj, q)
            % find the closest value in the domain
            
            % TODO: currently, in discrete RV's it assumes that the query
            %       will be in the domain of the empirical density already
            %       created.  error protection here?
            
            if(q > obj.domain(end))
                domainVal = obj.domain(end);
                idx = length(obj.domain);
            elseif(q < obj.domain(1))
                domainVal = obj.domain(1);
                idx = 1;
            else
                tmp = abs(obj.domain-q);
                [~, idx] = min(tmp); %index of closest value
                domainVal = obj.domain(idx); %closest value
            end
        end
        
        function [val] = queryDensity(obj, q)
            idx = obj.findClosestDomainVal(q);
            val = obj.density(idx);
        end
        
        function [val] = queryDistribution(obj, q)
            idx = obj.findClosestDomainVal(q);
            val = rvEmpiricalInfo.adjustVal(obj.distribution(idx));
        end
        
        function [idx, distributionVal] = findClosestDistributionVal(obj, q)
            if(q > obj.distribution(end))
                distributionVal = obj.distribution(end);
                idx = length(obj.distribution);
            elseif(q < obj.distribution(1))
                distributionVal = obj.distribution(1);
                idx = 1;
            else
                tmp = abs(obj.distribution-q);
                [~, idx] = min(tmp); %index of closest value
                distributionVal = obj.distribution(idx); %closest value
            end
        end
        
        function [val] = invDistribution(obj, u)
            idx = obj.findClosestDistributionVal(u);
            val = obj.domain(idx);
        end
        
    end
    
end

