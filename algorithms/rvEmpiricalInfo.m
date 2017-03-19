classdef rvEmpiricalInfo
    %RVEMPIRICALINFO Defines an empirical distribution function
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
    %* 
    %**************************************************************************
    
    properties
        domain;
        density;
        distribution;
        isdiscrete;
        
        % properties that are calculated internally
        mean;
        variance;
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
        function obj = rvEmpiricalInfo(domain, density, distribution, isdiscrete)
            obj.domain = domain;
            obj.density = density;
            obj.distribution = distribution;
            obj.isdiscrete = isdiscrete;
            
            % calculate the mean and variance
            if(~isempty(density))
                vecCalc = obj.density;
            elseif(~isempty(distribution))
                vecCalc = obj.distribution;
            else
                error('Need to specify atleast a density or distribution!');
            end
            obj.mean = trapz(obj.domain,vecCalc.*obj.domain);
            obj.variance = trapz(obj.domain,vecCalc.*((obj.domain).^2))-(trapz(obj.domain,vecCalc.*obj.domain))^2;
        end
        
        function [idx, domainVal] = findClosestDomainVal(obj, q)
            % find the closest value in the domain
            
            % TODO: currently, in discrete RV's it assumes that the query
            %       will be in the domain of the empirical density already
            %       created.  error protection here?
            
            if(q > obj.domain(end))
                domainVal = obj.domain(end);
                idx = length(obj.domain);
            elseif(~obj.isdiscrete && q < obj.domain(1))
                domainVal = obj.domain(1);
                idx = 1;
            elseif(obj.isdiscrete && q < obj.domain(2))
                domainVal = obj.domain(2);
                idx = 2;
            else
                tmp = abs(obj.domain-q);
                [~, idx] = min(tmp); %index of closest value
                domainVal = obj.domain(idx); %closest value
            end
        end
        
        function [val] = pdf(obj, q)
            idx = obj.findClosestDomainVal(q);
            val = obj.density(idx);
        end
        
        function [val] = cdf(obj, q)
            idx = obj.findClosestDomainVal(q);
            val = rvEmpiricalInfo.adjustVal(obj.distribution(idx));
        end
        
        function [val] = icdf(obj, u)
            idx = obj.findClosestDistributionVal_(u);
            val = obj.domain(idx);
        end
        
        function [idx, distributionVal] = findClosestDistributionVal_(obj, q)
            if(q > obj.distribution(end))
                distributionVal = obj.distribution(end);
                idx = length(obj.distribution);
            elseif(~obj.isdiscrete && q < obj.distribution(1))
                distributionVal = obj.distribution(1);
                idx = 1;
            elseif(obj.isdiscrete && q < obj.distribution(2))
                % we do idx=2 b/c we keep the 0 entry for discrete also
                distributionVal = obj.distribution(2);
                idx = 2;
            else
                tmp = abs(obj.distribution-q);
                [~, idx] = min(tmp); %index of closest value
                distributionVal = obj.distribution(idx); %closest value
            end
        end
        
    end
    
end

