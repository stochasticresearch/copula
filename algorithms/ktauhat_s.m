classdef ktauhat_s < handle
    %ktauhat_s Definition of the "streaming" ktauhat, a more efficient
    %implementation of ktauhat for computing the RSDM than the naive
    %implementation given in ktauhat.m
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
    
    properties(SetAccess = immutable)
        M;
        u;
        v;
        CLOSE_TO_ZERO_THRESH;
        ONE_OVER_CLOSE_TO_ZERO_THRESH;
    end
    
    properties(SetAccess = private)
        iiBegin;
        iiEnd;
        K;
        uu;
        vv;
        uMap;
        vMap;
        mm;
        mmGroups;
        mm_choose_2;
        closeToZeroThresh;
        
        % for rewind capability
        % TODO: is there a better (more efficient & cleaner way) to do
        %       rewind?
        iiBegin_prev;
        iiEnd_prev;
        K_prev;
        uu_prev;
        vv_prev;
        uMap_prev;
        vMap_prev;
        mm_prev;
        mmGroups_prev;
        mm_choose_2_prev;
        closeToZeroThresh_prev;
    end
    
    methods(Access = private)
        function [metric] = consume__(obj)
            obj.iiEnd = obj.iiEnd + 1;
            obj.mm = obj.mm + 1;
            obj.mm_choose_2 = obj.mm_choose_2 + obj.mm - 1;

            uTmp = obj.u(obj.iiBegin:obj.iiEnd);
            vTmp = obj.v(obj.iiBegin:obj.iiEnd);
            uCurrSamp = obj.u(obj.iiEnd);
            vCurrSamp = obj.v(obj.iiEnd);

            uDiff = uCurrSamp-uTmp(end-1:-1:1);
            vDiff = vCurrSamp-vTmp(end-1:-1:1);
            uPosTmp = sum(uDiff>0 & vDiff~=0); vPosTmp = sum(vDiff>0 & uDiff~=0);
            uNegTmp = sum(uDiff<0 & vDiff~=0); vNegTmp = sum(vDiff<0 & uDiff~=0);
            
            uAddVal = vPosTmp-vNegTmp;
            if(uPosTmp<uNegTmp)
                uAddVal = -1*uAddVal;
            end
            obj.K = obj.K+uAddVal;
            
            obj.uMap(uCurrSamp) = obj.uMap(uCurrSamp) + 1;
            obj.uu = obj.uu + obj.uMap(uCurrSamp) - 1;
            obj.vMap(vCurrSamp) = obj.vMap(vCurrSamp) + 1;
            obj.vv = obj.vv + obj.vMap(vCurrSamp) - 1;

            % attempt to automatically determine if we have hybrid-data or
            % all discrete/continuous
            % the below block of code is a recursive replacement of the
            % "closeToZero" function in ktauhat
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(~mod(obj.mm, obj.ONE_OVER_CLOSE_TO_ZERO_THRESH))
                obj.mmGroups = obj.mmGroups + 1;
                obj.closeToZeroThresh = obj.closeToZeroThresh + obj.mmGroups - 1;
            end
            uuCloseToZero = obj.uu <= obj.closeToZeroThresh;
            vvCloseToZero = obj.vv <= obj.closeToZeroThresh;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % hybrid scenario
            if( (uuCloseToZero && obj.vv>0) || (obj.uu>0 && vvCloseToZero) )
                tt = max(obj.uu, obj.vv);
                den = sqrt(obj.mm_choose_2-tt)*sqrt(obj.mm_choose_2-tt);
            else
                % all continuous/discrete case
                % we assume u is independent var, v is the dependent var
                den = sqrt(obj.mm_choose_2-obj.uu)*sqrt(obj.mm_choose_2-obj.vv);
            end
            
            if(obj.K==0 || den==0)
                metric = 0;
            else
                metric = obj.K/den;
            end
        end
    end
    
    methods
        function obj = ktauhat_s(x, y)
            % set immutable properties
            obj.M = length(x);
            obj.CLOSE_TO_ZERO_THRESH = 0.02;      % if we are > 2% of length in terms of combinations;
            obj.ONE_OVER_CLOSE_TO_ZERO_THRESH = 1/obj.CLOSE_TO_ZERO_THRESH;
            
            % convert data to an integer based representation
            u = sort(x); [~,uu] = ismember(x,u);
            v = sort(y); [~,vv] = ismember(y,v);
            [obj.u,I] = sort(uu);
            obj.v = vv(I);
            
            % to properly process hybrid and discrete data, we need to sort
            % each "block" of v ... i.e., all the y's with the same u
            obj.v = sortSubblocks(obj.u, obj.v);
            
            % initialize all the required variables for processing
            obj.resetState();
        end
        
        function [metric] = consume(obj, numPts)
            obj.iiBegin_prev = obj.iiBegin;
            obj.iiEnd_prev = obj.iiEnd;
            obj.K_prev = obj.K;
            obj.uu_prev = obj.uu;
            obj.vv_prev = obj.vv;
            obj.uMap_prev = obj.uMap;
            obj.vMap_prev = obj.vMap;
            obj.mm_prev = obj.mm;
            obj.mm_choose_2_prev = obj.mm_choose_2;
            obj.closeToZeroThresh_prev = obj.closeToZeroThresh;
            obj.mmGroups_prev = obj.mmGroups;
            
            for ii=1:numPts
                if(obj.iiEnd<obj.M)
                    metric = obj.consume__();
                else
                    break;
                end
            end
        end
        
        function [metric] = consumeAll(obj)
            metric = obj.consume(obj.M);
        end
        
        function [] = rewind(obj)
            obj.iiBegin = obj.iiBegin_prev;
            obj.iiEnd = obj.iiEnd_prev;
            obj.K = obj.K_prev;
            obj.uu = obj.uu_prev;
            obj.vv = obj.vv_prev;
            obj.uMap = obj.uMap_prev;
            obj.vMap = obj.vMap_prev;
            obj.mm = obj.mm_prev;
            obj.mm_choose_2 = obj.mm_choose_2_prev;
            obj.closeToZeroThresh = obj.closeToZeroThresh_prev;
            obj.mmGroups = obj.mmGroups_prev;
        end
        
        function [] = clearState(obj)
            % called when we want to start a new rectangle
            
            % reset all private properties defined above
            obj.iiBegin = obj.iiEnd + 1;
            obj.iiEnd = obj.iiBegin;
            obj.uu = 0;
            obj.vv = 0;
            obj.K = 0;
            obj.mm = 1;
            obj.mmGroups = 0;
            obj.mm_choose_2 = 0;
            obj.closeToZeroThresh = 0;
            
            % reinitialize the maps 
            obj.uMap = zeros(1,length(obj.u));
            obj.vMap = zeros(1,length(obj.v));
            
            % initialize the uMap and vMap by add the first element into it
            obj.uMap(obj.u(obj.iiBegin)) = 1;
            obj.vMap(obj.v(obj.iiBegin)) = 1;
        end
        
        function [] = resetState(obj)
            % called when we want to "reset".  We do this when we
            % re-process the data at a different scan rate
            
            % initialize remaining properties
            obj.iiBegin = 1;
            obj.iiEnd = 1;
            obj.K = 0;
            obj.uu = 0;
            obj.vv = 0;
            obj.mm = 1;     % we start w/ 1 sample, and consume the next, etc...
            obj.mmGroups = 0;
            obj.mm_choose_2 = 0;
            obj.closeToZeroThresh = 0;
            
            % initialize the maps
            % we now just use arrays b/c we converted the data to be
            % integer representation only
            obj.uMap = zeros(1,length(obj.u));
            obj.vMap = zeros(1,length(obj.v));
            
            % add the first sample to the maps
            obj.uMap(obj.u(obj.iiBegin)) = 1;
            obj.vMap(obj.v(obj.iiBegin)) = 1;
        end
        
    end
end