classdef ktauhat_s < handle
    %ktauhat_s Definition of the "streaming" ktauhat, a more efficient
    %implementation of ktauhat for computing the RSDM than the naive
    %implementation given in ktauhat.m
    %
    % TODO:
    %  [ ] - All maps could be made faster by preallocation w/ the keys
    %        from a call to unique(x) and unique(y). There may be a way to
    %        get the unique elements from x and y after the sort call,
    %        without calling unique, which attempts re-sorts the array.
    %        WARNING: you would then need to change the call from "isKey"
    %        to see if the value in the map is 0 or greater than 0.
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
    end
    
    properties(SetAccess = private)
        iiBegin;
        iiEnd;
        K;
        uu;
        vv;
        uMap;
        vMap;
        uOvlpMap;
        vOvlpMap;
        mm;
        mm_choose_2;
        
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
        uOvlpMap_prev;
        vOvlpMap_prev;
        mm_prev;
        mm_choose_2_prev;
    end
    
    methods(Access = private)
        function [metric] = consume__(obj)
            obj.iiEnd = obj.iiEnd + 1;
            obj.mm = obj.mm + 1;
            obj.mm_choose_2 = obj.mm_choose_2 + obj.mm - 1;

            uTmp = obj.u(obj.iiBegin:obj.iiEnd);
            vTmp = obj.v(obj.iiBegin:obj.iiEnd);
            uCurrSamp = obj.u(obj.iiEnd); uPrevSamp = obj.u(obj.iiEnd-1);
            vCurrSamp = obj.v(obj.iiEnd); vPrevSamp = obj.v(obj.iiEnd-1);

            uDiff = uTmp(end)-uTmp(end-1:-1:1);
            vDiff = vTmp(end)-vTmp(end-1:-1:1);
            uLastSampDiff = uCurrSamp-uPrevSamp;
            vLastSampDiff = vCurrSamp-vPrevSamp;
            uPosTmp = sum(uDiff>0 & vDiff~=0); vPosTmp = sum(vDiff>0 & uDiff~=0);
            uNegTmp = sum(uDiff<0 & vDiff~=0); vNegTmp = sum(vDiff<0 & uDiff~=0);
            
            uAddVal = vPosTmp-vNegTmp;
            if(uPosTmp<uNegTmp)
                uAddVal = -1*uAddVal;
            end
            obj.K = obj.K+uAddVal;
            
            % count overlaps -- for now, we do this in separate if/else
            % statements to make the code clear, but in C++ this can be
            % optimized
            if(uLastSampDiff==0 && vLastSampDiff~=0)
                if(~isKey(obj.uOvlpMap,uCurrSamp) && ~isKey(obj.uOvlpMap, uPrevSamp))
                    obj.uOvlpMap(uCurrSamp) = 0;
                end
                obj.uOvlpMap(uCurrSamp) = obj.uOvlpMap(uCurrSamp) + 1;
            end
            if(vLastSampDiff==0 && uLastSampDiff~=0)
                if(~isKey(obj.vOvlpMap,vCurrSamp) && ~isKey(obj.vOvlpMap, vPrevSamp))
                    obj.vOvlpMap(vCurrSamp) = 0;
                end
                obj.vOvlpMap(vCurrSamp) = obj.vOvlpMap(vCurrSamp) + 1;
            end

            % TODO: C++ implementation of this must be faster in order to
            % achieve real speed-up w/ RSDM-S rather than RSDM.  We might
            % need to do something clever where we convert the inputs to
            % ranks, and then index into a vector rather than a map
            % container to obtain speed-ups necessary
            obj.uMap(uCurrSamp) = obj.uMap(uCurrSamp) + 1;
            obj.uu = obj.uu + obj.uMap(uCurrSamp) - 1;
            obj.vMap(vCurrSamp) = obj.vMap(vCurrSamp) + 1;
            obj.vv = obj.vv + obj.vMap(vCurrSamp) - 1;

            % attempt to automatically determine if we have hybrid-data or
            % all discrete/continuous
            uuCloseToZero = obj.closeToZero(obj.uu);
            vvCloseToZero = obj.closeToZero(obj.vv);
            
            % hybrid scenario
            if( (uuCloseToZero && obj.vv>0) || (obj.uu>0 && vvCloseToZero) )
                hybridFlag = 1;
                % determine which variable is continuous, and which is
                % discrete
                if(uuCloseToZero)
                    continuousRvIndicator = 0;
                else
                    continuousRvIndicator = 1;
                end
                % compute the correction factor
                cf = obj.correctionFactor(continuousRvIndicator);
                tt = max(obj.uu, obj.vv)-cf;
%                 fprintf('>>> cf=%0.02f, tt=%0.02f\n', cf, tt);
            else
                hybridFlag = 0;
            end
            
            if(hybridFlag)
                den = sqrt(obj.mm_choose_2-tt)*sqrt(obj.mm_choose_2-tt);
            else
                % all continuous/discrete case
                % we assume u is independent var, v is the dependent var
                den = sqrt(obj.mm_choose_2-obj.uu)*sqrt(obj.mm_choose_2-obj.vv);
            end
%             fprintf('>>> K=%d uu=%d vv=%d u(closeToZero)=%d v(closeToZero)=%d\n', ...
%                 obj.K, obj.uu, obj.vv, uuCloseToZero, vvCloseToZero);

            if(obj.K==0 || den==0)
                metric = 0;
            else
                metric = obj.K/den;
            end
        end
        
        
        function [out] = closeToZero(obj, in)
            out = 1;
            mmFloor = floor(obj.mm*obj.CLOSE_TO_ZERO_THRESH);
            
            if(mmFloor>=2)
                % TODO: remove this nchoosek call for high performance!!
                cmpVal = nchoosek(mmFloor,2);
            else
                cmpVal = 0;
            end

            if(in>cmpVal)
                out = 0;
            end
        end
        
        function [cf] = correctionFactor(obj,continuousRvIndicator)
            if(continuousRvIndicator==0)
                valVec = cell2mat(obj.uOvlpMap.values());
            else
                valVec = cell2mat(obj.vOvlpMap.values());
            end
%             continuousRvIndicator
%             valVec
            meanVal = floor(mean(valVec));
            lenMap = length(valVec);
            if(isnan(meanVal) || meanVal<2)
                cf = 0;
            else
                cf = nchoosek(meanVal,2)*lenMap;
            end
        end
    end
    
    methods
        function obj = ktauhat_s(x, y)
            % set immutable properties
            obj.M = length(x);
            obj.CLOSE_TO_ZERO_THRESH = 0.02;      % if we are > 2% of length in terms of combinations;
            [obj.u,I] = sort(x);
            obj.v = y(I);
            
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
            obj.uOvlpMap_prev = obj.uOvlpMap;
            obj.vOvlpMap_prev = obj.vOvlpMap;
            obj.mm_prev = obj.mm;
            obj.mm_choose_2_prev = obj.mm_choose_2;
            
            for ii=1:numPts-1
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
            obj.uOvlpMap = obj.uOvlpMap_prev;
            obj.vOvlpMap = obj.vOvlpMap_prev;
            obj.mm = obj.mm_prev;
            obj.mm_choose_2 = obj.mm_choose_2_prev;
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
            obj.mm_choose_2 = 0;
            
            % reinitialize the maps 
            % WARNING!: overwriting objects here ... is this bad coding practice?
            obj.uMap = containers.Map(obj.u, zeros(1,length(obj.u)));
            obj.vMap = containers.Map(obj.v, zeros(1,length(obj.v)));
            
            remove(obj.uOvlpMap, obj.uOvlpMap.keys());
            remove(obj.vOvlpMap, obj.vOvlpMap.keys());
            
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
            obj.mm_choose_2 = 0;
            
            % initialize the maps
            obj.uMap = containers.Map(obj.u, zeros(1,length(obj.u)));
            obj.vMap = containers.Map(obj.v, zeros(1,length(obj.v)));
            obj.uOvlpMap = containers.Map('KeyType', 'double', 'ValueType', 'uint64');
            obj.vOvlpMap = containers.Map('KeyType', 'double', 'ValueType', 'uint64');
            
            % add the first sample to the maps
            obj.uMap(obj.u(obj.iiBegin)) = 1;
            obj.vMap(obj.v(obj.iiBegin)) = 1;
        end
        
    end
end