function [y] = sortSubblocks(x,y)
% x - expected to be sorted already
% y - the var for which subblocks are to be sorted

uDiff = diff(x);
uIdxs = find(uDiff==0);
vBlocks = cell(1,length(uIdxs));
vBlocksIdx = 1;
for ii=1:length(uIdxs)
    if(ii<length(uIdxs))
        if((uIdxs(ii)+1)==uIdxs(ii+1))
            % fast forward
            iStart = ii;
            while(ii<length(uIdxs))
                if((uIdxs(ii)+1)~=uIdxs(ii+1))
                    break;
                end
                ii = ii + 1;
            end
            vBlock = uIdxs(iStart):1:uIdxs(ii)+1;
            vBlocks{vBlocksIdx} = vBlock;
            vBlocksIdx = vBlocksIdx + 1;
            ii = ii + 1;
            if(ii>length(uIdxs))
                break;
            end
        else
            vBlock = [uIdxs(ii) uIdxs(ii)+1];
            vBlocks{vBlocksIdx} = vBlock;
            vBlocksIdx = vBlocksIdx + 1;
        end
    else
        vBlock = [uIdxs(ii) uIdxs(ii)+1];
        vBlocks{vBlocksIdx} = vBlock;
        vBlocksIdx = vBlocksIdx + 1;
    end
end
numVBlocks = vBlocksIdx - 1;
for ii=1:numVBlocks
    vBlock = vBlocks{ii};
    y(vBlock) = sort(y(vBlock));
end

end