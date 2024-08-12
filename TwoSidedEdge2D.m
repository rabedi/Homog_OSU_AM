classdef TwoSidedEdge2D
    properties
        match; % 1 (v1,2 of one to v1,2 of other | 2 v1,2 of one to v2,1 of the other)
        % inside edge info
        insideEdge = Edge2D;
        insideEdge_grain_ID;
        insideEdge_grain_indexInDomain = 0;
        % refers to index in twoSidedEdgeGroups
        insideEdge_groupIndex_InGrain = 0;
        
        % outside edge info (not meaningful for edges on domain boundary)
        outsideEdge = Edge2D;
        outsideEdge_grain_ID;
        outsideEdge_grain_indexInDomain = 0;
        outsideEdge_groupIndex_InGrain = 0;
    end
    methods
    end
end