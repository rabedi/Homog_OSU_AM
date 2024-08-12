classdef Edge2D
    properties
        edge_node1 = Node;
        edge_node2 = Node;
        % list of [grain, element (index & id), and edge number]
        % corresponding to this edge
        containingGrain;
        containingElementIsTri;
        containingElementIndex;
        containingElementID;
        containingEdgeNumber;
        
        
        % boolean that helps processing these objects
        processed = 0;
    end
    methods
        % 1 means v1 -> v1 neighbor, 2 -> 2 | 2: 1 -> 2 and 2 -> 1
        function match = DoEdgesMath(obj, otherEdge, tol)
            if (nargin < 3)
                tol = 1e-6;
            end
            match = 0;
            delNode1_other_node1 = norm(obj.edge_node1.crd - otherEdge.edge_node1.crd);
            if (delNode1_other_node1 < tol) % first nodes match
                if (norm(obj.edge_node2.crd - otherEdge.edge_node2.crd) < tol)
                    match = 1;
                end
                return;
            end
            % so node 1's don't match, see if node 1 match node 2 of other
            delNode1_other_node2 = norm(obj.edge_node1.crd - otherEdge.edge_node2.crd);
            if (delNode1_other_node2 < tol) % first nodes match
                if (norm(obj.edge_node2.crd - otherEdge.edge_node1.crd) < tol)
                    match = 2;
                end
                return;
            end
        end
    end
end