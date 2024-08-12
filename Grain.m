classdef Grain
    properties
        id = 1;
        
        nodeNum = 0;
        nodeIds; % vector of numNodes size
        nodeCrds; % numNodes x 2 matrix
        grainxm, grainxM, grainym, grainyM; % grain bounding box
        
        % for node location i, nodeIncidentElementIDs(i) contains the
        % IDs of elements that are incident to it
        nodeIncident_tri_ElementIDs;
        nodeIncident_quad_ElementIDs;
        
        triElements_num = 0;
        triElementsConnectivity; % [num x 3] matrix of nodal connectivty
        triElementsIDs;
        
        quadElements_num = 0;
        quadElementsConnectivity; % [num x 4] matrix of nodal connectivty
        quadElementsIDs;
        
        gb_elements; % gb_elements{I}{J} is the surface boundary definition I, part J,
        % for example for grain 1: I = 1 corresponds to boundary of grain 1 to grain 2 (_SURF-MASTER_1_2)
        % and there are two segments for it (SURF-MASTER_1_2_S1 and
        % SURF-MASTER_1_2_S2) so J goes from 1 to 2
        
        num_gb_neighbors = 0; % how many grains are neighbor to this grain % 3 because it's neighbor to 2, 9, 17
        gb_neighbor_ids = []; % in this example this will be [2, 9, 17]
        num_segments_p_neighbors = cell(0); % indexed "similar" to gp_elements {I}(J)
        % this is the number of segments for each neighborhood. For first neighbor 2, there are 2 segments
        % SURF-MASTER_1_2_S1, SURF-MASTER_1_2_S2, for 9 and 17 there are 3, so this vector would be [2, 3, 3]
        continuousNodeNumbering = 1;

        % to read node on sides of a VE
        bottom_side = 1;
        right_side = 2;
        top_side = 3;
        left_side = 4;
        side_names = {'BOTTOM', 'RIGHT', 'TOP', 'LEFT'};
        % index by left, right, ..., contains nodes there
        side_nodes = cell(0);
        % left mode node can also be saved
        left_mode_node = [];
        
        
        %%%%%%%%%% edge neighborhood information
        twoSidedEdgeGroups = cell(0); %SetOfTwoSidedEdges2D;
        domainBoundaryNodeIDs;
        domainBoundaryNodeCrds;
    end
    methods
        function objout = read(obj, fid, idIn)
            objout = obj;
            objout.id = idIn;
            buf = fscanf(fid, '%s', 1);
            while (strcmp(buf, '*End') == 0)
                sz = length(buf);
                if (sz == 0)
                    continue;
                end
                hasStar = (buf(1) == '*');
                if (hasStar && (sz > 1) && (buf(2) == '*')) % this is a comment ...
                    wholeCommentLine = fgetl(fid);
                    buf = fscanf(fid, '%s', 1);
                end
                hasStar = (buf(1) == '*');
                if (~hasStar)
                    buf
                    fprintf(1, 'a bug in reading ... as each block should start with *\n');
                    pause;
                end
                if (strcmp(buf, '*Node') == 1) % reading the node block
                    objout.nodeNum = 0;
                    buf2 = fscanf(fid, '%s', 1);
                    while (strcmp(buf2, '**') == 0)
                        objout.nodeNum = objout.nodeNum + 1;
                        ln = length(buf2);
                        idVal = str2num(buf2(1:ln - 1));
                        objout.nodeIds(objout.nodeNum) = idVal;
                        
                        x = fscanf(fid, '%f', 1);
                        buf2 = fscanf(fid, '%s', 1); % ,
                        y = fscanf(fid, '%f', 1);
                        buf2 = fscanf(fid, '%s', 1); % ,
                        
                        objout.nodeCrds(objout.nodeNum, 1) = x;
                        objout.nodeCrds(objout.nodeNum, 2) = y;
                        
                        % beginning of the next line to check if it's an ID
                        % or end of block
                        buf2 = fscanf(fid, '%s', 1);
                    end
                    objout.grainxm = min(objout.nodeCrds(:, 1));
                    objout.grainxM = max(objout.nodeCrds(:, 1));
                    objout.grainym = min(objout.nodeCrds(:, 2));
                    objout.grainyM = max(objout.nodeCrds(:, 2));
                    buf = fscanf(fid, '%s', 1);
                elseif (strcmp(buf ,'*ELEMENT,') == 1) % reading the elements
                    buf2 = fscanf(fid, '%s', 1);
                    if (strcmp(buf2, 'TYPE=CPS3,') == 1) % triangles
                        nnPerEl = 3;
                    else
                        nnPerEl = 4;
                    end
                    remainderOfTheLine = fgetl(fid);
                    elConnectivity = [];
                    
                    numEl = 0;
                    buf2 = fscanf(fid, '%s', 1); % beginning of the next line
                    while (buf2(1) ~= '*')
                        numEl = numEl + 1;
                        ln = length(buf2);
                        idVal = str2num(buf2(1:ln - 1));
                        elIDs(numEl) = idVal;
                        
                        for j = 1:nnPerEl
                            elConnectivity(numEl, j) = fscanf(fid, '%f', 1);
                            if (j ~= nnPerEl)
                                comma = fscanf(fid, '%s', 1);
                            end
                        end
                        buf2 = fscanf(fid, '%s', 1); % beginning of the next line
                    end
                    if (nnPerEl == 3)
                        objout.triElements_num = numEl;
                        objout.triElementsConnectivity = elConnectivity;
                        objout.triElementsIDs = elIDs;
                    else
                        objout.quadElements_num = numEl;
                        objout.quadElementsConnectivity = elConnectivity;
                        objout.quadElementsIDs = elIDs;
                    end
                    buf = buf2;
                elseif (strcmp(buf, '*Surface,') == 1)
                    remainderOfTheLine = fgetl(fid);
                    I = objout.num_gb_neighbors;
                    numSeg = sum(objout.num_segments_p_neighbors{I} > 0);
                    for i = 1:numSeg
                        remainderOfTheLine = fgetl(fid);
                    end
                    buf = fscanf(fid, '%s', 1);
                elseif (strcmp(buf, '*Elset,') == 1)
                    a = fscanf(fid, '%s', 1);
                    b = cell(0);
                    for i = 1:5
                        [b{i}, a] = strtok(a, '_');
                    end
                    otherGrain = str2num(b{4});
                    partNo = b{5};
                    partNo = str2num(partNo(2:length(partNo)));
                    remainderOfTheLine = fgetl(fid);
                    buf2 = fscanf(fid, '%s', 1);
                    numEl = 0;
                    els = [];
                    while ((strcmp(buf2, '*Elset,') == 0) && (strcmp(buf2, '*Surface,') == 0))
                        numEl = numEl + 1;
                        els(numEl) = str2num(buf2);
                        buf2 = fscanf(fid, '%s', 1);
                    end
                    I = length(objout.gb_elements);
                    J = partNo;
                    if (length(find(objout.gb_neighbor_ids == otherGrain)) == 0)
                        I = I + 1;
                        objout.num_gb_neighbors = objout.num_gb_neighbors + 1;
                        objout.gb_neighbor_ids(I) = otherGrain;
                        objout.num_segments_p_neighbors{I} = zeros(3, 1);
                        objout.gb_elements{I} = cell(3, 1);
                    end
                    objout.gb_elements{I}{J} = els;
                    objout.num_segments_p_neighbors{I}(J) = length(els);
                    buf = buf2;
                elseif (strcmp(buf, '*SOLID') == 1)
                    remainderOfTheLine = fgetl(fid);
                    remainderOfTheLine = fgetl(fid);
                    remainderOfTheLine = fgetl(fid);
                    buf = fscanf(fid, '%s', 1);
                elseif (strcmp(buf, '*Nset,') == 1)
                    buf = fscanf(fid, '%s', 1);
                    if (contains(buf, 'LEFT_MOST'))
                        objout.left_mode_node = fscanf(fid, '%d');
                        buf = fscanf(fid, '%s', 1);
                    else
                        if (contains(buf, 'BOTTOM'))
                            side = objout.bottom_side;
                        elseif (contains(buf, 'TOP'))
                            side = objout.top_side;
                        elseif (contains(buf, 'LEFT'))
                            side = objout.left_side;
                        elseif (contains(buf, 'RIGHT'))
                            side = objout.right_side;
                        else
                            fprintf(1, 'Invalid string: %s\n', buf);
                            pause;
                        end
                        buf = fscanf(fid, '%s', 1);
                        nid = str2num(buf);
                        cnt = 0;
                        vc = [];
                        while (length(nid) > 0)
                            cnt = cnt + 1;
                            vc(cnt) = nid;
                            buf = fscanf(fid, '%s', 1);
                            nid = str2num(buf);
                        end
                        objout.side_nodes{side} = vc;
                    end
%                    buf = fscanf(fid, '%s', 1);
                end
            end
            buf = fscanf(fid, '%s', 1);
        end
        function write(obj, fid)
            fprintf(fid, '*Node\n');
            for i = 1:obj.nodeNum
                fprintf(fid, '\t\t%d,\t\t%f,\t\t%f,\n', obj.nodeIds(i), obj.nodeCrds(i, 1), obj.nodeCrds(i, 2));
            end
            fprintf(fid, '**\n');
            if (obj.triElements_num > 0)
                fprintf(fid, '*ELEMENT, TYPE=CPS3, ELSET=MATERIAL_%dELSET\n', obj.id);
                for i = 1:obj.triElements_num
                    fprintf(fid, '\t\t%d,\t\t%d,\t\t%d,\t\t%d\n', obj.triElementsIDs(i), obj.triElementsConnectivity(i, 1), obj.triElementsConnectivity(i, 2), obj.triElementsConnectivity(i, 3));
                end
                fprintf(fid, '**\n');
            end
            if (obj.quadElements_num > 0)
                fprintf(fid, '*ELEMENT, TYPE=CPS4, ELSET=MATERIAL_%dELSET\n', obj.id);
                for i = 1:obj.quadElements_num
                    fprintf(fid, '\t\t%d,\t\t%d,\t\t%d,\t\t%d,\t\t%d\n', obj.quadElementsIDs(i), obj.quadElementsConnectivity(i, 1), obj.quadElementsConnectivity(i, 2), obj.quadElementsConnectivity(i, 3), obj.quadElementsConnectivity(i, 4));
                end
                fprintf(fid, '**\n');
            end
            
            for i = 1:length(obj.side_nodes)
                nodeSet = obj.side_nodes{i};
                szSet = length(nodeSet);
                if (szSet > 0)
                    fprintf(fid, '*Nset, nset=%s_NODES_%d\n',obj.side_names{i}, obj.id);
                    for j = 1:szSet
                        fprintf(fid, '%d,\t', nodeSet(j));
                        if ((mod(j, 16) == 0) && (j ~= szSet))
                            fprintf(fid, '\n');
                        end
                    end
                    fprintf(fid, '\n**\n');
                end
            end
            sz_boundary_nodes = length(obj.domainBoundaryNodeIDs);
            for i = 1:sz_boundary_nodes
                nodeID = obj.domainBoundaryNodeIDs(i);
                fprintf(fid, '*Nset, nset=BOUNDARY_NODES_%d_%d\n', i, obj.id);
                fprintf(fid, '%d\n', nodeID);
            end
            if (length(obj.left_mode_node) > 0)
                    fprintf(fid, '*Nset, nset=LEFT_MOST_NODE_%d\n', obj.id);
                    fprintf(fid, '%d\n', obj.left_mode_node);
            end
            for I = 1:obj.num_gb_neighbors
                neighborGrainID = obj.gb_neighbor_ids(I);
                if (neighborGrainID <= 0)
                    continue;
                end
                num_segments_p_neighbor = obj.num_segments_p_neighbors{I};
                neighborGrainID = obj.gb_neighbor_ids(I);
                sm = sum(num_segments_p_neighbor);
                if (sm == 0)
                    continue;
                end
                if (obj.id < neighborGrainID)
                    word = ['_SURF-MASTER_'];
                else
                    word = ['_SURF-SLAVE_'];
                end
                word = [word, num2str(obj.id), '_', num2str(neighborGrainID)];
                % example word: _SURF-MASTER_1_17
                segName = cell(0);
                for J = 1:3
                    if (num_segments_p_neighbor(J) == 0)
                        continue;
                    end
                    segName{J} = [word, '_S', num2str(J)];
                    fprintf(fid, '*Elset, elset=%s, internal\n', segName{J});
                    els = obj.gb_elements{I}{J};
                    for k = 1:length(els)
                        fprintf(fid, '\t\t%d,', els(k));
                        if (mod(k, 16) == 0)
                            fprintf(fid, '\n');
                        end
                    end
                    if (mod(k, 16) ~= 0)
                        fprintf(fid, '\n');
                    end
                end
                wrd = word(2:length(word));
                fprintf(fid, '*Surface, type=ELEMENT, name=%s\n', wrd);
                for J = 1:3
                    if (num_segments_p_neighbor(J) == 0)
                        continue;
                    end
                    fprintf(fid, '%s, S%d\n', segName{J}, J);
                end
            end
            fprintf(fid, '*SOLID SECTION, ELSET=MATERIAL_%dELSET, MATERIAL=MATERIAL_%d\n', obj.id, obj.id);
            fprintf(fid, '1.,\n');
            fprintf(fid, '**\n');
            fprintf(fid, '*End Part\n');
        end
        function fillMesh(obj, clr_tri, clr_quad)
            if nargin < 2
                clr_tri = 'b';
            end
            if nargin < 3
                clr_quad = 'r';
            end
           xs = zeros(3, obj.triElements_num);
           ys = zeros(3, obj.triElements_num);
           for i = 1:obj.triElements_num
               for v = 1:3
                   vid = obj.triElementsConnectivity(i, v);
                   if (obj.continuousNodeNumbering)
                       vloc = vid;
                   else
                       vloc = find(obj.nodeIds == vid);
                   end
                   xs(v, i) = obj.nodeCrds(vloc, 1);
                   ys(v, i) = obj.nodeCrds(vloc, 2);
               end
           end
           fill(xs, ys, clr_tri);
           hold on;
           
           xs = zeros(4, obj.quadElements_num);
           ys = zeros(4, obj.quadElements_num);
           for i = 1:obj.quadElements_num
               for v = 1:4
                   vid = obj.quadElementsConnectivity(i, v);
                   if (obj.continuousNodeNumbering)
                       vloc = vid;
                   else
                       vloc = find(obj.nodeIds == vid);
                   end
                   xs(v, i) = obj.nodeCrds(vloc, 1);
                   ys(v, i) = obj.nodeCrds(vloc, 2);
               end
           end
           fill(xs, ys, clr_quad);
        end
        
        function plotMesh(obj)
           xs = [];
           ys = [];
           cntr = 0;
           for i = 1:obj.triElements_num
               for v = 1:3
                   vid = obj.triElementsConnectivity(i, v);
                   vloc = find(obj.nodeIds == vid);
                   vx(v) = obj.nodeCrds(vloc, 1);
                   vy(v) = obj.nodeCrds(vloc, 2);
               end
               cntr = cntr + 1; xs(cntr) = vx(1);    ys(cntr) = vy(1);  
               cntr = cntr + 1; xs(cntr) = vx(2);    ys(cntr) = vy(2);  
               cntr = cntr + 1; xs(cntr) = vx(3);    ys(cntr) = vy(3);  
               cntr = cntr + 1; xs(cntr) = vx(1);    ys(cntr) = vy(1);  
               cntr = cntr + 1; xs(cntr) = nan;    ys(cntr) = nan;
           end
           
           for i = 1:obj.quadElements_num
               for v = 1:4
                   vid = obj.quadElementsConnectivity(i, v);
                   vloc = find(obj.nodeIds == vid);
                   vx(v) = obj.nodeCrds(vloc, 1);
                   vy(v) = obj.nodeCrds(vloc, 2);
               end
               cntr = cntr + 1; xs(cntr) = vx(1);    ys(cntr) = vy(1);  
               cntr = cntr + 1; xs(cntr) = vx(2);    ys(cntr) = vy(2);  
               cntr = cntr + 1; xs(cntr) = vx(3);    ys(cntr) = vy(3);  
               cntr = cntr + 1; xs(cntr) = vx(4);    ys(cntr) = vy(4);  
               cntr = cntr + 1; xs(cntr) = vx(1);    ys(cntr) = vy(1);  
               cntr = cntr + 1; xs(cntr) = nan;    ys(cntr) = nan;
           end
           plot(xs, ys);
          
%             % Gaige
%             figure(1)            
%             %dealing with tri elements...
%             vs_triElements = obj.nodeCrds;
%             triElements_X_Cords = vs_triElements(:,1);
%             triElements_Y_Cords = vs_triElements(:,2);                        
%             triplot(obj.triElementsConnectivity,triElements_X_Cords,triElements_Y_Cords);
%             hold on
%             %dealing with quad elements...
%             vs_quadElements = obj.nodeCrds;
%             QuadElements_X_Cords = vs_quadElements(:,1);
%             QuadElements_Y_Cords = vs_quadElements(:,2);                                                                    
%             test22 = quadplot(obj.quadElementsConnectivity,QuadElements_X_Cords,QuadElements_Y_Cords); %testing ... remember to change obj.quad back to mappedQuad same for above with tri...
%             hold off
        end
        
        
        % it includes any element (tri or quad) whose centroid IS ( xm < <
        % xM), ym, yM. Once valid elements are known, put those in objout
        % any node that is not used will be written
        function [objout, hasIntersection] = Extract(obj, xmin, xmax, ymin, ymax, renumbering)
            delx = xmax - xmin;
            dely = ymax - ymin;
            tol = 1e-5 * min(delx, dely);
            xm = xmin - tol;
            xM = xmax + tol;
            ym = ymin - tol;
            yM = ymax + tol;
            % use these limits for your bounding box.
            objout = Grain;

   %                 grainxm, , grainym, grainyM; % grain bounding box

            if ((obj.grainxM < xm) || (obj.grainxm > xM) || ...
                    (obj.grainyM < ym) || (obj.grainym > yM))
                hasIntersection = 0;
                return;
            end
            
            objout.id = obj.id;
            
            
            elementLoc = 0;
            % making a temporary storage for the element IDs that are
            % connected to the nodes of the original mesh
            nodeIncident_tri_ElementIDs_Temp = cell(obj.nodeNum, 1);
            triElementLocations2Keep = [];
            triElementIDs2Keep = [];
            for i = 1:obj.triElements_num
                elementID = obj.triElementsIDs(i);
                vids = obj.triElementsConnectivity(i, :);
                cent = [0 0];
                node_locs = zeros(3, 1);
                for j = 1:3
                    [found, node_loc] = find(obj.nodeIds == vids(j));
                    node_locs(j) = node_loc;
                    if (found == 0)
                        fprintf(1, 'node %d not found', vids(j));
                    end
                    v_x = obj.nodeCrds(node_loc, 1);
                    v_y = obj.nodeCrds(node_loc, 2);
                    cent(1) = cent(1) + v_x;
                    cent(2) = cent(2) + v_y;
                end
                cent = cent / 3;
                % now check if the centroid is in the box
                x_cent = cent(1);
                if ((xm > x_cent) || (xM < x_cent))
                    continue;
                end
                y_cent = cent(2);
                if ((ym > y_cent) || (yM < y_cent))
                    continue;
                end
                elementLoc = elementLoc + 1;
                triElementLocations2Keep(elementLoc) = i;
                triElementIDs2Keep(elementLoc) = elementID;
                
                % now update the 3 node connections to this element
                for j = 1:3
                    node_loc = node_locs(j);
                    nodeIncident_tri_ElementIDs_Temp{node_loc} = [nodeIncident_tri_ElementIDs_Temp{node_loc}, elementID];
                end
            end
            
            % quad elements
            elementLoc = 0;
            % making a temporary storage for the element IDs that are
            % connected to the nodes of the original mesh
            nodeIncident_quad_ElementIDs_Temp = cell(obj.nodeNum, 1);
            quadElementLocations2Keep = [];
            quadElementIDs2Keep = [];
            for i = 1:obj.quadElements_num
                elementID = obj.quadElementsIDs(i);
                vids = obj.quadElementsConnectivity(i, :);
                cent = [0 0];
                node_locs = zeros(4, 1);
                for j = 1:4
                    [found, node_loc] = find(obj.nodeIds == vids(j));
                    node_locs(j) = node_loc;
                    if (found == 0)
                        fprintf(1, 'node %d not found', vids(j));
                    end
                    v_x = obj.nodeCrds(node_loc, 1);
                    v_y = obj.nodeCrds(node_loc, 2);
                    cent(1) = cent(1) + v_x;
                    cent(2) = cent(2) + v_y;
                end
                cent = cent / 4;
                % now check if the centroid is in the box
                x_cent = cent(1);
                if ((xm > x_cent) || (xM < x_cent))
                    continue;
                end
                y_cent = cent(2);
                if ((ym > y_cent) || (yM < y_cent))
                    continue;
                end
                elementLoc = elementLoc + 1;
                quadElementLocations2Keep(elementLoc) = i;
                quadElementIDs2Keep(elementLoc) = elementID;
                
                % now update the 4 node connections to this element
                for j = 1:4
                    node_loc = node_locs(j);
                    nodeIncident_quad_ElementIDs_Temp{node_loc} = [nodeIncident_quad_ElementIDs_Temp{node_loc}, elementID];
                end
            end
            objout.triElements_num = length(triElementLocations2Keep);
            objout.quadElements_num = length(quadElementLocations2Keep);
            hasIntersection = ((objout.triElements_num > 0) || (objout.quadElements_num > 0));
            if (hasIntersection == 0)
                return;
            end
            % forming nodes and elements
            objout.nodeNum = 0;
            objout.grainxm = inf;
            objout.grainxM = -inf;
            objout.grainym = inf;
            objout.grainyM = -inf;
            
            for i = 1:obj.nodeNum
                if ((length(nodeIncident_tri_ElementIDs_Temp{i}) == 0) && ...
                        (length(nodeIncident_quad_ElementIDs_Temp{i}) == 0))
                    continue;
                end
                x = obj.nodeCrds(i, 1);
                y = obj.nodeCrds(i, 2);
                
                % acceptable node for the new mesh. Add node count for it
                objout.nodeNum = objout.nodeNum + 1;
                objout.nodeIds(objout.nodeNum) = obj.nodeIds(i);
                objout.nodeCrds(objout.nodeNum, 1) = x;
                objout.nodeCrds(objout.nodeNum, 2) = y;
                
                objout.grainxm = min(objout.grainxm, x);
                objout.grainxM = max(objout.grainxM, x);
                objout.grainym = min(objout.grainym, y);
                objout.grainyM = max(objout.grainyM, y);
            end
            % node is done
            
            % taking care of tri elements
            objout.triElementsConnectivity = zeros(objout.triElements_num, 3);
            objout.triElementsIDs = triElementIDs2Keep;
                        
            for i = 1:objout.triElements_num
                locInObj = triElementLocations2Keep(i);
                for j = 1:3
                    objout.triElementsConnectivity(i, j) = obj.triElementsConnectivity(locInObj, j);
                end
            end
            
            % taking care of quad elements
            objout.quadElementsConnectivity = zeros(objout.quadElements_num, 3);
            objout.quadElementsIDs = quadElementIDs2Keep;
                        
            for i = 1:objout.quadElements_num
                locInObj = quadElementLocations2Keep(i);
                for j = 1:4
                    objout.quadElementsConnectivity(i, j) = obj.quadElementsConnectivity(locInObj, j);
                end
            end
            
            %%%%% forming the neighborhood
            elementIDsRemainingInside = [objout.triElementsIDs, objout.quadElementsIDs];
            objout.num_gb_neighbors = obj.num_gb_neighbors;
            for I = 1:obj.num_gb_neighbors
                gb_element = obj.gb_elements{I};
                objout.gb_neighbor_ids(I) = obj.gb_neighbor_ids(I);
                num_segments_p_neighbor = obj.num_segments_p_neighbors{I};

                for J = 1:3
                    gb_element_4OneS = gb_element{J};
                    num_segments_p_neighbor_4OneS = num_segments_p_neighbor(J);
                    if (num_segments_p_neighbor_4OneS == 0)
                        objout.gb_elements{I}{J} = [];
                        objout.num_segments_p_neighbors{I}(J) = 0;
                        continue;
                    end
                    booleanSet = ismember(gb_element_4OneS, elementIDsRemainingInside);
                    cntr4RemainingElements = 0;
                    gb_element_4OneS_AfterIntersection = [];
                    for k = 1:num_segments_p_neighbor_4OneS
                        if (booleanSet(k) == 0)
                            continue;
                        end
                        % the element remains in the intersection
                        cntr4RemainingElements = cntr4RemainingElements + 1;
                        gb_element_4OneS_AfterIntersection(cntr4RemainingElements) = gb_element_4OneS(k);
                    end
                    objout.gb_elements{I}{J} = gb_element_4OneS_AfterIntersection;
                    objout.num_segments_p_neighbors{I}(J) = length(objout.gb_elements{I}{J});
                end
            end
            objout.continuousNodeNumbering = 0;

            %%%%%%%%%%%% renumbering
            if (~renumbering)
                return;
            end
            objout.continuousNodeNumbering = 1;
            vIDMax = max(objout.nodeIds);
            current_node_ID_to_new = zeros(vIDMax, 1);
            for i = 1:objout.nodeNum
                loc = objout.nodeIds(i);
                current_node_ID_to_new(loc) = i;
            end
            eIDMax = 0;
            if (objout.triElements_num > 0)
                eIDMax = max(objout.triElementsIDs);
            end
            if (objout.quadElements_num > 0)
                eIDMax = max(eIDMax, max(objout.quadElementsIDs));
            end
            current_element_ID_to_new = zeros(eIDMax, 1);
            for i = 1:objout.triElements_num
                loc = objout.triElementsIDs(i);
                current_element_ID_to_new(loc) = i;
            end
            for i = 1:objout.quadElements_num
                loc = objout.quadElementsIDs(i);
                current_element_ID_to_new(loc) = i + objout.triElements_num;
            end
            
            % renumbering nodes
            for i = 1:objout.nodeNum
                objout.nodeIds(i) = i;
            end
            
            % renumbering elements and nodes in them
            for i = 1:objout.triElements_num
                objout.triElementsIDs(i) = i;
                for j = 1:3
                    vidOld = objout.triElementsConnectivity(i, j);
                    vidNew = current_node_ID_to_new(vidOld);
                    objout.triElementsConnectivity(i, j) = vidNew;                 
                end
            end
            for i = 1:objout.quadElements_num
                objout.quadElementsIDs(i) = i + objout.triElements_num;
                for j = 1:4
                    vidOld = objout.quadElementsConnectivity(i, j);
                    vidNew = current_node_ID_to_new(vidOld);
                    objout.quadElementsConnectivity(i, j) = vidNew;                 
                end
            end
            
            % now renumbering elements in boundary element sets
            for I = 1:objout.num_gb_neighbors
                for i = 1:length(objout.gb_elements{I})
                    for k = 1:length(objout.gb_elements{I}{i})
                        eIDOld = objout.gb_elements{I}{i}(k);
                        eIDNew = current_element_ID_to_new(eIDOld);
                        objout.gb_elements{I}{i}(k) = eIDNew;
                    end
                end
            end
        end
        function current_node_ID_to_new = Get_current_node_ID_to_new(obj)
            if (obj.continuousNodeNumbering == 1)
                current_node_ID_to_new =1:length(obj.nodeIds);
                return;
            end
            vIDMax = max(obj.nodeIds);
            current_node_ID_to_new = zeros(vIDMax, 1);
            for i = 1:obj.nodeNum
                loc = obj.nodeIds(i);
                current_node_ID_to_new(loc) = i;
            end
        end
        function listDomainEdges = getDomainEdges(obj)
            listDomainEdges = cell(0);
            current_node_ID_to_new = obj.Get_current_node_ID_to_new();
            numNodes = length(obj.nodeIds);
            edges_bet_nInd1_nInd2 = cell(numNodes, 1);
            edges_bet_nInd1_nInd2_other_node_id = cell(numNodes, 1);
            % going over tri elements
            for ei = 1:obj.triElements_num
                for edgeNo = 1:3
                    vi1 = edgeNo;
                    vi2 = edgeNo + 1;
                    if (vi2 == 4)
                        vi2 = 1;
                    end
                    vloc1 = current_node_ID_to_new(obj.triElementsConnectivity(ei, vi1));
                    vloc2 = current_node_ID_to_new(obj.triElementsConnectivity(ei, vi2));
                    if (vloc2 < vloc1)
                        tmp = vloc1;
                        vloc1 = vloc2;
                        vloc2 = tmp;
                    end
                    edg = Edge2D;
                    edg.containingElementIsTri = 1;
                    edg.containingElementIndex = ei;
                    edg.containingEdgeNumber = edgeNo;
                    
                    [suc, loc] = find(edges_bet_nInd1_nInd2_other_node_id{vloc1} == vloc2);
                    if (isempty(suc))
                        loc = length(edges_bet_nInd1_nInd2_other_node_id{vloc1}) + 1;
                        edges_bet_nInd1_nInd2{vloc1}{loc} = cell(0);
                    end
                    edges_bet_nInd1_nInd2_other_node_id{vloc1}(loc) = vloc2;
                    szTmp = length(edges_bet_nInd1_nInd2{vloc1}{loc});
                    edges_bet_nInd1_nInd2{vloc1}{loc}{szTmp + 1} = edg;
                end
            end
                % quad ones
            for ei = 1:obj.quadElements_num
                for edgeNo = 1:4
                    vi1 = edgeNo;
                    vi2 = edgeNo + 1;
                    if (vi2 == 5)
                        vi2 = 1;
                    end
                    vloc1 = current_node_ID_to_new(obj.quadElementsConnectivity(ei, vi1));
                    vloc2 = current_node_ID_to_new(obj.quadElementsConnectivity(ei, vi2));
                    if (vloc2 < vloc1)
                        tmp = vloc1;
                        vloc1 = vloc2;
                        vloc2 = tmp;
                    end
                    edg = Edge2D;
                    edg.containingElementIsTri = 0;
                    edg.containingElementIndex = ei;
                    edg.containingEdgeNumber = edgeNo;
                    
                    [suc, loc] = find(edges_bet_nInd1_nInd2_other_node_id{vloc1} == vloc2);
                    if (isempty(suc))
                        loc = length(edges_bet_nInd1_nInd2_other_node_id{vloc1}) + 1;
                        edges_bet_nInd1_nInd2{vloc1}{loc} = cell(0);
                    end
                    edges_bet_nInd1_nInd2_other_node_id{vloc1}(loc) = vloc2;
                    szTmp = length(edges_bet_nInd1_nInd2{vloc1}{loc});
                    edges_bet_nInd1_nInd2{vloc1}{loc}{szTmp + 1} = edg;
                end
            end
            % now see which edges have only one attached element
            cntr = 0;
            for vloc1 = 1:numNodes
                connectedNodes = edges_bet_nInd1_nInd2_other_node_id{vloc1};
                num_connectedNodes = length(connectedNodes);
                for j = 1:num_connectedNodes
                    vloc2 = connectedNodes(j);
                    set_connected_edges = edges_bet_nInd1_nInd2{vloc1}{j};
                    sz = length(set_connected_edges);
                    if (sz ~= 1)
                        continue;
                    end
                    edg = set_connected_edges{1};
                    vloc = vloc1;
                    edg.edge_node1.nodeGrainID = obj.id;
                    edg.edge_node1.n_index = vloc;
                    edg.edge_node1.n_id = obj.nodeIds(vloc);
                    edg.edge_node1.crd = obj.nodeCrds(vloc, :);

                    vloc = vloc2;
                    edg.edge_node2.nodeGrainID = obj.id;
                    edg.edge_node2.n_index = vloc;
                    edg.edge_node2.n_id = obj.nodeIds(vloc);
                    edg.edge_node2.crd = obj.nodeCrds(vloc, :);
                    
                    isTri = edg.containingElementIsTri(1);
                    elIndex = edg.containingElementIndex;
                    if (isTri)
                        edg.containingElementID = obj.triElementsIDs(elIndex);
                    else
                        edg.containingElementID = obj.quadElementsIDs(elIndex);
                    end
                    edg.containingGrain = obj.id;
                    
                    cntr = cntr + 1;
                    listDomainEdges{cntr} = edg;
                end
            end
        end
        function plotMeshRestricted(obj) % seems the connectivities are too high for the number of elements
            %so this attemps to map the elements by renumbering the
            %connectivities of both tri and quad elements.            
            figure(2)
            vs_triElements = obj.nodeCrds;
            triElements_X_Cords = vs_triElements(:,1);
            triElements_Y_Cords = vs_triElements(:,2);
            uniq_TriElementConnectivityVals = sort(unique(obj.triElementsConnectivity));
            new_TriElementConnectivityVals = [1:length(uniq_TriElementConnectivityVals)]';
            M = containers.Map(uniq_TriElementConnectivityVals,new_TriElementConnectivityVals);
            triElementsConnectivityDimen = length(obj.triElementsConnectivity);
            triConnectivityRow = triElementsConnectivityDimen(1);
            triConnectivityClmns = triElementsConnectivityDimen(2);
            mappedTriElementConnectivity = [];
            for kr = 1:triConnectivityRow
                for kc = 1:triConnectivityClmns
                    mappedTriElementConnectivity(kr,kc) = M(obj.triElementsConnectivity(kr,kc));
                end
            end 
            triplot(mappedTriElementConnectivity,triElements_X_Cords,triElements_Y_Cords)
            hold on
            %plotting quad elements...
            vs_quadElements = obj.nodeCrds;
            QuadElements_X_Cords = vs_quadElements(:,1);
            QuadElements_Y_Cords = vs_quadElements(:,2);
            uniq_quadElementConnectivityVals = sort(unique(obj.quadElementsConnectivity));
            new_quadElementConnectivityVals = [1:length(uniq_quadElementConnectivityVals)]';
            M1 = containers.Map(uniq_quadElementConnectivityVals,new_quadElementConnectivityVals);
            quadElementsConnectivityDimen = length(obj.quadElementsConnectivity);
            quadConnectivityRow = quadElementsConnectivityDimen(1);
            quadConnectivityClmns = quadElementsConnectivityDimen(2);
            mappedQuadElementConnectivity = [];
            for jr = 1:quadConnectivityRow
                for jc = 1:quadConnectivityClmns
                    mappedQuadElementConnectivity(jr,jc) = M1(obj.quadElementsConnectivity(jr,jc));
                end
            end 
          quadplot(mappedQuadElementConnectivity,QuadElements_X_Cords,QuadElements_Y_Cords)
          hold off
        end
        function [vecSegmentXs, vecSegmentYs, vecSegmentNodeIDs] = GetDomainEdges(obj)
            if ((length(obj.twoSidedEdgeGroups) == 0) || (length(obj.twoSidedEdgeGroups{1}.twoSidedEdges) == 0))
                vecSegmentXs = cell(0);
                vecSegmentYs = cell(0);
                vecSegmentNodeIDs = cell(0);
                return;
            end
            sz = length(obj.twoSidedEdgeGroups{1}.twoSidedEdges);
            vecSegmentXs = cell(sz, 1);
            vecSegmentYs = cell(sz, 1);
            vecSegmentNodeIDs = cell(sz, 1);
            for i = 1:sz
                vecSegmentX(1) = obj.twoSidedEdgeGroups{1}.twoSidedEdges{i}.insideEdge.edge_node1.crd(1);
                vecSegmentX(2) = obj.twoSidedEdgeGroups{1}.twoSidedEdges{i}.insideEdge.edge_node2.crd(1);
                vecSegmentY(1) = obj.twoSidedEdgeGroups{1}.twoSidedEdges{i}.insideEdge.edge_node1.crd(2);
                vecSegmentY(2) = obj.twoSidedEdgeGroups{1}.twoSidedEdges{i}.insideEdge.edge_node2.crd(2);
                vecSegmentNodeID(1) = obj.twoSidedEdgeGroups{1}.twoSidedEdges{i}.insideEdge.edge_node1.n_id;
                vecSegmentNodeID(2) = obj.twoSidedEdgeGroups{1}.twoSidedEdges{i}.insideEdge.edge_node2.n_id;

                vecSegmentXs{i} = vecSegmentX;
                vecSegmentYs{i} = vecSegmentY;
                vecSegmentNodeIDs{i} = vecSegmentNodeID;
            end
        end
    end
end