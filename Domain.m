classdef Domain
    properties
        grains = cell(0);
        all_grainIDs;
        num_grains;
        grainID2Index = [];
        grainID_2NeighborGrainIDs = cell(0);
        domainxm;
        domainxM;
        domainym;
        domainyM;
        % this boolean should belong elsewhere, but for convenience added
        % here. It will offset coordinate so all SVEs have lower left
        % (around) 0, 0
        intersectedDomain_crd_offset = 1;
        
        intersected_xlim;
        intersected_ylim;
    end
    methods
        function objout = readBeforeMaterialAssembly(obj, fid)
            objout = obj;
            buf = fscanf(fid, '%s', 1);
            while ((strcmp(buf, '*Assembly,') == 0) && (strcmp(buf, 'ASSEMBLY') == 0))
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
                    if (strcmp(buf, 'ASSEMBLY') == 1)
                        objout.num_grains = length(obj.grains);
                        return;
                    end
                    buf
                    fprintf(1, 'a bug in reading ... as each block should start with *\n');
                    pause;
                end
                if (strcmp(buf, '*HEADING') == 1) % reading new grain
                    buf = fscanf(fid, '%s', 1);
                elseif (strcmp(buf, '*Part,') == 1) % reading new grain
                    buf = fscanf(fid, '%s', 1);
                    [a, b] = strtok(buf, 'GRAIN');
                    [c, d] = strtok(b, '-');
                    grainID = -str2num(d);
                    index = length(objout.grains) + 1;
                    objout.grainID2Index(grainID) = index;
                    objout.grains{index} = Grain;
                    objout.grains{index} = objout.grains{index}.read(fid, grainID);
                    objout.all_grainIDs(index) = grainID;
                    objout.grainID_2NeighborGrainIDs{grainID} = objout.grains{index}.gb_neighbor_ids;
                    fprintf(1, 'grain %d read\n', grainID);
                    buf = fscanf(fid, '%s', 1);
                else
                    buf = fscanf(fid, '%s', 1);
                end
            end
            objout.num_grains = length(objout.grains);
            [objout.domainxm, objout.domainxM, objout.domainym, objout.domainyM] = ...
            objout.ComputeDomainBounds();
        end
        function [xm, xM, ym, yM] = ComputeDomainBounds(obj)
            xm = obj.grains{1}.grainxm;
            xM = obj.grains{1}.grainxM;
            ym = obj.grains{1}.grainym;
            yM = obj.grains{1}.grainyM;
            for gi = 2:obj.num_grains
                xm = min(xm, obj.grains{gi}.grainxm);
                xM = max(xM, obj.grains{gi}.grainxM);
                ym = min(ym, obj.grains{gi}.grainym);
                yM = max(yM, obj.grains{gi}.grainyM);
            end
        end
        function vecGrainBoundaries = getVecGrainBoundaries(obj)
            sz = length(obj.grains);
            vecGrainBoundaries = cell(sz, 1);
            for i = 1:sz
                vecGrainBoundaries{i} = obj.grains{i}.getDomainEdges();
            end
        end
        function objout = FormAllGrainEdgeNeighborhoodInfo(obj, tol)
            if nargin < 2
                tol = 1e-6;
            end
            objout = obj;
            % setting domain boundary edges
            for grindi = 1:objout.num_grains
                 objout.grains{grindi}.twoSidedEdgeGroups = cell(1, 1);
                 objout.grains{grindi}.twoSidedEdgeGroups{1} = SetOfTwoSidedEdges2D;
            end
            % get all neighborhood info
        	vecGrainBoundaries = objout.getVecGrainBoundaries();
            for grindi = 1:objout.num_grains
                grIDi = objout.all_grainIDs(grindi);
                possible_neighborGrain_ids = objout.grainID_2NeighborGrainIDs{grIDi};
%                if (isempty(possible_neighborGrain_ids))
%                    possible_neighborGrain_ids = objout.all_grainIDs;
%                end
                
                num_grain_edge_i = length(vecGrainBoundaries{grindi});
                
                
                % looping over possible neighbor grain neighbor ids
                num_p_neighbors = length(possible_neighborGrain_ids);
                for grain_out_cntr = 1:num_p_neighbors 
                    grIDo = possible_neighborGrain_ids(grain_out_cntr);
                    grindo = objout.grainID2Index(grIDo);
                    
                    num_grain_edge_o = length(vecGrainBoundaries{grindo});
                    
                    % looping over edges in and see if they match the
                    % neighbor
                    for edge_i = 1:num_grain_edge_i
                        edgeIn = vecGrainBoundaries{grindi}{edge_i};
                        if (edgeIn.processed == 1)
                            continue;
                        end
                        
                        toBreak = 0;
                        for edge_o = 1:num_grain_edge_o
                            edgeOut = vecGrainBoundaries{grindo}{edge_o};
                            if (edgeOut.processed == 1)
                                continue;
                            end
                            % see if the two edges match
                            match = edgeIn.DoEdgesMath(edgeOut, tol);
                            if (match > 0) % the two edges match, so add them to corresponding lists
                                % first see if the neighborgroup set is in
                                % inside
                                sz_g_i = length(objout.grains{grindi}.twoSidedEdgeGroups);
                                groupIndexInside = -1;
                                for k = 1:sz_g_i
                                    if (objout.grains{grindi}.twoSidedEdgeGroups{k}.otherSide_grainID == grIDo)
                                        groupIndexInside = k;
                                        break;
                                    end
                                end
                                if (groupIndexInside == -1) % not found on the other side
                                    groupIndexInside = sz_g_i + 1;
                                    objout.grains{grindi}.twoSidedEdgeGroups{groupIndexInside} = SetOfTwoSidedEdges2D;
                                    objout.grains{grindi}.twoSidedEdgeGroups{groupIndexInside}.otherSide_grainID = grIDo;
                                end
                             
                                % second see if the other grain already has
                                % stored this grain (inside) as neightbor
                                sz_g_o = length(objout.grains{grindo}.twoSidedEdgeGroups);
                                groupIndexOutside = -1;
                                for k = 1:sz_g_o
                                    if (objout.grains{grindo}.twoSidedEdgeGroups{k}.otherSide_grainID == grIDi)
                                        groupIndexOutside = k;
                                        break;
                                    end
                                end
                                if (groupIndexOutside == -1) % not found
                                    groupIndexOutside = sz_g_o + 1;
                                    objout.grains{grindo}.twoSidedEdgeGroups{groupIndexOutside} = SetOfTwoSidedEdges2D;
                                    objout.grains{grindo}.twoSidedEdgeGroups{groupIndexOutside}.otherSide_grainID = grIDi;
                                end

                                objout.grains{grindi}.twoSidedEdgeGroups{groupIndexInside}.otherSide_grain_IndexInNeighboringSet = groupIndexOutside;
                                objout.grains{grindo}.twoSidedEdgeGroups{groupIndexOutside}.otherSide_grain_IndexInNeighboringSet = groupIndexInside;
                                
                                sz_two_edgesInside = length(objout.grains{grindi}.twoSidedEdgeGroups{groupIndexInside}.twoSidedEdges);
                                twoEdgeIndex_Inside = sz_two_edgesInside + 1;
                                sz_two_edgesOutside = length(objout.grains{grindo}.twoSidedEdgeGroups{groupIndexOutside}.twoSidedEdges);
                                twoEdgeIndex_Outside = sz_two_edgesOutside + 1;

                                % now adding edge to both sides
                                tmpTwoSidedEdge2D = TwoSidedEdge2D;
                                tmpTwoSidedEdge2D.match = match;
                                tmpTwoSidedEdge2D.insideEdge = vecGrainBoundaries{grindi}{edge_i};
                                tmpTwoSidedEdge2D.outsideEdge = vecGrainBoundaries{grindo}{edge_o};
                                tmpTwoSidedEdge2D.insideEdge_grain_ID = grIDi;
                                tmpTwoSidedEdge2D.outsideEdge_grain_ID = grIDo;
                                tmpTwoSidedEdge2D.insideEdge_grain_indexInDomain = grindi;
                                tmpTwoSidedEdge2D.outsideEdge_grain_indexInDomain = grindo;
                                tmpTwoSidedEdge2D.insideEdge_groupIndex_InGrain = twoEdgeIndex_Inside;
                                tmpTwoSidedEdge2D.outsideEdge_groupIndex_InGrain = twoEdgeIndex_Outside;
                                objout.grains{grindi}.twoSidedEdgeGroups{groupIndexInside}.twoSidedEdges{twoEdgeIndex_Inside} = tmpTwoSidedEdge2D;

                                % now from the other side ...
                                tmpTwoSidedEdge2D = TwoSidedEdge2D;
                                tmpTwoSidedEdge2D.match = match;
                                tmpTwoSidedEdge2D.insideEdge = vecGrainBoundaries{grindo}{edge_o};
                                tmpTwoSidedEdge2D.outsideEdge = vecGrainBoundaries{grindi}{edge_i};
                                tmpTwoSidedEdge2D.insideEdge_grain_ID = grIDo;
                                tmpTwoSidedEdge2D.outsideEdge_grain_ID = grIDi;
                                tmpTwoSidedEdge2D.insideEdge_grain_indexInDomain = grindo;
                                tmpTwoSidedEdge2D.outsideEdge_grain_indexInDomain = grindi;
                                tmpTwoSidedEdge2D.insideEdge_groupIndex_InGrain = twoEdgeIndex_Outside;
                                tmpTwoSidedEdge2D.outsideEdge_groupIndex_InGrain = twoEdgeIndex_Inside;
                                objout.grains{grindo}.twoSidedEdgeGroups{groupIndexOutside}.twoSidedEdges{twoEdgeIndex_Outside} = tmpTwoSidedEdge2D;

                                % now marking facet as being processed
                                vecGrainBoundaries{grindi}{edge_i}.processed = 1;
                                vecGrainBoundaries{grindo}{edge_o}.processed = 1;
                                toBreak = 1;
                            end
                            if (toBreak)
                                break;
                            end
                        end
                    end
                end
                for edge_i = 1:num_grain_edge_i
                    if (vecGrainBoundaries{grindi}{edge_i}.processed == 1)
                        continue;
                    end
                   groupIndexInside = 1; % for domain boundary
                   sz_two_edgesInside = length(objout.grains{grindi}.twoSidedEdgeGroups{groupIndexInside}.twoSidedEdges);
                   twoEdgeIndex_Inside = sz_two_edgesInside + 1;
                   tmpTwoSidedEdge2D = TwoSidedEdge2D;
                   tmpTwoSidedEdge2D.match = 1;
                   tmpTwoSidedEdge2D.insideEdge = vecGrainBoundaries{grindi}{edge_i};
                   tmpTwoSidedEdge2D.outsideEdge = vecGrainBoundaries{grindi}{edge_i};
                   tmpTwoSidedEdge2D.insideEdge_grain_ID = grIDi;
                   tmpTwoSidedEdge2D.outsideEdge_grain_ID = -1;
                   tmpTwoSidedEdge2D.insideEdge_grain_indexInDomain = grindi;
                   tmpTwoSidedEdge2D.outsideEdge_grain_indexInDomain = -1;
                   tmpTwoSidedEdge2D.insideEdge_groupIndex_InGrain = twoEdgeIndex_Inside;
                   tmpTwoSidedEdge2D.outsideEdge_groupIndex_InGrain = -1;
                   objout.grains{grindi}.twoSidedEdgeGroups{groupIndexInside}.twoSidedEdges{twoEdgeIndex_Inside} = tmpTwoSidedEdge2D;

                   % add nodes to domain boundary list
                   nd{1} = vecGrainBoundaries{grindi}{edge_i}.edge_node1;
                   nd{2} = vecGrainBoundaries{grindi}{edge_i}.edge_node2;
                   for kk = 1:2
                        fnd = find(objout.grains{grindi}.domainBoundaryNodeIDs == nd{kk}.n_id);
                        if (isempty(fnd))
                            szz = length(objout.grains{grindi}.domainBoundaryNodeIDs);
                            indn = szz + 1;
                            objout.grains{grindi}.domainBoundaryNodeIDs(indn) = nd{kk}.n_id;
                            objout.grains{grindi}.domainBoundaryNodeCrds(indn, 1) = nd{kk}.crd(1);
                            objout.grains{grindi}.domainBoundaryNodeCrds(indn, 2) = nd{kk}.crd(2);
                        end
                   end
                   
                   vecGrainBoundaries{grindi}{edge_i}.processed = 1;
                end
            end
        end
        function fill_visualize(obj, grainID)
            if nargin < 2
                grainID = -1;
            end
            if(grainID < 0)
                grainList = obj.all_grainIDs;
            else
                grainList = obj.grainID_2NeighborGrainIDs{grainID};
            end
            for gi = 1:obj.num_grains
                gid = obj.all_grainIDs(gi);
                found = (gid == grainID);
                if (~found)
                    found = find(grainList == gid);
                end
                if (isempty(found))
                    continue;
                end
                clr = rand(1,3);
                obj.grains{gi}.fillMesh(clr, clr);
                hold on;
                xAve = mean(obj.grains{gi}.nodeCrds(:, 1));
                yAve = mean(obj.grains{gi}.nodeCrds(:, 2));
                text(xAve, yAve, num2str(obj.all_grainIDs(gi)));
                hold on;
            end
        end
        function objout = Extract(obj, xmin, xmax, ymin, ymax, renumbering, tol, scales)
            if nargin < 6
                renumbering = 1;
            end
            if nargin < 7
                tol = 1e-6;
            end
            if nargin < 8
                scales = [1, 1];
            end
            objout = Domain;
            objout.intersected_xlim = [xmin, xmax];
            objout.intersected_ylim = [ymin, ymax];
            
            objout.num_grains = 0;
            for i = 1:obj.num_grains
            	[grainOut, hasIntersection] = obj.grains{i}.Extract(xmin, xmax, ymin, ymax, renumbering);
                if (~hasIntersection)
                    continue;
                end
                if (obj.intersectedDomain_crd_offset)
                    numNodes = grainOut.nodeNum;
                    for j = 1:numNodes
                        grainOut.nodeCrds(j, 1) = grainOut.nodeCrds(j, 1) - xmin;
                        grainOut.nodeCrds(j, 2) = grainOut.nodeCrds(j, 2) - ymin;
                    end
                end
                objout.num_grains = objout.num_grains + 1;
                objout.grains{objout.num_grains} = grainOut;
                objout.all_grainIDs(objout.num_grains) = grainOut.id;
                objout.grainID2Index(grainOut.id) = objout.num_grains;
                objout.grainID_2NeighborGrainIDs{grainOut.id} = grainOut.gb_neighbor_ids;
            end
            if (objout.num_grains == 0)
                return;
            end
            % reduce the number of neighborhoods if needed
            for gi = 1:objout.num_grains
                id = objout.grains{gi}.id;
                gb_neighbor_ids_BK = objout.grainID_2NeighborGrainIDs{id};
                objout.grainID_2NeighborGrainIDs{id} = [];
                sz = length(gb_neighbor_ids_BK);
                cntr = 0;
                for j = 1:sz
                    grainNeiborID = gb_neighbor_ids_BK(j);
                    found = find(objout.all_grainIDs == grainNeiborID);
                    if (found)
                        cntr = cntr + 1;
                        objout.grainID_2NeighborGrainIDs{id}(cntr) = grainNeiborID;
                    else
                        objout.grains{gi}.gb_neighbor_ids(j) = -1;
                    end
                end
            end

            scale_x = scales(1);
            scale_y = scales(2);
            for i = 1:objout.num_grains
                [m, n] = size(objout.grains{i}.nodeCrds);
                if (scale_x ~= 1.0)
                    for j = 1:m
                        objout.grains{i}.nodeCrds(j, 1) = scale_x * objout.grains{i}.nodeCrds(j, 1);
                    end
                    objout.grains{i}.grainxm = scale_x * objout.grains{i}.grainxm;
                    objout.grains{i}.grainxM = scale_x * objout.grains{i}.grainxM;
                end
                if (scale_y ~= 1.0)
                    for j = 1:m
                        objout.grains{i}.nodeCrds(j, 2) = scale_y * objout.grains{i}.nodeCrds(j, 2);
                    end
                    objout.grains{i}.grainym = scale_y * objout.grains{i}.grainym;
                    objout.grains{i}.grainyM = scale_y * objout.grains{i}.grainyM;
                end
            end
            if (scale_x ~= 1.0)
                    objout.intersected_xlim = scale_x * objout.intersected_xlim;
            end
            if (scale_y ~= 1.0)
                    objout.intersected_ylim = scale_y * objout.intersected_ylim;
            end
            
            
            % for inter-grain neighborhoods
            objout = objout.FormAllGrainEdgeNeighborhoodInfo(tol);
        end
 
        function [vecSegmentXss, vecSegmentYss, vecSegmentNodeIDss, grainIDs] = GetDomainEdges(obj)
            cntr = 0;
            for gi = 1:obj.num_grains
                [vecSegmentXs, vecSegmentYs, vecSegmentNodeIDs] = obj.grains{gi}.GetDomainEdges();
                if (length(vecSegmentXs) == 0)
                    continue;
                end
                cntr = cntr + 1;
                vecSegmentXss{cntr} = vecSegmentXs;
                vecSegmentYss{cntr} = vecSegmentYs;
                vecSegmentNodeIDss{cntr} = vecSegmentNodeIDs;
                grainIDs(cntr) = obj.all_grainIDs(gi);
            end
        end
        function [vecSegmentXss, vecSegmentYss, vecSegmentNodeIDss, grainIDs] = PlotGetDomainEdges(obj)
            [vecSegmentXss, vecSegmentYss, vecSegmentNodeIDss, grainIDs] = obj.GetDomainEdges();
            sz = length(vecSegmentXss);
            mrkrs = {'x', 'o', 'p', 's', '*', '+', '>', '<'};
            for i = 1:sz
                szb = length(vecSegmentXss{i});
                cntr = 0;
                xs = [];
                ys = [];
                for j = 1:szb
                    cntr = cntr + 1;
                    xs(cntr) = vecSegmentXss{i}{j}(1);
                    ys(cntr) = vecSegmentYss{i}{j}(1);
                    cntr = cntr + 1;
                    xs(cntr) = vecSegmentXss{i}{j}(2);
                    ys(cntr) = vecSegmentYss{i}{j}(2);
                    cntr = cntr + 1;
                    xs(cntr) = nan;
                    ys(cntr) = nan;
                end
                clr = rand(1,3);
                mrkI = mod(i - 1, 8) + 1;
                mrk = mrkrs{mrkI};
                plot(xs, ys, 'Color', clr, 'Marker', mrk);
                hold on;
            end
        end
        function writeBeforeBC(obj, fid, configuration)
            if nargin < 3
                configuration = OSU_SpecificOutputConfiguration;
            end
            configuration.print_segment(fid, '__part_0_header');
            for gi = 1:obj.num_grains
                gID = obj.all_grainIDs(gi);
                fprintf(fid, '*Part, name=GRAIN-%d\n', gID);
                obj.grains{gi}.write(fid);
            end
            % printing assembly
            configuration.print_segment(fid, '__part_1_assembly');
            for gi = 1:obj.num_grains
                gID = obj.all_grainIDs(gi);
                fprintf(fid, '*Instance, name=GRAIN-%d, part=GRAIN-%d\n', gID, gID);
                fprintf(fid, '*End Instance\n');
                if (gi < obj.num_grains)
                    fprintf(fid, '**\n');
                end
            end
            configuration.print_segment(fid, '__part_2_assembly');
            configuration.print_segment(fid, '__part_3_material');
            for gi = 1:obj.num_grains
                gID = obj.all_grainIDs(gi);
                configuration.print_grain_material_properties(gID, fid);
            end
            configuration.print_segment(fid, '__part_4_interface');
            for gi = 1:obj.num_grains
                set2Print = obj.grains{gi}.gb_neighbor_ids;
                gid = obj.grains{gi}.id;
                szz = length(set2Print);
                for j = 1:szz
                    otherGrainID = set2Print(j);
                    if (otherGrainID < gid)
                        continue;
                    end
                    % just to make sure it's in the final domain, for
                    % example, potentially the case that the boundary of domain right
                    % is between grains ...
                    found = find(obj.grainID_2NeighborGrainIDs{gid} == otherGrainID);
                    if (found)
                        fprintf(fid, 'GRAIN-%d.SURF-MASTER_%d_%d ,  GRAIN-%d.SURF-SLAVE_%d_%d , COHPROP\n', gid, gid, otherGrainID, otherGrainID, otherGrainID, gid);
                    end
                end
            end
            configuration.print_segment(fid, '__part_5_load_steps');
            configuration.print_segment(fid, '__part_6_BC');
        end
        function writeStrainBC(obj, fid, eps_xx, eps_yy, eps_xy, center_4DispCal)
            if nargin < 5
                center_4DispCal = [];
            end
            if (length(center_4DispCal) == 0)
                if ((length(obj.intersected_xlim) == 2) && (length(obj.intersected_ylim) == 2))
                    if (~obj.intersectedDomain_crd_offset)
                        center_4DispCal(1) = 0.5 * (obj.intersected_xlim(1) + obj.intersected_xlim(2));
                        center_4DispCal(2) = 0.5 * (obj.intersected_ylim(1) + obj.intersected_ylim(2));
                    else
                        center_4DispCal(1) = 0.5 * (obj.intersected_xlim(2) - obj.intersected_xlim(1));
                        center_4DispCal(2) = 0.5 * (obj.intersected_ylim(2) - obj.intersected_ylim(1));
                    end
                else
                    center_4DispCal = [0, 0];
                end
            end            
            for gi = 1:obj.num_grains
                gID = obj.all_grainIDs(gi);
                sz = length(obj.grains{gi}.domainBoundaryNodeIDs);
                for i = 1:sz;
                    nodeID = obj.grains{gi}.domainBoundaryNodeIDs(i);
                    node_x = obj.grains{gi}.domainBoundaryNodeCrds(i, 1);
                    node_y = obj.grains{gi}.domainBoundaryNodeCrds(i, 2);
                    fprintf(fid, 'GRAIN-%d.BOUNDARY_NODES_%d_%d, 1, 1', gID, i, gID);
                    x = node_x - center_4DispCal(1);
                    y = node_y - center_4DispCal(2);
                    u_x = eps_xx * x + eps_xy * y;
                    u_y = eps_xy * x + eps_yy * y;
                    fprintf(fid, ', %10.6f\n', u_x);
                    fprintf(fid, 'GRAIN-%d.BOUNDARY_NODES_%d_%d, 2, 2', gID, i, gID);
                    fprintf(fid, ', %10.6f\n', u_y);
                end
            end
        end
        function writeAfterBC(obj, fid, configuration)
            if nargin < 3
                configuration = OSU_SpecificOutputConfiguration;
            end
        configuration.print_segment(fid, '__part_7_output_request');
        end
        
        function writeDifferentBC_SVEs(obj, fileNameBaseWithSVENo, configuration, center_4DispCal)
            if nargin < 3
                configuration = OSU_SpecificOutputConfiguration;
            end
            if nargin < 4
                center_4DispCal = [];
            end
            [eps_xx_eps_yy_eps_xyS, addedNames] = configuration.getStrains();
            num_loadCases = length(addedNames);
            %
            % printing beginning part
            for filei = 1:num_loadCases
                fileNames{filei} = [fileNameBaseWithSVENo, addedNames{filei}, '.inp'];
            end
            fid = fopen(fileNames{1}, 'w');
        	obj.writeBeforeBC(fid, configuration);
            fclose(fid);
            for filei = 2:num_loadCases
                source = fileNames{1};
                destination = fileNames{filei};
                [status,msg,msgID] = copyfile(source, destination);
            end
            for filei = 1:num_loadCases
                fid = fopen(fileNames{filei}, 'a');
                eps_xx = eps_xx_eps_yy_eps_xyS(filei, 1);
                eps_yy = eps_xx_eps_yy_eps_xyS(filei, 2);
                eps_xy = eps_xx_eps_yy_eps_xyS(filei, 3);
                obj.writeStrainBC(fid, eps_xx, eps_yy, eps_xy, center_4DispCal);
                obj.writeAfterBC(fid, configuration);
                fclose(fid);
            end
        end
    end
end
