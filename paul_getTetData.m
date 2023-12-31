function data = paul_getTetData(T,X,lite,force)

    if(nargin == 2)
        lite = true;
        force = false;
    elseif(nargin == 3)
        force = false;
    end
    
    % lite means don't load all properties. force means don't set to lite
    % even if size is far beyond recommended. will take a long time to run.
    if ((size(T,1) > 330000 || size(X,1) > 61000) && ~force)
        'SETTING TO LITE. TOO MANY VERTICES/TETS.'
        lite = 1;
    end

    data.tetrahedra = T;
    data.vertices = X;
    data.numVertices = size(X,1);
    data.numTetrahedra = size(T,1);

    % 123 214 134 324
    tri = [T(:,1) T(:,2) T(:,3) ; T(:,2) T(:,1) T(:,4) ;
           T(:,1) T(:,3) T(:,4) ; T(:,3) T(:,2) T(:,4) ];

    [~,IA,IC] = unique(sort(tri,2),'rows');
    data.triangles = tri(IA,:);
    data.numTriangles = size(data.triangles,1);
    data.tetsToTriangles = reshape(IC,data.numTetrahedra,4);

    counts = accumarray(IC,ones(size(IC)),[data.numTriangles 1]);
    data.isBoundaryTriangle = double(counts == 1); % 16653:16655
    bt = data.triangles(find(data.isBoundaryTriangle),:);
    data.isBoundaryVertex = zeros(1,data.numVertices); data.isBoundaryVertex(unique(sort(bt(:))))=1;
    data.boundaryTriangles = find(data.isBoundaryTriangle);

    tt = data.triangles;
    data.triangleNormals = -cross(X(tt(:,2),:)-X(tt(:,1),:),X(tt(:,3),:)-X(tt(:,1),:));
    data.triangleAreas = sqrt(sum(data.triangleNormals.^2,2))/2;
    data.triangleNormals = bsxfun(@rdivide,data.triangleNormals,2*data.triangleAreas);

    n = data.numTetrahedra;

    data.tetBarycenters = .25*(X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:)+X(T(:,4),:));
    data.triangleBarycenters = (X(data.triangles(:,1),:)+X(data.triangles(:,2),:)+X(data.triangles(:,3),:))/3;

    E1 = X(T(:,2),:) - X(T(:,1),:);
    E2 = X(T(:,3),:) - X(T(:,1),:);
    E3 = X(T(:,4),:) - X(T(:,1),:);
    data.tetVolumes = abs(sum(cross(E1,E2).*E3,2)/6);

    %%
    % E is pairs of vertices for each triangle
    E = [data.triangles(:,1) data.triangles(:,2) ; data.triangles(:,2) data.triangles(:,3) ; data.triangles(:,3) data.triangles(:,1)];
    [~,IA,IC] = unique(sort(E,2),'rows');
    data.edges = E(IA,:);
    flipped = find(data.edges(:,1)>data.edges(:,2)); data.edges(flipped,:)=[data.edges(flipped,2) data.edges(flipped,1)];
    data.numEdges = size(data.edges,1);
    data.trianglesToEdges = reshape(IC,data.numTriangles,3);

    % get tetsToEdges
    E = [data.tetrahedra(:,1) data.tetrahedra(:,2) ; ...
        data.tetrahedra(:,1) data.tetrahedra(:,3) ; ...
        data.tetrahedra(:,1) data.tetrahedra(:,4) ; ...
        data.tetrahedra(:,2) data.tetrahedra(:,3) ; ...
        data.tetrahedra(:,2) data.tetrahedra(:,4) ; ...
        data.tetrahedra(:,3) data.tetrahedra(:,4)];

    [~,IA,IC] = unique(sort(E,2),'rows');
    TetEdges = E(IA,:);
    data.tetsToEdges = reshape(IC,data.numTetrahedra,6);

    % check that tets to edges is right. can make check random for speed. Or
    % prove but im not sure...
    tetsToEdgesWithDup = [reshape(data.trianglesToEdges(data.tetsToTriangles,1),size(data.tetsToTriangles)) ...
        reshape(data.trianglesToEdges(data.tetsToTriangles,2),size(data.tetsToTriangles)) ...
        reshape(data.trianglesToEdges(data.tetsToTriangles,3),size(data.tetsToTriangles))];
    
    %for j = 1:data.numTetrahedra
    %    assert(sum(sort(unique(tetsToEdgesWithDup(j,:))) == sort(data.tetsToEdges(j,:)))==6);
    %end

    data.isBoundaryEdge = zeros(data.numEdges,1);
    for i=1:3
        data.isBoundaryEdge(data.trianglesToEdges(data.boundaryTriangles,i)) = 1;
    end
    data.boundaryEdges = find(data.isBoundaryEdge);
    adj = sparse(repmat((1:n)',4,1),data.tetsToTriangles(:),ones(4*n,1));

    % get trianglesToTets
    if(~lite)
        numberedAdj = adj.*[1:size(adj,1)]';
        collapsedAdj = numberedAdj(find(numberedAdj~=0));
        trianglesToTets = mat2cell(full(collapsedAdj), full(sum(adj)));
        data.trianglesToTets = trianglesToTets;
        
        %% compute nonboundaryTriToTets
        nonBoundaryTrianglesToTets = cell2mat(data.trianglesToTets(find(~data.isBoundaryTriangle)));
        nonBoundaryTrianglesToTets = reshape(nonBoundaryTrianglesToTets,2,numel(nonBoundaryTrianglesToTets)/2)';
        data.nonBoundaryTrianglesToTets =nonBoundaryTrianglesToTets ;
        
        dualGraphAdj = sparse(nonBoundaryTrianglesToTets(:,1),nonBoundaryTrianglesToTets(:,2),ones(size(nonBoundaryTrianglesToTets(:,1),1),1),data.numTetrahedra,data.numTetrahedra);
        dualGraphAdj = dualGraphAdj + dualGraphAdj' ~= 0;
        data.dualGraphAdjacency = dualGraphAdj;
        data.dualGraph = graph(data.dualGraphAdjacency);
    end
    
    % verts to tets
    tetsToVertsIndicator = sparse(repmat((1:data.numTetrahedra)',4,1),data.tetrahedra(:),ones(4*data.numTetrahedra,1),data.numTetrahedra, data.numVertices);
    data.tetsToVertsIndicator = tetsToVertsIndicator;
    
    % get edgesToTets
    tetsToEdgesIndicator = sparse(repmat((1:data.numTetrahedra)',6,1),data.tetsToEdges(:),ones(6*data.numTetrahedra,1),data.numTetrahedra, data.numEdges);
    data.tetsToEdgesIndicator = tetsToEdgesIndicator;
    if(~lite)
        
        numberedTetsToEdgesIndicator = tetsToEdgesIndicator.*[1:size(tetsToEdgesIndicator,1)]';
        collapsednumberedTetsToEdgesIndicator = numberedTetsToEdgesIndicator(find(numberedTetsToEdgesIndicator~=0));

        edgesToTets = mat2cell(full(collapsednumberedTetsToEdgesIndicator), full(sum(tetsToEdgesIndicator)));
        data.edgesToTets = edgesToTets;
    end
    
    % get edgesToTrianglesUnoriented
    trianglesToEdgesIndicator = sparse(repmat((1:data.numTriangles)',3,1),data.trianglesToEdges(:),repmat((1:data.numTriangles)',3,1));
    data.edgesToTrianglesIndicator = trianglesToEdgesIndicator';
    if(~lite)
        
        collapsedTrianglesToEdgesIndicator = trianglesToEdgesIndicator(find(trianglesToEdgesIndicator~=0));
    
        edgesToTriangles = mat2cell(full(collapsedTrianglesToEdgesIndicator), full(sum(trianglesToEdgesIndicator~=0)));
        data.edgesToTrianglesUnoriented = cellfun(@transpose,edgesToTriangles,'un',0);
    end
    
    %nonboundary triangle indices
    triangleSubindices = find(~data.isBoundaryTriangle);
    adj = adj(:,triangleSubindices);
    [Itet,Jtri,~] = find(adj);
    Jtri = triangleSubindices(Jtri);
    triTet = sortrows([Jtri Itet]);

    assert(isequal(triTet(1:2:end,1),triTet(2:2:end,1)));

    % get boundaryTets
    data.isBoundaryTet = sum(data.isBoundaryTriangle(data.tetsToTriangles),2)~=0;

    %tet1, edge, tet2
    entries = [];
    for j=1:3
        e = data.trianglesToEdges(triTet(1:2:end,1),j);
        possibleEntries = [];
        possibleEntries = [possibleEntries; triTet(1:2:end,2) e triTet(2:2:end,2) ];
        possibleEntries = [possibleEntries; triTet(2:2:end,2) e triTet(1:2:end,2) ];
        %idx = ~data.isBoundaryEdge(e);
        %possibleEntries = possibleEntries([idx;idx],:);
        entries = [entries;possibleEntries];
    end
    mtx = sortrows(entries);
    
    % mtx = sortrows([plusTet commonEdge minusTet]);
    nextTet1 = sparse(mtx(1:2:end,1),mtx(1:2:end,2),mtx(1:2:end,3),data.numTetrahedra,data.numEdges);
    nextTet2 = sparse(mtx(2:2:end,1),mtx(2:2:end,2),mtx(2:2:end,3),data.numTetrahedra,data.numEdges);

    % tets cycles
    edgeCycles = cell(data.numEdges,1);
    ne = data.numEdges;
    isBoundaryEdge = data.isBoundaryEdge;
    if(~lite)
%         tic
%         for i=1:data.numEdges
%             if mod(i,1000) == 1
%                 fprintf('Edge %d of %d...\n',i,ne);
%             end
% 
%             if isBoundaryEdge(i)
%                continue
%             end
% 
%             rtri = data.edgesToTrianglesUnoriented{i}(1);
%             adjTets = data.edgesToTets{i};
%             tets = data.trianglesToTets{rtri};
%             [~,ind] = ismember(tets,adjTets);
%             subG = subgraph(data.dualGraph, adjTets);
% %             subA = data.dualGraphAdjacency(adjTets,adjTets);
% %             subG = graph(subA); subG = subG.rmedge(ind(1),ind(2));
%             subG = subG.rmedge(ind(1),ind(2));
%             cycle = shortestpath(subG,ind(1),ind(2));
%             cycle = adjTets(cycle);
%             
%             assert(length(cycle) > 2);
%             edgeCycles{i} = cycle;
%         end
%         toc
%         tic
        for i=1:data.numEdges
            if mod(i,1000) == 1
                fprintf('Edge %d of %d...\n',i,ne);
            end

            if isBoundaryEdge(i)
               continue
            end

            cycle = find(nextTet1(:,i),1);
            cycle = [cycle nextTet1(cycle,i)];
            while (cycle(end) ~= cycle(1))
                if nextTet1(cycle(end),i) == cycle(end-1)
                    cycle = [cycle nextTet2(cycle(end),i)];
                else
                    cycle = [cycle nextTet1(cycle(end),i)];
                end
            end
            assert(length(cycle) > 2);
            edgeCycles{i} = cycle;
        end
%         toc
        

        for i=1:data.numEdges
            if mod(i,1000) == 1
                fprintf('Edge %d of %d...\n',i,ne);
            end

            if ~isBoundaryEdge(i)
               continue
            end

            unorderedTets = data.edgesToTets{i};
            if(numel(unorderedTets)==1)
                cycle = unorderedTets;
                edgeCycles{i} = cycle;
                continue;
            end
            % since tets can have multiple faces on boundary, generally, a boundary edge could have
            % more than 2 'boundary' tets.
            boundaryTets = unorderedTets(find(data.isBoundaryTet(unorderedTets)));

            % get the boundary tets that start and end the cycle for edge i.
            linedUpEdgeMatch = reshape(data.trianglesToEdges(data.tetsToTriangles(boundaryTets,:)',:)',12,numel(boundaryTets))' == i;
            TriangleWithBoundaryEdge = [sum(linedUpEdgeMatch(:,1:3),2) sum(linedUpEdgeMatch(:,4:6),2) sum(linedUpEdgeMatch(:,7:9),2) sum(linedUpEdgeMatch(:,10:12),2)] > 0;
            BoundaryTriangles = data.isBoundaryTriangle(data.tetsToTriangles(boundaryTets,:));
            TrueBoundaryTets = boundaryTets(find(sum(TriangleWithBoundaryEdge & BoundaryTriangles, 2)>0));
            assert(numel(TrueBoundaryTets)==2); % these boundary tets are the start and end of the cycle for edge i

            cycle = TrueBoundaryTets(1);

            while numel(cycle) ~= numel(unorderedTets)
                currentTet = cycle(end);
                triangs = data.tetsToTriangles(currentTet,:);
                inttri = triangs(find(~data.isBoundaryTriangle(triangs)));
                exttri = triangs(find(data.isBoundaryTriangle(triangs)));
                neighborTets = unique([reshape([data.trianglesToTets{inttri,:}],1,2*numel(inttri)) data.trianglesToTets{exttri}]);
                neighborTets(find(neighborTets == currentTet)) = [];
                viableNextTets = intersect(unorderedTets,neighborTets); % one or two tets
                assert(numel(viableNextTets)<=2);
                if(sum(find(cycle==viableNextTets(1)))>0)
                    cycle = [cycle viableNextTets(2)];
                else
                    cycle = [cycle viableNextTets(1)];
                end
                edgeCycles{i} = cycle;
            end
            assert(data.isBoundaryTet(cycle(end)));

        end

        data.edgeCycles = edgeCycles;

        isOrientedEdge = zeros(data.numEdges, 1);
        isOrientedTriangle = zeros(data.numTriangles, 1);
        startTriangle = find(~data.isBoundaryTriangle,1);
        startEdge = data.trianglesToEdges(startTriangle,1);
        isOrientedEdge(startEdge) = true;
        trianglesToTraverse = data.edgesToTrianglesUnoriented{startEdge};
        while numel(trianglesToTraverse)~=0
            currentTriangle = trianglesToTraverse(1);

            trianglesToTraverse = trianglesToTraverse(2:end);
            if(isOrientedTriangle(currentTriangle))
                continue;
            end

            verts = data.triangles(currentTriangle, :);
            edges = data.trianglesToEdges(currentTriangle, :);
            edgeverts = data.edges(edges,:);

            %if data.edgecycles{k} has only one tet, orientation is not well
            %defined.
            badOrientation = [numel(data.edgeCycles{edges(1)}) numel(data.edgeCycles{edges(2)}) numel(data.edgeCycles{edges(3)})]==1;

            orientedEdgeInd = find(isOrientedEdge(edges)' & ~badOrientation,1);
            assert(numel(orientedEdgeInd)==1);
            remainingEdges = edges(find([1:3]~=orientedEdgeInd));
            tets = data.trianglesToTets{currentTriangle};

            edgeOrientations = [findOrderOfAInB(edgeverts(1,:), verts) findOrderOfAInB(edgeverts(2,:), verts) findOrderOfAInB(edgeverts(3,:), verts)];
            forwardEdgeOrientation = edgeOrientations(orientedEdgeInd);

            tetOrientations = [findOrderOfAInB(tets, data.edgeCycles{edges(1)}) findOrderOfAInB(tets, data.edgeCycles{edges(2)}) findOrderOfAInB(tets, data.edgeCycles{edges(3)})];

            forwardTetOrientation = tetOrientations(orientedEdgeInd);
            ind1 = find(edges==remainingEdges(1));
            ind2 = find(edges==remainingEdges(2));

            if(~isOrientedEdge(edges(ind1)) && ...
                (   xor(tetOrientations(ind1), edgeOrientations(ind1)) ~= xor(forwardTetOrientation, forwardEdgeOrientation) ...
                    || numel(data.edgeCycles{edges(ind1)}) == 1 ))
                data.edgeCycles{edges(ind1)} = fliplr(data.edgeCycles{edges(ind1)});
            end

            if(~isOrientedEdge(edges(ind2)) && ...
                (xor(tetOrientations(ind2), edgeOrientations(ind2)) ~= xor(forwardTetOrientation, forwardEdgeOrientation) ...
                || numel(data.edgeCycles{edges(ind2)}) == 1))
                data.edgeCycles{edges(ind2)} = fliplr(data.edgeCycles{edges(ind2)});
            end

            % double check orientation is right now.
            tetOrientations = [findOrderOfAInB(tets, data.edgeCycles{edges(1)}) findOrderOfAInB(tets, data.edgeCycles{edges(2)}) findOrderOfAInB(tets, data.edgeCycles{edges(3)})];
            assert(badOrientation(ind1) || ~(xor(tetOrientations(ind1), edgeOrientations(ind1)) ~= xor(forwardTetOrientation, forwardEdgeOrientation)) || numel(data.edgeCycles{edges(ind1)}) == 1);
            assert(badOrientation(ind2) || ~(xor(tetOrientations(ind2), edgeOrientations(ind2)) ~= xor(forwardTetOrientation, forwardEdgeOrientation)) || numel(data.edgeCycles{edges(ind2)}) == 1);

            isOrientedEdge(edges(ind2))=true;
            isOrientedEdge(edges(ind1))=true;

            isOrientedTriangle(currentTriangle)=true;

            % badorientation edges don't need to have their triangles added to
            % the traversal list. They will be added when one of their non-bad edges is
            % oriented. If they are never added, that means all 3 of their
            % edges were badorienation, meaning had only 1 adjacent tet. This
            % corresponds to only one scenario: the tet is a single alone tet.
            % we're not interested in such a simple case.
            moreToTraverse = unique([data.edgesToTrianglesUnoriented{edges(find(~badOrientation))}]);
            moreToTraverse = moreToTraverse(~isOrientedTriangle(moreToTraverse));
            trianglesToTraverse = [trianglesToTraverse moreToTraverse];
        end
    
        %% verify that cycles are correctly oriented for interior edges
        for triIter = [1:data.numTriangles]
            edges = data.trianglesToEdges(triIter, :);
            if data.isBoundaryTriangle(triIter) || sum(data.isBoundaryEdge(edges))>0
                continue
            end

            edges = data.trianglesToEdges(triIter, :);
            edgeverts = data.edges(edges,:);
            verts = data.triangles(triIter, :);

            orient1 = [findOrderOfAInB(edgeverts(1,:), verts) findOrderOfAInB(edgeverts(2,:), verts) findOrderOfAInB(edgeverts(3,:), verts)];

            tets = data.trianglesToTets{triIter};

            tets1 = data.edgeCycles{edges(1)};
            tets2 = data.edgeCycles{edges(2)};
            tets3 = data.edgeCycles{edges(3)};

            orient2 = [findOrderOfAInB(tets, tets1) findOrderOfAInB(tets, tets2) findOrderOfAInB(tets, tets3)];

            assert(sum(orient1 == orient2)==3 | sum(orient1 == orient2)==0);
        end
    
        
        
        for i = 1:numel(data.edgeCycles)
            isbe = data.isBoundaryEdge(i);
            tetCycle = data.edgeCycles{i};
            neighboringTris = data.edgesToTrianglesUnoriented{i};
            ntris = numel(tetCycle)+1;
            tets2tris = data.tetsToTriangles(tetCycle,:);
            
            if(isbe)
                if(numel(tetCycle)==1)
                    % trying to carry orientation from tets. 1 tet alone
                    % has no orientation. gonna be lazy and hope users
                    % subdivide.
                    display('trouble! cannot have tet with 2 boundary faces!');
                end
                
                startTris = intersect(tets2tris(1,:),neighboringTris);
                startTri = startTris(find(data.isBoundaryTriangle(startTris)));
                triCycle = [startTri ARemoveB(startTris,startTri)];
                
                for j = 2:size(tets2tris,1)
                    tris = intersect(tets2tris(j,:),neighboringTris);
                    nextTri = ARemoveB(tris, triCycle(j));
                    triCycle = [triCycle nextTri];
                end
            else
                startTris = intersect(tets2tris(1,:),neighboringTris);
                nextTri = intersect(tets2tris(1,:),tets2tris(2,:));
                triCycle = [ARemoveB(startTris,nextTri) nextTri];
                for j = 2:(size(tets2tris,1)-1)
                    tris = intersect(tets2tris(j,:),neighboringTris);
                    nextTri = ARemoveB(tris, triCycle(j));
                    triCycle = [triCycle nextTri];
                end
            end
            data.edgeTriCycles{i}=triCycle;
        end
        
        % verify edge TriCycles
%         for i = 1:numel(data.edgeTriCycles)
%             triCycle = data.edgeTriCycles{i};
%             vs = data.vertices(data.edges(i,:),:); 
%             %plot3(vs(:,1),vs(:,2),vs(:,3),'r'); 
%             hold on; axis equal;
%             triloop = data.triangleBarycenters(triCycle,:);
%             %plot3(triloop(:,1),triloop(:,2),triloop(:,3),'r');
%             triloop = .2*triloop + .9*sum(vs)/2;
%             triverts = data.vertices(data.triangles(triCycle,:),:);
%             scatter3(triverts(:,1),triverts(:,2),triverts(:,3));
%             plot3(triloop(:,1),triloop(:,2),triloop(:,3),'b');
%             %hold off;
%         end
        
    end
    
    %% separate nonboundary and boundary edges
    data.NonBoundaryEdges = data.edges(find(~data.isBoundaryEdge),:);
    data.BoundaryEdges = data.edges(find(data.isBoundaryEdge),:);
    
    
    %% compute primal spanning tree of the volume. spans vertices
    % data.PrimalVolumeVertexSpanningTree = PrimalVolumeVertexSpanningTree(data.edges);
    
    % compute dual spanning tree of the volume. spans tets
    % data.DualVolumeVertexSpanningTree = DualVolumeVertexSpanningTree(data);
    
    % compute primal spanning tree of volume WITHOUT BOUNDARY EDGES
    %Inds=find(~data.isBoundaryEdge); Inds = Inds(PrimalVolumeVertexSpanningTree(data.NonBoundaryEdges));
    %data.BoundaryLessPrimalSpanningTreeRelToEdges = Inds;
    
    assert(all(data.edges(:,1)<=data.edges(:,2)));
    
    
    
    % compute incidence matrix
    data.primalIncidenceMatrix = sparse([1:data.numEdges 1:data.numEdges]', [data.edges(:,1);data.edges(:,2)], [ones(data.numEdges,1) -ones(data.numEdges,1)]    ,data.numEdges,data.numVertices);
    
    % compute primal vert/edge laplacian
    data.primalOneLaplacian = data.primalIncidenceMatrix'*data.primalIncidenceMatrix;
    
end
