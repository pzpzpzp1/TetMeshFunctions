function EdgesRemaining = paul_dualSurfaceTest(accumCurveEdges)

    [X,T]=paul_loadTetGenMesh('./meshes/sphere_61k');
    data = paul_getTetData(T,X,0,1);

    %% create test singular primal curve

    if accumCurveEdges == 0
        accumCurveEdges={};
        [curveEdges, curveV, f] = chooseCurve(data, 0, 0, []); accumCurveEdges{1} = curveEdges;
        [curveEdges2, curveV2, f] = chooseCurve(data, f, 'b', accumCurveEdges); accumCurveEdges{2} = curveEdges2;

        close;
    elseif  accumCurveEdges == 1
        load curveV.mat;
        load curveV2.mat;
        load curveEdges.mat;
        load curveEdges2.mat;
    else
        % column vectors list of edge indices for the mesh.
        curveEdges = accumCurveEdges{1};
        curveEdges2 = accumCurveEdges{2};
        curveV = unique(data.edges(curveEdges,:))';
        curveV2 = unique(data.edges(curveEdges2,:))';
    end

    %{ 
    % just a peice to use for display
    f = figure(); hold off; scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),.01,'b'); hold on;
    for i = 1:numel(accumCurveEdges)
        VisualizeEdges(accumCurveEdges{i}, data, '.', f)
    end
    %}

    %% compute tree with boundary, and without edges adj to curve, and with curve, and with interfering curve.
    v2v = sparse(data.edges(:,1),data.edges(:,2),1:size(data.edges,1)); v2v(data.numVertices+1,data.numVertices+1)=0;
    es = [v2v(curveV,:) v2v(:,curveV)'];
    edgesAdjToCurve = full(sort(unique(es(find(es~=0)))));
    exclude = unique([find(data.isBoundaryEdge); edgesAdjToCurve]);
    exclude = edgesAdjToCurve; % include the boundary in the tree. empirically, decreases non manifold areas at boundary.
    include = curveEdges2;
    vinds = 1:data.numVertices; vinds(curveV)=[];
    startVert = vinds(1);
    data.BoundaryLessPrimalSpanningTreeRelToEdges = PrimalVolumeVertexSpanningTreeMoreOps(data.edges, exclude, include, []);
    assert(numel(intersect(data.BoundaryLessPrimalSpanningTreeRelToEdges, curveEdges))==0); % by construction, curve shouldn't be in the tree yet.
    data.BoundaryLessPrimalSpanningTreeRelToEdges = [data.BoundaryLessPrimalSpanningTreeRelToEdges curveEdges'];
    % VisualizeGraph(data.edges, data.vertices, 'r', 1, data.BoundaryLessPrimalSpanningTreeRelToEdges)    

    % Naive method to avoid 2nd curve. This method fails. Adding after the fact is
    % unreliable and can result in no solution found even when there is one.
    % data.BoundaryLessPrimalSpanningTreeRelToEdges = [data.BoundaryLessPrimalSpanningTreeRelToEdges curveEdges2']; 
    
    
    %% Try to close dual edges until no more can be closed unless they are adjacent to primal curve
    dualEdges = 1:data.numTriangles;
    isInCurve = sparse(curveV, ones(numel(curveV),1), ones(numel(curveV),1)); isInCurve(data.numVertices+1)=0;
    dualEdgesCloseToCurve = find(sum(isInCurve(data.triangles),2)~=0);
    dualEdgesToClose = dualEdges; dualEdgesToClose(dualEdgesCloseToCurve)=[];

    % find the dual edges that are one edge away from closing.
    GrowingTree = data.BoundaryLessPrimalSpanningTreeRelToEdges;
    inTreeIndicator = sparse(GrowingTree(:), ones(numel(GrowingTree),1),ones(numel(GrowingTree),1)); inTreeIndicator(data.numEdges+1)=0;
    readyToClose = find(sum(inTreeIndicator(data.trianglesToEdges(dualEdgesToClose,:)),2)==2);
    % readytoClose is an index into data.triangles(dualEdgesToClose,:));

    % keep closing triangles until there are no more that can be closed
    while(numel(readyToClose) ~= 0)
        trisToClose = data.trianglesToEdges(dualEdgesToClose(readyToClose),:); % vertinds!
        edgeIndsToCloseTris = find(~inTreeIndicator(data.trianglesToEdges(dualEdgesToClose(readyToClose),:)));
        edgesToCloseTris = trisToClose(edgeIndsToCloseTris);

        GrowingTree = unique([GrowingTree edgesToCloseTris']); inTreeIndicator(edgesToCloseTris)=1;

        dualEdgesToClose = ARemoveB(dualEdgesToClose, dualEdgesToClose(readyToClose));
        readyToClose = find(sum(inTreeIndicator(data.trianglesToEdges(dualEdgesToClose,:)),2)==2);
    end

    % compute dual surface boundaries. (dual edges = triangles)
    % assert these triangles are all adj to boundary or curveV
    DualSurfaceBoundaries = find(sum(inTreeIndicator(data.trianglesToEdges),2)==2);
    for i = 1:numel(DualSurfaceBoundaries)
        tri = data.triangles(DualSurfaceBoundaries(i),:);
        inter = intersect([find(data.isBoundaryVertex) curveV], tri);
        if(numel(inter) == 0)
            DualSurfaceBoundaries(i) % triangle
            assert(0==1);
        end
    end
    size(GrowingTree)

    %% Visualize remaining dual surface
    EdgesRemaining = 1:data.numEdges;
    EdgesRemaining(GrowingTree)=[];
    f = VisualizeDualSurface(data, EdgesRemaining, 'p');
    VisualizeEdges(curveEdges, data, '-',f);

    %% Visualize non manifold locations: this should be very low. possibly empty
    indexNonManifoldTriangles = find(sum(inTreeIndicator(data.trianglesToEdges),2)==0);
    indexNonManifoldBoundaryTriangles = indexNonManifoldTriangles(find(data.isBoundaryTriangle(indexNonManifoldTriangles)));
    m = reshape(cell2mat(data.trianglesToTets(indexNonManifoldBoundaryTriangles)),1,numel(indexNonManifoldBoundaryTriangles));
    p1 = data.tetBarycenters(m,:);
    p2 = data.triangleBarycenters(indexNonManifoldBoundaryTriangles,:);
    interp = [0:.05:1]; final = [];
    for i = 1:numel(interp)
        inter = interp(i) * p1 + p2 * (1-interp(i)); final = [final; inter];
    end
    figure(f); hold on; axis equal; 
    scatter3(final(:,1),final(:,2),final(:,3),1,'red');

    indexNonManifoldIntTriangles = indexNonManifoldTriangles(find(~data.isBoundaryTriangle(indexNonManifoldTriangles)));
    m = reshape(cell2mat(data.trianglesToTets(indexNonManifoldIntTriangles)),2,numel(indexNonManifoldIntTriangles));
    p1 = data.tetBarycenters(m(1,:)',:);
    p2 = data.tetBarycenters(m(2,:)',:);
    interp = [0:.05:1]; final = [];
    for i = 1:numel(interp)
        inter = interp(i) * p1 + p2 * (1-interp(i)); final = [final; inter];
    end
    figure(f); hold on; axis equal; 
    scatter3(final(:,1),final(:,2),final(:,3),2,'red');
    VisualizeEdges(curveEdges, data, '-', f);
    VisualizeEdges(curveEdges2, data, '.', f);

    %% Visualize boundary edges of the dual surface (that dont hit boundary of volume) Usually no contribution.
    DualSurfaceBoundaries = find(sum(inTreeIndicator(data.trianglesToEdges),2)==2);
    DualSurfaceBoundaries = ARemoveB(DualSurfaceBoundaries', find(data.isBoundaryTriangle)');
    m = reshape(cell2mat(data.trianglesToTets(DualSurfaceBoundaries)),2,numel(DualSurfaceBoundaries));
    p1 = data.tetBarycenters(m(1,:)',:);
    p2 = data.tetBarycenters(m(2,:)',:);
    interp = [0:.05:1]; final = [];
    for i = 1:numel(interp)
        inter = interp(i) * p1 + p2 * (1-interp(i)); final = [final; inter];
    end
    figure(f); hold on; axis equal; 
    scatter3(final(:,1),final(:,2),final(:,3),2,'red');
    VisualizeEdges(curveEdges, data, '-', f);

    %% correspond to dual surface intersecting volume boundary.
    boundaryEdgesNotInTree = ARemoveB(find(data.isBoundaryEdge)', GrowingTree);
    btris = [];
    for i = 1:numel(boundaryEdgesNotInTree)
        triangleFan = data.edgesToTrianglesUnoriented{boundaryEdgesNotInTree(i)};
        % btris should always be n x 2
        btris = [btris; triangleFan(find(data.isBoundaryTriangle(triangleFan)))];
    end
    p1 = data.triangleBarycenters(btris(:,1),:); p2 = data.triangleBarycenters(btris(:,2),:);
    interp = 0:.05:1; final = [];
    for i = 1:numel(interp)
        final = [final; p1*interp(i)+p2*(1-interp(i))];
    end
    scatter3(final(:,1),final(:,2),final(:,3),2,'blue');

    % boundary edges that aren't in tree make dual surface boundaries on the
    % primal surface. Method seems to leave none of these.
    boundaryTrianglesWhoseDualEdgesAreBoundaries = find(sum(inTreeIndicator(data.trianglesToEdges(find(data.isBoundaryTriangle),:)),2)==2);
    assert(numel(boundaryTrianglesWhoseDualEdgesAreBoundaries)==0);

    % Sanity check. Number of edges adj to curve that arent the curve that end up in the
    % tree. This is provably 0. consider any triangle adj to the curve. it has
    % at least 2 edges in the edgesAdjToCurve set, which are not an will not
    % ever be in the tree. so these edges will never get added to the tree.
    % Also these edges are explicitly not allowed to close.
    assert(numel(ARemoveB(intersect(edgesAdjToCurve,GrowingTree)',curveEdges'))==0);

    %% Simple Visualize remaining surface
    %{
    EdgesRemaining = 1:data.numEdges;
    EdgesRemaining(GrowingTree)=[];
    f=VisualizeDualSurface(data, EdgesRemaining, '.');
    VisualizeEdges(curveEdges, data, '.',f);
    %}


end