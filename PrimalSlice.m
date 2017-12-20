
%% run from TetMeshFunctions
flow = 0
if(flow == 1)
    path2HMeshFuncs='./../HexMeshSingularitySheetingTest';
    addpath(path2HMeshFuncs);

    % HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/bunny_ours.vtk']);
    % HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/bunny.vtk']);
    % HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/rocker_arm.vtk']);
    % HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/polycut2013/fertility.vtk']);
    % HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/sculpture-B_ours.vtk']);

    % HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/rod_ours.vtk']);
    % HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/fertility.vtk']);
    HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/double_torus.vtk']);
    %HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/ellipsoid-B.vtk']);

    TMesh = HexToTet(HMesh);
    X = TMesh.V2P;
    T = TMesh.T2V;
    data = paul_getTetData(T,X,0);
    
    %% load singular edges
    % scatter3(data.vertices(find(TMesh.isSingularVertex),1),data.vertices(find(TMesh.isSingularVertex),2),data.vertices(find(TMesh.isSingularVertex),3),5,'red');
    e1 = sum(data.edges(:,1)==find(TMesh.isSingularVertex),2)~=0;
    e2 = sum(data.edges(:,2)==find(TMesh.isSingularVertex),2)~=0;
    data.isSingularEdge = e1 & e2;
    data.isSingularVertex = TMesh.isSingularVertex;
    SEdges = find(data.isSingularEdge);
    
    %{
    hold on;
    for i = 1:size(SEdges)
        pts = data.vertices(data.edges(SEdges(i),:),:);
        plot3(pts(:,1),pts(:,2),pts(:,3))
    end
    %}
else
    %[X,T]=paul_loadTetGenMesh('./meshes/sphere_61k');
    %[X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\sphere\spherer.1');
    %[X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\sphere\spherer.1');
    [X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\anchor0\anchor_0.1');
    
    data = paul_getTetData(T,X,0);
    
    SEdges = [-1];
    
    %accumCurveEdges={};
    %[curveEdges, curveV, f] = chooseCurve(data, 0, 0, []); accumCurveEdges{1} = curveEdges;
    %SEdges = [curveEdges];
    
    %[curveEdges2, curveV2, f] = chooseCurve(data, f, 'b', accumCurveEdges); accumCurveEdges{2} = curveEdges2;
    %SEdges = [curveEdges; curveEdges2];
end

m = data.nonBoundaryTrianglesToTets;
g = graph(m(:,1),m(:,2),(1:size(data.nonBoundaryTrianglesToTets,1))');
%g.Edges.Weights =  randi(100,size(data.nonBoundaryTrianglesToTets,1),1);
%g.Edges.Labels = (1:size(data.nonBoundaryTrianglesToTets,1))';
dualspantree = minspantree(g);
nonBTrisInd = find(~data.isBoundaryTriangle);
TriInds = nonBTrisInd(dualspantree.Edges.Weight);
inTreeIndicator = sparse(TriInds, ones(size(TriInds,1),1), ones(size(TriInds,1),1), data.numTriangles, 1);

%dualspantree.Edges.Labels(35,:)
%dualspantree.Edges.EndNodes(35,:)
%data.nonBoundaryTrianglesToTets(dualspantree.Edges.Labels(35),:)
%data.trianglesToTets{nonBTrisInd(dualspantree.Edges.Labels(35))}'

%{
hold on;
dspv2v = dualspantree.Edges.EndNodes;
for i = 1:size(dspv2v)
    tets = dspv2v(i,:);
    tetps = data.tetBarycenters(tets,:);
    plot3(tetps(:,1),tetps(:,2),tetps(:,3),'r','linewidth',2);
    pause(.01)
end
%}

%{
hold on;
dspv2v = find(inTreeIndicator);
for i = 1:size(dspv2v)
    tri = dspv2v(i,:);
    tets = data.trianglesToTets{tri};
    assert(numel(tets)==2);
    tetps = data.tetBarycenters(tets,:);
    plot3(tetps(:,1),tetps(:,2),tetps(:,3),'r','linewidth',2);
    pause(.01)
end
%}

e2tIndicator = data.edgesToTrianglesIndicator ~= 0;
totalTrianglesPerEdge = sum(e2tIndicator ,2);

trianglesToClose = -1;
while(numel(trianglesToClose)~=0)
    closedTrianglesPerEdge = e2tIndicator*inTreeIndicator;
    DualLoopsToClose = find((totalTrianglesPerEdge - closedTrianglesPerEdge)==1);
    DualLoopsToClose = ARemoveB(DualLoopsToClose', SEdges')';
    trianglesToClose = find(sum(e2tIndicator(DualLoopsToClose, :),1)~=0);
    
    triangsintree = sum(inTreeIndicator);
    inTreeIndicator(trianglesToClose) = 1;
    
    triangsaddedToTree = sum(inTreeIndicator)-triangsintree;
    loopsVtriangs = numel(DualLoopsToClose)-triangsaddedToTree;
    %1
    if(loopsVtriangs~=0)
        'fewer triangles added to tree than loops closed!! possible contradictions.'
        %pause;
    end
end

% plot singular edges
if(SEdges~=-1)
    f = VisualizeEdges(SEdges, data, '-', 0, [1 0 0]);
else; f = figure; hold on; end;
% boundary verts
scatter3(data.vertices(find(data.isBoundaryVertex),1),data.vertices(find(data.isBoundaryVertex),2),data.vertices(find(data.isBoundaryVertex),3),1,'blue');
% plot nonmanifold edges
nonMEdges = find((totalTrianglesPerEdge - closedTrianglesPerEdge)>2);
f = VisualizeEdges(nonMEdges, data, '-',f,[1 0 1]);
% plot surface intersect with boundary
surfaceTriInds = find(~inTreeIndicator);
surfaceTriInds = ARemoveB(surfaceTriInds',find(data.isBoundaryTriangle)');
surfaceEdges = data.trianglesToEdges(surfaceTriInds, :);
surfaceEdges = intersect(surfaceEdges(:)', find(data.isBoundaryEdge)');
f = VisualizeEdges(unique(surfaceEdges), data, '-',f,[0 0 1]);

% plot primal slice surface
hold on;
%scatter3(data.vertices(find(data.isBoundaryVertex),1),data.vertices(find(data.isBoundaryVertex),2),data.vertices(find(data.isBoundaryVertex),3),5,'blue')
surfaceTriInds = find(~inTreeIndicator);
surfaceTriInds = ARemoveB(surfaceTriInds',find(data.isBoundaryTriangle)');
for i = 1:numel(surfaceTriInds)
    vinds = data.triangles(surfaceTriInds(i),:);
    polyPtch = data.vertices(vinds, :);
    ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .8);
    pause(.01)
end


