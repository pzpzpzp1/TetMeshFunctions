Visualize = 0;


if(exist('cache/data.mat')==2)
    'WARNING: cache contains data. will be overwritten.'
    pause;
end

%1. start from hex mesh (with singularities) 
%2. start from tet mesh (no singularities)
flow = 1;
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
    %HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/double_torus.vtk']);
    HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/ellipsoid-B.vtk']);

    TMesh = HexToTet(HMesh);
    X = TMesh.V2P;
    T = TMesh.T2V;
    data = paul_getTetData(T,X,0);
    
    %% load singular edges
%     % scatter3(data.vertices(find(TMesh.isSingularVertex),1),data.vertices(find(TMesh.isSingularVertex),2),data.vertices(find(TMesh.isSingularVertex),3),5,'red');
%     e1 = sum(data.edges(:,1)==find(TMesh.isSingularVertex),2)~=0;
%     e2 = sum(data.edges(:,2)==find(TMesh.isSingularVertex),2)~=0;
%     data.isSingularEdge = e1 & e2;
%     data.isSingularVertex = TMesh.isSingularVertex;
%     SEdges = find(data.isSingularEdge);
    
    %% match TMesh.SE2V with data.edges.
    flipped = data.edges(:,1)>data.edges(:,2); data.edges(flipped,:) = [data.edges(flipped,2) data.edges(flipped,1)];
    ia = find(ismember(data.edges,TMesh.SE2V,'rows'));
    data.isSingularEdge = sparse(ia,ones(size(ia)),ones(size(ia)),size(data.edges,1),1);
    SEdgeValences = sparse(ia, ones(size(ia)), TMesh.singularValence, size(data.edges,1),1);
    SEdges = find(data.isSingularEdge);
    
%     hold on; axis equal;
%     for se = 1:size(SEdges,1)
%         v = data.vertices(data.edges(SEdges(se),:),:);
%         plot3(v(:,1),v(:,2),v(:,3),'r')
%     end
%     
%     hold on; axis equal;
%     for se = 1:size(TMesh.SE2V,1)
%         v = TMesh.V2P(TMesh.SE2V(se,:),:);
%         plot3(v(:,1),v(:,2),v(:,3),'r')
%     end
    
else
    %[X,T]=paul_loadTetGenMesh('./meshes/sphere_61k');
    %[X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\sphere\spherer.1');
    %[X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\sphere\spherer.1');
    [X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\anchor0\anchor_0.1');
    
    data = paul_getTetData(T,X,0);
    
    SEdges = [];
    
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
        % this is made to count how many dual loops where closed by adding
        % triangles to the tree. Each dual loop represents constraints that
        % had to be satisfied. So having more than 1 constrain on the value
        % of one triangle means it's possible that no value of that
        % triangle satisfies all constraints.
        'fewer triangles added to tree than loops closed!! possible contradictions in the system?'
        %pause;
    end
end

surfaceTriInds = find(~inTreeIndicator);
surfaceTriInds = ARemoveB(surfaceTriInds',find(data.isBoundaryTriangle)');

%% get H1 generators.
%dualspantree = minspantree(g);
%nonBTrisInd = find(~data.isBoundaryTriangle);
%TriInds = nonBTrisInd(dualspantree.Edges.Weight);
%inTreeIndicator = sparse(TriInds, ones(size(TriInds,1),1), ones(size(TriInds,1),1), data.numTriangles, 1);
root = randi(data.numTetrahedra);
surfaceToPuncture = surfaceTriInds;
%f = VisualizeEdges(SEdges, data, '-', 0, [1 0 0]);
H1DualEdgeGenerators = {};
H1Generators = {}; MetaSurfaceClosed{1} = []; StartingTri = {}; genpos = 1;
while(numel(surfaceToPuncture)~=0)
    tri = surfaceToPuncture(randi(numel(surfaceToPuncture)));
    StartingTri{genpos}=tri;
    tets = data.trianglesToTets{tri};
    p1 = shortestpath(dualspantree,root,tets(1));
    p2 = shortestpath(dualspantree,root,tets(2));
    % loop of tets
    p3 = [p1 fliplr(p2)];
    H1Generators{genpos}=p3; 
    % loop of triangles
    tetloop2triloopInds = repmat(2:numel(p3),2,1); tetloop2triloopInds = tetloop2triloopInds(:);
    tetloop2triloopInds = [1 tetloop2triloopInds(1:end-1)'];
    triloop2tets = reshape(p3(tetloop2triloopInds),2,[])'; 
    flipinds = triloop2tets(:,1) < triloop2tets(:,2); triloop2tets(flipinds,:) = [triloop2tets(flipinds,2) triloop2tets(flipinds,1)];
    flipinds = data.nonBoundaryTrianglesToTets(:,1) < data.nonBoundaryTrianglesToTets(:,2); data.nonBoundaryTrianglesToTets(flipinds,:) = [data.nonBoundaryTrianglesToTets(flipinds,2) data.nonBoundaryTrianglesToTets(flipinds,1)];
    [isit, matchind] = ismember(triloop2tets, data.nonBoundaryTrianglesToTets,'rows');
    H1DualEdgeGenerators{genpos}=matchind;
    %plot3(data.tetBarycenters(p3,1),data.tetBarycenters(p3,2),data.tetBarycenters(p3,3),'b');
    %% update surfaceToPuncture
    inTreeIndicator(tri)=1;
    MetaSurfaceClosed{genpos} = [MetaSurfaceClosed{genpos} tri];
    trianglesToClose = -1; 
    while(numel(trianglesToClose)~=0)
        closedTrianglesPerEdge = e2tIndicator*inTreeIndicator;
        DualLoopsToClose = find((totalTrianglesPerEdge - closedTrianglesPerEdge)==1);
        DualLoopsToClose = ARemoveB(DualLoopsToClose', SEdges')';
        trianglesToClose = find(sum(e2tIndicator(DualLoopsToClose, :),1)~=0);
        triangsintree = sum(inTreeIndicator);
        newTrisAddedInds = find(~inTreeIndicator(trianglesToClose));
        inTreeIndicator(trianglesToClose(newTrisAddedInds)) = 1;
        MetaSurfaceClosed{genpos} = [MetaSurfaceClosed{genpos} trianglesToClose(newTrisAddedInds)];
    end
    MetaSurfaceClosed{genpos} = unique(MetaSurfaceClosed{genpos});
    surfaceToPuncture = find(~inTreeIndicator);
    surfaceToPuncture = ARemoveB(surfaceToPuncture',find(data.isBoundaryTriangle)');
    genpos = genpos + 1; MetaSurfaceClosed{genpos} = [];
end
numel(H1Generators)



MetaVertices = HMesh.MetaVertices;
save('cache/SEdges.mat','SEdges');
save('cache/data.mat','data');
save('cache/surfaceTriInds.mat','surfaceTriInds');
save('cache/MetaVertices.mat','MetaVertices');
save('cache/SEdgeValences.mat','SEdgeValences');
save('cache/H1Generators.mat','H1Generators');
save('cache/MetaSurfaceClosed.mat','MetaSurfaceClosed');
save('cache/StartingTri.mat','StartingTri');
save('cache/H1DualEdgeGenerators.mat','H1DualEdgeGenerators');


if(Visualize)
    %% Visualize stuff
    % plot singular edges
    if(numel(SEdges)~=0)
        f = VisualizeEdges(SEdges, data, '-', 0, [1 0 0]);
    else; f = figure; hold on; end;
    % boundary verts
    scatter3(data.vertices(find(data.isBoundaryVertex),1),data.vertices(find(data.isBoundaryVertex),2),data.vertices(find(data.isBoundaryVertex),3),1,'blue');
    % plot nonmanifold edges
    nonMEdges = find((totalTrianglesPerEdge - closedTrianglesPerEdge)>2);
    f = VisualizeEdges(nonMEdges, data, '-',f,[1 0 1]);
    % plot surface intersect with boundary
    surfaceEdges = data.trianglesToEdges(surfaceTriInds, :);
    surfaceEdges = intersect(surfaceEdges(:)', find(data.isBoundaryEdge)');
    f = VisualizeEdges(unique(surfaceEdges), data, '-',f,[0 0 1]);

    % plot primal slice surface
    hold on;
    %scatter3(data.vertices(find(data.isBoundaryVertex),1),data.vertices(find(data.isBoundaryVertex),2),data.vertices(find(data.isBoundaryVertex),3),5,'blue')
    for i = 1:numel(surfaceTriInds)
        vinds = data.triangles(surfaceTriInds(i),:);
        polyPtch = data.vertices(vinds, :);
        ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .8);
        pause(.01)
    end
end

