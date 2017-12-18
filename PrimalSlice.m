
%% run from TetMeshFunctions
path2HMeshFuncs='./../HexMeshSingularitySheetingTest';
addpath(path2HMeshFuncs);

% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/bunny_ours.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/bunny.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/rocker_arm.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/polycut2013/fertility.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/sculpture-B_ours.vtk']);

% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/rod_ours.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/fertility.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/double_torus.vtk']);
% HMesh = LoadVTK(0, [path2HMeshFuncs '/meshes/FF/ellipsoid-B.vtk']);

TMesh = HexToTet(HMesh);
X = TMesh.V2P;
T = TMesh.T2V;

%[X,T]=paul_loadTetGenMesh('./meshes/sphere_61k');
%[X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\sphere\spherer.1');
[X,T]=paul_loadTetGenMesh('C:\Users\Administrator\Documents\jsolomon\octahedral_frames\meshes\sphere\spherer.1');

%data = paul_getTetData(T,X,1);
data = paul_getTetData(T,X,0);

% data.nonBoundaryTriangles = 
% i = randi(size(data.nonBoundaryTrianglesToTets,1));
% data.nonBoundaryTrianglesToTets(i)
% data.trianglesToEdges


m = data.nonBoundaryTrianglesToTets;
g = graph(m(:,1),m(:,2));
g.Edges.Weights =  randi(100,size(data.nonBoundaryTrianglesToTets,1),1);
g.Edges.Labels = (1:size(data.nonBoundaryTrianglesToTets,1))';
dualspantree = minspantree(g);
nonBTrisInd = find(~data.isBoundaryTriangle);
TriInds = nonBTrisInd(dualspantree.Edges.Labels);
inTreeIndicator = sparse(TriInds, ones(size(TriInds,1),1), ones(size(TriInds,1),1), data.numTriangles, 1);

%data.nonBoundaryTrianglesToTets(dualspantree.Edges.Labels(5),:)
%data.trianglesToTets{nonBTrisInd(dualspantree.Edges.Labels(5))}

%{
hold on;
dspv2v = dualspantree.Edges.EndNodes;
for i = 1:size(dspv2v)
    tets = dspv2v(i,:);
    tetps = data.tetBarycenters(tets,:);
    plot3(tetps(:,1),tetps(:,2),tetps(:,3));
    pause(.01)
end
%}

e2tIndicator = data.edgesToTrianglesIndicator ~= 0;
totalTrianglesPerEdge = sum(e2tIndicator ,2);

trianglesToClose = -1;
while(numel(trianglesToClose)~=0)
    closedTrianglesPerEdge = sum(e2tIndicator.*inTreeIndicator',2);
    DualLoopsToClose = find((totalTrianglesPerEdge - closedTrianglesPerEdge)==1);
    trianglesToClose = find(sum(e2tIndicator(DualLoopsToClose, :),1)~=0);
    inTreeIndicator(trianglesToClose) = 1;
end

% triangInds = 1:data.numTriangles;
% intTriToTriIndex = triangInds(find(~data.isBoundaryTriangle));

%tbs = data.triangleBarycenters(find(~inTreeIndicator),:);
%scatter3(tbs(:,1),tbs(:,2),tbs(:,3),1);

hold on;
scatter3(data.vertices(find(data.isBoundaryVertex),1),data.vertices(find(data.isBoundaryVertex),2),data.vertices(find(data.isBoundaryVertex),3),1,'green')
surfaceTriInds = find(~inTreeIndicator);
for i = 1:numel(surfaceTriInds)
    vinds = data.triangles(surfaceTriInds(i),:);
    polyPtch = data.vertices(vinds, :);
    ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .05);
end



