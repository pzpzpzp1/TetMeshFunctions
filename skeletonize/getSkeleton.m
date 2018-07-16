close all; clear all;

addpath('skeletonize')
path2HMeshFuncs='./../HexMeshSingularitySheetingTest';
addpath(path2HMeshFuncs);

%HMesh = LoadVTK(1, [path2HMeshFuncs '/meshes/FF/ellipsoid-B.vtk']);
HMesh = LoadVTK(1, [path2HMeshFuncs '/meshes/HTri2Hex/mesh01.vtk']);
 
rotate3d on

MetaVertices = HMesh.MetaVertices;
TMesh = HexToTet(HMesh);
X = TMesh.V2P;
T = TMesh.T2V;
[X2, T2] = preprocess_data(X, T);
data = paul_getTetData(T2,X2,1);

% subdivide at edges
[found, SEdges] = FindEdgeIndexOfVertPairs(TMesh.SE2V, data.edges);
SEdges = SEdges(~data.isBoundaryEdge(SEdges));
adjTetsPerEdge = (sum(data.tetsToEdgesIndicator(:,SEdges),2));
assert(sum(adjTetsPerEdge > 1)==0); % no tet can have more than one singular edge.
[X2, T2] = SubdivideTMeshAtTets(X2,T2,find(adjTetsPerEdge>0));
data = paul_getTetData(T2,X2,1);

% subdivide at vertices
[found, SEdges] = FindEdgeIndexOfVertPairs(TMesh.SE2V, data.edges);
SEdges = SEdges(~data.isBoundaryEdge(SEdges));
sverts = unique(data.edges(SEdges,:));
adjTetsPerVert = sum(data.tetsToVertsIndicator(:,sverts),2);
[X2, T2] = SubdivideTMeshAtTets(X2,T2,find(adjTetsPerVert>0));
data = paul_getTetData(T2,X2,1);

% subdivide at vertices
[found, SEdges] = FindEdgeIndexOfVertPairs(TMesh.SE2V, data.edges);
SEdges = SEdges(~data.isBoundaryEdge(SEdges));
sverts = unique(data.edges(SEdges,:));
coredTets = find(sum(data.tetsToVertsIndicator(:,sverts),2));
for i = coredTets'
    tris = data.tetsToTriangles(i,:);
    for tri = tris
        polyPtch = data.vertices(data.triangles(tri,:), :);
        ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); %alpha(ptc, .2);
    end
%     edges = data.tetsToEdges(i,:);
%     for e = edges
%         verts = data.vertices(data.edges(e,:),:);
%         plot3(verts(:,1),verts(:,2),verts(:,3),'r');
%     end
end

figure; rotate3d on; hold on; axis equal; title('remaining After coring line');
remainingTets = ARemoveB([1:data.numTetrahedra], coredTets);
for i = remainingTets
    tris = data.tetsToTriangles(i,:);
    for tri = tris
        polyPtch = data.vertices(data.triangles(tri,:), :);
        ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); %alpha(ptc, .2);
    end
end 
        
        