basepath = [userpath '\..\'];
addpath(genpath([basepath 'jsolomon']));

%filename = 'jsolomon\octahedral_frames\meshes\torus\torus1k.1';
%filename = 'jsolomon\octahedral_frames\meshes\moomoo\moomoo.1';
%filename = 'jsolomon\octahedral_frames\meshes\torus\torus_39k';
%filename = 'jsolomon\octahedral_frames\meshes\elk\elk18k.1';
%filename = 'jsolomon\octahedral_frames\meshes\cylinder\cylinder_18k.1';
%filename = 'jsolomon\octahedral_frames\meshes\cube\cube_tri.1';
%filename = 'jsolomon\octahedral_frames\meshes\sphere\spherer.1';
filename = 'jsolomon\octahedral_frames\meshes\sphere\sphere_61k';

%elkFrame = 'jsolomon\octahedral_frames\comparison_data\FF_ray_sokolov\FFfull\elk_18k_frame.txt';
%filename = 'jsolomon\octahedral_frames\comparison_data\FF_ray_sokolov\FFfull\elk_18k';
[X,T]=loadTetGenMesh([basepath filename]);
data = paul_getTetData(T,X);

%{
VisualizeGraph(data.edges, data.vertices, 'r', .00001, data.PrimalVolumeVertexSpanningTree);
VisualizeGraph(data.nonBoundaryTrianglesToTets, data.tetBarycenters, 'b', .00001, data.DualVolumeVertexSpanningTree);
scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),1);
%}

%% load elk frame
if(false)
    error('This part is turned off. Not implemented to load frames into mesh.');
    assert(strcmp(filename,'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\elk\elk18k.1'));
    fid = fopen(elkFrame,'r');
    frames = fscanf(fid,'%f',9*data.numTetrahedra);
    orthogMatrices = reshape(frames,3,3*data.numTetrahedra);
end

%% create test singular primal curve
accumCurveEdges={};
[curveEdges, curveV, f] = chooseCurve(data, 0, 0, []); accumCurveEdges{1} = curveEdges;
[curveEdges2, curveV2, f] = chooseCurve(data, f, 'b', accumCurveEdges); accumCurveEdges{2} = curveEdges2;
%[curveEdges3, curveV3, f] = chooseCurve(data, f, 'b', accumCurveEdges); accumCurveEdges{3} = curveEdges3;
close;

%{ 
% just a peice to use for display
f = figure(); hold off; scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),.01,'b'); hold on;
for i = 1:numel(accumCurveEdges)
    VisualizeEdges(accumCurveEdges{i}, data, '.', f)
end
%}
%% compute tree without boundary, and without edges adj to curve, and with curve.
v2v = sparse(data.edges(:,1),data.edges(:,2),1:size(data.edges,1)); v2v(data.numVertices+1,data.numVertices+1)=0;
es = [v2v(curveV,:) v2v(:,curveV)'];
edgesAdjToCurve = full(sort(unique(es(find(es~=0)))));
exclude = unique([find(data.isBoundaryEdge); edgesAdjToCurve]);
exclude = edgesAdjToCurve; % include the boundary in the tree. empirically, decreases non manifold areas at boundary.
include = curveEdges2;
vinds = 1:data.numVertices; vinds(curveV)=[];
startVert = vinds(1);
%data.BoundaryLessPrimalSpanningTreeRelToEdges = PrimalVolumeVertexSpanningTreeMoreOps(data.edges, exclude, [], startVert);
data.BoundaryLessPrimalSpanningTreeRelToEdges = PrimalVolumeVertexSpanningTreeMoreOps(data.edges, exclude, include, []);
assert(numel(intersect(data.BoundaryLessPrimalSpanningTreeRelToEdges, curveEdges))==0); % by construction, curve shouldn't be in the tree yet.
data.BoundaryLessPrimalSpanningTreeRelToEdges = [data.BoundaryLessPrimalSpanningTreeRelToEdges curveEdges'];
% VisualizeGraph(data.edges, data.vertices, 'r', 1, data.BoundaryLessPrimalSpanningTreeRelToEdges)    

% add in 2nd curve to avoid it. this method fails. adding after the fact is
% unreliable and can result in no solution found even when there is one.
% data.BoundaryLessPrimalSpanningTreeRelToEdges = [data.BoundaryLessPrimalSpanningTreeRelToEdges curveEdges2']; 




%% Try to close dual edges until no more can be closed
dualEdges = 1:data.numTriangles;
% get dual edges adjacent to the boundary. i.e. triangles with any boundary
% vert. Or dual edges adjacent to the curve.
dualEdgesCloseToBoundary = []; %find(sum(data.isBoundaryVertex(data.triangles),2)~=0);
isInCurve = sparse(curveV, ones(numel(curveV),1), ones(numel(curveV),1)); isInCurve(data.numVertices+1)=0;
dualEdgesCloseToCurve = find(sum(isInCurve(data.triangles),2)~=0);
dualEdgesToKeepOpen = unique([dualEdgesCloseToBoundary; dualEdgesCloseToCurve])';
dualEdgesToClose = dualEdges; dualEdgesToClose(dualEdgesToKeepOpen)=[];

% find the dual edges that are one edge away from closing.
%GrowingTree = data.PrimalVolumeVertexSpanningTree;
GrowingTree = data.BoundaryLessPrimalSpanningTreeRelToEdges;
%GrowingTree = unique([GrowingTree find(data.isBoundaryEdge)' curveEdges']);
%GrowingTree = unique([GrowingTree curveEdges']);

inTreeIndicator = sparse(GrowingTree(:), ones(numel(GrowingTree),1),ones(numel(GrowingTree),1)); inTreeIndicator(data.numEdges+1)=0;
readyToClose = find(sum(inTreeIndicator(data.trianglesToEdges(dualEdgesToClose,:)),2)==2);
% readytoClose is an index into data.triangles(dualEdgesToClose,:));

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

%% Visualize remaining surface
EdgesRemaining = 1:data.numEdges;
EdgesRemaining(GrowingTree)=[];
f = VisualizeDualSurface(data, EdgesRemaining, 'p');
VisualizeEdges(curveEdges, data, '-',f);

%% Visualize non manifold locations: this should be very low. possibly empty depending on variations.
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
% primal surface. Method seems to leave none of these?
boundaryTrianglesWhoseDualEdgesAreBoundaries = find(sum(inTreeIndicator(data.trianglesToEdges(find(data.isBoundaryTriangle),:)),2)==2);
assert(numel(boundaryTrianglesWhoseDualEdgesAreBoundaries)==0);

% Sanity check. Number of edges adj to curve that arent the curve that end up in the
% tree. This is provably 0. consider any triangle adj to the curve. it has
% at least 2 edges in the edgesAdjToCurve set, which are not an will not
% ever be in the tree. so these edges will never get added to the tree.
assert(numel(ARemoveB(intersect(edgesAdjToCurve,GrowingTree)',curveEdges'))==0);

%% Simple Visualize remaining surface
%{
EdgesRemaining = 1:data.numEdges;
EdgesRemaining(GrowingTree)=[];
f=VisualizeDualSurface(data, EdgesRemaining, '.');
VisualizeEdges(curveEdges, data, '.',f);
%}












