basepath = [userpath '\..\'];
addpath(genpath([basepath 'jsolomon']));

%filename = 'jsolomon\octahedral_frames\meshes\torus\torus1k.1';
%filename = 'jsolomon\octahedral_frames\meshes\moomoo\moomoo.1';
filename = 'jsolomon\octahedral_frames\meshes\torus\torus_39k';
%filename = 'jsolomon\octahedral_frames\meshes\elk\elk18k.1';
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
f = figure; hold on; axis equal; donechoosing  = false;
while ~donechoosing 
    hold off; scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),.01,'b'); hold on;
    nonboundaryedges = data.edges(find(~data.isBoundaryEdge),:);
    Adj = sparse(nonboundaryedges(:,1),nonboundaryedges(:,2),1:(size(nonboundaryedges,1)));
    startedge = randi(size(nonboundaryedges,1));
    v = data.vertices(nonboundaryedges(startedge,:)',:);
    endpoints = nonboundaryedges(startedge,:);
    plot3(v(:,1),v(:,2),v(:,3),'r'); desdir = v(2,:)-v(1,:); desdir=desdir/norm(desdir);
    notdone = any(~data.isBoundaryVertex(endpoints));
    curve = startedge;
    curveV = nonboundaryedges(startedge,:);
    while notdone
        for endpointInd = 1:2    
            if(~data.isBoundaryVertex(endpoints(endpointInd)))
                e1 = endpoints(endpointInd);

                candidateNv1 = find(Adj(e1,:));
                candidateNv2 = find(Adj(:,e1));

                candidateEdges1 = Adj(e1,candidateNv1);
                candidateEdges2 = Adj(candidateNv2,e1);

                candidateNv = [candidateNv1 candidateNv2'];
                candidateEdges = [candidateEdges1 candidateEdges2'];

                RemoveV = [curveV];
                removeIndicator = (sum(RemoveV == candidateNv',2) > 0)';
                candidateNv(removeIndicator ) = [];
                candidateEdges(removeIndicator ) = [];
                
                choice = randi(numel(candidateNv));
                dists = data.vertices(candidateNv',:) - repmat(data.vertices(e1,:),numel(candidateNv),1);
                dists = dists./sqrt(sum(dists.*dists,2))
                factor = (endpointInd-1)*2-1;
                [a b] = max(factor * (sum(dists.*repmat(desdir, numel(candidateNv),1),2)));
                choice = b;
                
                curveV = [curveV candidateNv(choice)];
                curve = [curve candidateEdges(choice)];

                endpoints(endpointInd)=candidateNv(choice);
                v = data.vertices([e1 candidateNv(choice)]',:);
                plot3(v(:,1),v(:,2),v(:,3),'r');
                %pause(1);
            end
        end

       notdone = any(~data.isBoundaryVertex(endpoints));
    end
    scatter3(data.vertices(endpoints',1),data.vertices(endpoints',2),data.vertices(endpoints',3),1);
    assert(all(data.isBoundaryVertex(endpoints)));
    donechoosing = input('Is this curve good enough? Answer 0 or 1. ');
end
% *** note curve above is indexed into nonboundaryedges not data.edges.
inds = find(~data.isBoundaryEdge); curveEdges = inds(curve); % reindex curve to data.edges.
fprintf("Curve accepted \n");

%% Try to close dual edges until no more can be closed
method = -1;
dualEdges = 1:data.numTriangles;
% get dual edges adjacent to the boundary. i.e. triangles with any boundary
% vert. Or dual edges adjacent to the curve.
dualEdgesCloseToBoundary = []; %find(sum(data.isBoundaryVertex(data.triangles),2)~=0);
isInCurve = sparse(curveV, ones(numel(curveV),1), ones(numel(curveV),1)); isInCurve(data.numVertices+1)=0;
dualEdgesCloseToCurve = find(sum(isInCurve(data.triangles),2)~=0);
dualEdgesToKeepOpen = unique([dualEdgesCloseToBoundary; dualEdgesCloseToCurve])';
dualEdgesToClose = dualEdges; dualEdgesToClose(dualEdgesToKeepOpen)=[];

% find the dual edges that are one edge away from closing.
GrowingTree = data.PrimalVolumeVertexSpanningTree;
GrowingTree = data.BoundaryLessPrimalSpanningTreeRelToEdges';
GrowingTree = unique([GrowingTree find(data.isBoundaryEdge)']);

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
DualSurfaceBoundaries = find(sum(inTreeIndicator(data.trianglesToEdges),2)==2);

%% assert triangles are all adj to boundary or curveV
for i = 1:numel(DualSurfaceBoundaries)
    tri = data.triangles(DualSurfaceBoundaries(i),:);
    inter = intersect([find(data.isBoundaryVertex) curveV], tri);
    if(numel(inter) == 0)
        DualSurfaceBoundaries(i) % triangle
        assert(0==1);
    end
end

size(GrowingTree)
%% Sanity check for a correct algorithm: this should be empty
indexInDualEdgesToCloseForNonManifoldTriangles = find(sum(inTreeIndicator(data.trianglesToEdges),2)==0);
size(indexInDualEdgesToCloseForNonManifoldTriangles)

%% Visualize boundary edges of the dual surface (that dont hit boundary of volume)
DualSurfaceBoundaries = find(sum(inTreeIndicator(data.trianglesToEdges),2)==2);
DualSurfaceBoundaries = ARemoveB(DualSurfaceBoundaries', find(data.isBoundaryTriangle)');
m = reshape(cell2mat(data.trianglesToTets(DualSurfaceBoundaries)),2,numel(DualSurfaceBoundaries));
p1 = data.tetBarycenters(m(1,:)',:);
p2 = data.tetBarycenters(m(2,:)',:);
interp = [0:.1:1]; final = [];
for i = 1:numel(interp)
    inter = interp(i) * p1 + p2 * (1-interp(i)); final = [final; inter];
end
f = figure; hold on; axis equal; 
scatter3(final(:,1),final(:,2),final(:,3),.1);
boundaryverts = data.vertices(find(data.isBoundaryVertex),:);
scatter3(boundaryverts(:,1),boundaryverts(:,2),boundaryverts(:,3),1);
VisualizeEdges(curveEdges, data, '-', f);

%% Visualize remaining surface
EdgesRemaining = 1:data.numEdges;
EdgesRemaining(GrowingTree)=[];
f = VisualizeDualEdges(data, EdgesRemaining, 'p');
VisualizeEdges(curveEdges, data, '-',f);

%% Simple Visualize remaining surface
EdgesRemaining = 1:data.numEdges;
EdgesRemaining(GrowingTree)=[];
f=VisualizeDualEdges(data, EdgesRemaining, '.');
VisualizeEdges(curveEdges, data, '.',f);













