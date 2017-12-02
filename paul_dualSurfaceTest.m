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
fprintf("Curve accepted \n");

%% Try to close dual edges until no more can be closed
dualEdges = 1:data.numTriangles;
% get dual edges adjacent to the boundary. i.e. triangles with any boundary
% vert. Or dual edges adjacent to the curve.
dualEdgesCloseToBoundary = find(sum(data.isBoundaryVertex(data.triangles),2)~=0);
isInCurve = sparse(curveV, ones(numel(curveV),1), ones(numel(curveV),1));
if(size(isInCurve,1)~=data.numVertices) isInCurve(data.numVertices)=0; end
dualEdgesCloseToCurve = find(sum(isInCurve(data.triangles),2)~=0);
%dualEdgesToClose = ARemoveB(dualEdges, unique([dualEdgesCloseToBoundary; dualEdgesCloseToCurve])');
dualEdgesToKeepOpen = unique([dualEdgesCloseToBoundary; dualEdgesCloseToCurve])';
dualEdgesToClose = dualEdges; dualEdgesToClose(dualEdgesToKeepOpen)=[];

% find the dual edges that are one edge away from closing.
GrowingTree = data.PrimalVolumeVertexSpanningTree;
inTreeIndicator = sparse(GrowingTree(:), ones(numel(GrowingTree),1),ones(numel(GrowingTree),1)); inTreeIndicator(data.numEdges+1)=0;
readyToClose = find(sum(inTreeIndicator(data.triangles(dualEdgesToClose,:)),2)==2);

%VisualizeGraph(data.edges, data.vertices, 'r', .00001, data.PrimalVolumeVertexSpanningTree);
while(numel(readyToClose) ~= 0)
    trisToClose = data.triangles(dualEdgesToClose(readyToClose),:);
    edgeIndsToCloseTris = find(~inTreeIndicator(data.triangles(dualEdgesToClose(readyToClose),:)));
    edgesToCloseTris = trisToClose(edgeIndsToCloseTris);
    
    GrowingTree = [GrowingTree edgesToCloseTris']; InTreeIndicator(edgesToCloseTris)=1;
    dualEdgesToClose = ARemoveB(dualEdgesToClose, dualEdgesToClose(readyToClose));
    readyToClose = find(sum(inTreeIndicator(data.triangles(dualEdgesToClose,:)),2)==2);
    %{
    for i = 1:numel(edgesToCloseTris)
        e = data.edges(edgesToCloseTris(i),:);
        v = data.vertices(e',:);
        plot3(v(:,1),v(:,2),v(:,3),'r');
        pause(1);
    end
    %}
end

indexInDualEdgesToCloseForNonManifoldTriangles = find(sum(inTreeIndicator(data.triangles(dualEdgesToClose,:)),2)==3);


%% Visualize remaining surface
EdgesRemaining = 1:data.numEdges;
EdgesRemaining(GrowingTree)=[];
f = figure; hold on; axis equal;
scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),1);
for e = 1:numel(EdgesRemaining)
    edgeInd = EdgesRemaining(e);
    %edge = data.edges(edgeInd,:);
    tets = data.edgeCycles{edgeInd};
    polyPtch = data.tetBarycenters(tets',:);
    ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .05);
    pause(.001);
end




