addpath(genpath('../../jsolomon'));

filename = 'C:\Users\Administrator\Documents\jsolomon\meshes\torus\torus1k.1';
filename = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\torus\torus1k.1';
filename = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\moomoo\moomoo.1';
filename = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\torus\torus_39k';
filename = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\elk\elk18k.1';
elkFrame = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\comparison_data\FF_ray_sokolov\FFfull\elk_18k_frame.txt';
%filename = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\comparison_data\FF_ray_sokolov\FFfull\elk_18k';
[X,T]=loadTetGenMesh(filename);
data = paul_getTetData(T,X);

VisualizeGraph(data.edges, data.vertices, 'r', .00001, data.PrimalVolumeVertexSpanningTree);
VisualizeGraph(data.nonBoundaryTrianglesToTets, data.tetBarycenters, 'b', .00001, data.DualVolumeVertexSpanningTree);
scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),1);

%% load elk frame
if(false)
    error('This part is turned off. Not implemented to load frames into mesh.');
    assert(strcmp(filename,'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\elk\elk18k.1'));
    fid = fopen(elkFrame,'r');
    frames = fscanf(fid,'%f',9*data.numTetrahedra);
    orthogMatrices = reshape(frames,3,3*data.numTetrahedra);
end

%% create test singular primal curve
hold on; axis equal;
Adj = sparse(nonboundaryedges(:,1),nonboundaryedges(:,2),1:(size(nonboundaryedges,1)));
scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),1,'b');
nonboundaryedges = data.edges(find(~data.isBoundaryEdge),:);
startedge = randi(size(nonboundaryedges,1));
v = data.vertices(nonboundaryedges(startedge,:)',:);
endpoints = nonboundaryedges(startedge,:);
plot3(v(:,1),v(:,2),v(:,3),'r');
notdone = any(data.isBoundaryVertex(endpoints));
curve = startedge;
curveV = nonboundaryedges(startedge,:);
while notdone
    
    if(~data.isBoundaryVertex(endpoints(1)))
        e1 = endpoints(1);
        
        candidateNv1 = find(Adj(e1,:));
        candidateNv2 = find(Adj(:,e1));
        
        candidateEdges1 = Adj(e1,:);
        candidateEdges2 = Adj(:,e1);
        
        candidateNv = [candidateNv1 candidateNv1'];
        candidateEdges = [candidateEdges1 candidateEdges2'];
    end
    
    if(~data.isBoundaryVertex(endpoints(2)))
        e1 = endpoints(2);
        
        candidateNv1 = find(Adj(e1,:));
        candidateNv2 = find(Adj(:,e1));
        
        candidateEdges1 = Adj(e1,:);
        candidateEdges2 = Adj(:,e1);
        
        candidateNv = [candidateNv1 candidateNv1'];
        candidateEdges = [candidateEdges1 candidateEdges2'];
        
        RemoveV = [curveV];
        removeIndicator = (sum(RemoveV == candidateNv',2) > 0)';
        candidateNv(removeIndicator ) = [];
        candidateEdges(removeIndicator ) = [];
        
    end
    
   
   notdone = any(data.isBoundaryVertex(endpoints));
end


















