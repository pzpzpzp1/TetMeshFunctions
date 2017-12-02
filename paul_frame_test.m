addpath(genpath('../../jsolomon'));

filename = 'C:\Users\Administrator\Documents\jsolomon\meshes\torus\torus1k.1';
filename = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\torus\torus1k.1';
filename = 'C:\Users\pzpzp\Documents\jsolomon\octahedral_frames\meshes\moomoo\moomoo.1';
[X,T]=loadTetGenMesh(filename)
mesh = getMeshData(X,T);
data = paul_getTetData(T,X);

%% initialize random frames
% TODO: set boundary tet frames to be appropriate
ypr = rand(data.numTetrahedra,3)*2*pi;
data.ypr = ypr;
q = angle2quat(ypr(:,1), ypr(:,2), ypr(:,3));
data.quatFrames = q;

%% compute transitions per triangleface
%t12 = quatmultiply(q2,quatinv(q1));
%t12*q1==q2
% fakeboundary is a copy of trianglesToTets, but with boundary tets
% replaced by tet 1.
%fakeboundary = data.trianglesToTets; fakeboundary(find(fakeboundary==-1))=1;
%q1 = q(fakeboundary(1,:),:);
%q2 = q(fakeboundary(2,:),:);
%t12 = paul_quatmultiply(q1,quatinv(q2)');
%t12 = quatmultiply(q1,quatinv(q2));
%data.faceTransitions = t12;
%data.tet2tetTransitions = paul_quatmultiply(q,quatinv(q)');
%qa = data.tet2tetTransitions(:,:,1);
%qb = data.tet2tetTransitions(:,:,2);
%qc = data.tet2tetTransitions(:,:,3);
%qd = data.tet2tetTransitions(:,:,4);

tetIndicesPerInteriorTriangle = [data.trianglesToTets{~data.isBoundaryTriangle}];
q1 = q(tetIndicesPerInteriorTriangle(1,:),:);
q2 = q(tetIndicesPerInteriorTriangle(2,:),:);
t12 = quatmultiply(q1,quatinv(q2));

data.faceTransitions = zeros(data.numTriangles, 4);
data.faceTransitions(find(~data.isBoundaryTriangle),:)=t12;
data.faceTransitions(find(data.isBoundaryTriangle),1)=1; % angle2quat(0,0,0) is 1,0,0,0.
data.tet2tetTransitions = paul_quatmultiply(q,quatinv(q)');
qa = data.tet2tetTransitions(:,:,1);
qb = data.tet2tetTransitions(:,:,2);
qc = data.tet2tetTransitions(:,:,3);
qd = data.tet2tetTransitions(:,:,4);

% TODO: randomly add elements of O to transitions

%% calculate singularity type of non boundary edges.
singularityTypes = zeros(data.numEdges,4);
nonboundaryEdgeIndices = find(~data.isBoundaryEdge);
for nonboundaryEdgeIndicesIndex = 1:numel(nonboundaryEdgeIndices)
    edgeIter = nonboundaryEdgeIndices(nonboundaryEdgeIndicesIndex);
    cycle = data.edgeCycles{edgeIter};
    % indexing into tet2tetTransitions has flipped index. for qi -> tij -> qj,
    % tij is stored at index (j,i);
    % tt=1;permute(data.tet2tetTransitions(cycle(tt+1),cycle(tt),:),[3,1,2])-[qa(index(tt)) qb(index(tt)) qc(index(tt)) qd(index(tt))]';
    % ^ is always 0 since index is correct.
    % also, matlab's full index is vertical first then horizontal.
    index = full((cycle(2:end)) +(cycle(1:end-1)-1).*data.numTetrahedra);
    
    % [t1 t2 t3 t4 ...] to apply all of them we want tn...t4*t3*t2*t1
    cycleTrans = [qa(index); qb(index); qc(index); qd(index)];
    singularityType = angle2quat(0,0,0);
    for transiter = 1:size(cycleTrans,2)
        singularityType = quatmultiply(cycleTrans(:,transiter)', singularityType);
    end
    singularityTypes(edgeIter,:) = singularityType;
end

%% TODO: calculate singularity type of boundary edges.


%z = [data.edgeCycles{~data.isBoundaryEdge}];




%{
%data.edgeCycles is edges x tets(e)
matches = ((data.trianglesToEdges(:,1)==find(data.isBoundaryEdge)') | ...
(data.trianglesToEdges(:,2)==find(data.isBoundaryEdge)') | ...
(data.trianglesToEdges(:,3)==find(data.isBoundaryEdge)'));
fdegree = sum(matches);
[i1, i2] = find(matches);
matches = matches.*[1:data.numTriangles]';
nonzeros = matches(find(matches~=0));
boundaryEdgeToTriangles = mat2cell(nonzeros,fdegree);
data.boundaryEdgeToTriangles = boundaryEdgeToTriangles;
%data.trianglesToEdges(boundaryEdgeToTriangles{N},:) will be a list of
%triangles that all share one edge. it wont be edge N because boundary
%indexing is different from all triangle indexing.
%}

