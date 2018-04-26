
close all; 
clear;
addpath('transitionManipulation');

Visualize = 0;
VisualizeMV = 0;
subdivideX = 1;
if(exist('cache/data.mat')==2)
    load('cache/data.mat');
    load('cache/SEdges.mat');
    load('cache/SEdgeValences.mat');
    load('cache/MetaVertices.mat');
    HMesh.MetaVertices = MetaVertices;
else

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
        HMesh = LoadVTK(1, [path2HMeshFuncs '/meshes/FF/ellipsoid-B.vtk']); rotate3d on
        
        %% simple meshes
        %HMesh = LoadVTK(1, [path2HMeshFuncs '/meshes/HTri2Hex/mesh01.vtk']);

        MetaVertices = HMesh.MetaVertices;
        TMesh = HexToTet(HMesh);
        X = TMesh.V2P;
        T = TMesh.T2V;
        [X2, T2] = preprocess_data(X, T);
        %% subdivide neighborhood of singular vertices to make resolution good enough for corner paths.
        for i = 1:numel(MetaVertices)
            ['subdividing ' num2str(MetaVertices{i}.vind)]
            
            [X2, T2] = SubdivideTMeshAtTets(X2, T2, find(sum(T2==MetaVertices{i}.vind,2)));
            for j = 1:subdivideX;
                mvert = MetaVertices{i};
                [X2, T2, BrokenSingularEdges] = SubdivideTMeshAtVertex(X2, T2, mvert.vind, 0);
                
                % update singular structures.
                adjSE2V = find(sum(TMesh.SE2V == MetaVertices{i}.vind,2));
                [found, ind] = FindEdgeIndexOfVertPairs(BrokenSingularEdges(:,[1 3]), TMesh.SE2V(adjSE2V,:)); assert(sum(found)==size(adjSE2V,1));
                BrokenSingularEdges = BrokenSingularEdges(found,:);
                ind = ind(found);
                
                [found, ind] = FindEdgeIndexOfVertPairs(BrokenSingularEdges(:,[1 3]), TMesh.SE2V); assert(all(found));
                TMesh.singularValence = [TMesh.singularValence; TMesh.singularValence(ind); TMesh.singularValence(ind)];
                TMesh.singularValence(ind)=[];
                TMesh.SE2V = [TMesh.SE2V; BrokenSingularEdges(:,[1 2]); BrokenSingularEdges(:,[2 3])];
                TMesh.SE2V(ind,:)=[];
                
                for c = 1:numel(mvert.corners)
                    [found, ind] = FindEdgeIndexOfVertPairs(mvert.corners{c}.threeEdges2V, BrokenSingularEdges(:,[1 3]));
                    assert(all(found));
                    mvert.corners{c}.threeEdges2V = BrokenSingularEdges(ind,[1 2]);
                end
                MetaVertices{i} = mvert;
            end
        end
        [X2,T2,data] = preprocess_data(X2,T2);
        

        %% load singular edges
    %     % scatter3(data.vertices(find(TMesh.isSingularVertex),1),data.vertices(find(TMesh.isSingularVertex),2),data.vertices(find(TMesh.isSingularVertex),3),5,'red');
    %     e1 = sum(data.edges(:,1)==find(TMesh.isSingularVertex),2)~=0;
    %     e2 = sum(data.edges(:,2)==find(TMesh.isSingularVertex),2)~=0;
    %     data.isSingularEdge = e1 & e2;
    %     data.isSingularVertex = TMesh.isSingularVertex;
    %     SEdges = find(data.isSingularEdge);

        %% match TMesh.SE2V with data.edges.
        [~, ia] = FindEdgeIndexOfVertPairs(TMesh.SE2V, data.edges);
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
        SEdgeValences = sparse(size(data.edges,1),1);
        MetaVertices = [];
        %accumCurveEdges={};
        %[curveEdges, curveV, f] = chooseCurve(data, 0, 0, []); accumCurveEdges{1} = curveEdges;
        %SEdges = [curveEdges];

        %[curveEdges2, curveV2, f] = chooseCurve(data, f, 'b', accumCurveEdges); accumCurveEdges{2} = curveEdges2;
        %SEdges = [curveEdges; curveEdges2];
    end
    
    save('cache/SEdges.mat','SEdges');
    save('cache/data.mat','data');
    save('cache/SEdgeValences.mat','SEdgeValences');
    save('cache/MetaVertices.mat','MetaVertices');
    
end

%% construct corners and sectors and loops.
% strategy: Dual graph is dual to a sphere triangulation. all dual loops 
% in neighborhood of a singular edge should
% excluded from being able to be used. After a path is found, that path
% should be removed from the dual graph.
PathsIndex = sparse(data.numEdges, data.numEdges); % Index to path of tet2tet.
Paths = {}; pathpos = 1;    
for i=1:numel(MetaVertices)
    mvert = MetaVertices{i};
    vind = mvert.vind;
    
    % initialize dgraph adjacency
    adjEdges = find(sum(data.edges==vind,2)~=0);
    sindst = find(sum(data.edges(SEdges,:)==vind,2)); adjSE2V = data.edges(SEdges(sindst),:);
    [found, seinds] = FindEdgeIndexOfVertPairs(adjSE2V, data.edges);
    assert(all(found));
    adjFaces = find(sum(data.triangles==vind,2)~=0);
    adjTets = find(sum(data.tetrahedra==vind,2)~=0);
    intAdjFaces = adjFaces(~data.isBoundaryTriangle(adjFaces));
    cageEdges = unique(ARemoveB(data.tetsToEdges(adjTets,:),adjEdges));
    cageVerts = unique(ARemoveB(data.tetrahedra(adjTets,:),vind));
    tet2tet = reshape(cell2mat(data.trianglesToTets(intAdjFaces)),2,[])';
    dAdj = sparse(tet2tet(:,1),tet2tet(:,2),ones(size(tet2tet,1),1),data.numTetrahedra+data.numEdges,data.numTetrahedra+data.numEdges); dAdj = dAdj + dAdj'; dAdj = dAdj ~= 0;
    
    
    %zero out neighborhood paths of adj singular edges
    %add connections to edge indicator.
    for seind = seinds'
        tets2Remove = data.edgesToTets{seind};
        dAdj(tets2Remove, tets2Remove)=0;
        % weight these as very expensive so paths try to avoid it.
        dAdj(data.numTetrahedra+seind,tets2Remove)=100000;
        dAdj(tets2Remove,data.numTetrahedra+seind)=100000;
    end
    
    for c = 1:numel(mvert.corners)
        corner = mvert.corners{c};
        
        [found, cornerEdgeInds] = FindEdgeIndexOfVertPairs(corner.threeEdges2V, data.edges);
        shiftedInds = circshift(cornerEdgeInds,1);
        for cei = 1:3
            estart = cornerEdgeInds(cei);
            eend = shiftedInds(cei);
            
            if(PathsIndex(estart,eend)~=0)
                continue;
            end
            
            dGraph = graph(dAdj);
            P = shortestpath(dGraph,data.numTetrahedra+estart,data.numTetrahedra+eend);
            
            if(numel(P)==0)
                'Path is empty. Tet resolution in this neighborhood isnt good enough!!'
                assert(false);
            end
            
            % make sure previous path can't be used again.
            dAdj(P(2:end-1),:)=0; dAdj(:,P(2:end-1))=0;
            
            PathsIndex(estart,eend) = pathpos; Paths{pathpos} = P; pathpos = pathpos + 1;
            PathsIndex(eend,estart) = pathpos; Paths{pathpos} = fliplr(P); pathpos = pathpos + 1;
        end
    end

    if(VisualizeMV)
        % Visualize local singular vertex.
        figure; axis equal; hold on; rotate3d on; title('Singular Vertex Neighborhood');
        scatter3(data.vertices(mvert.vind,1),data.vertices(mvert.vind,2),data.vertices(mvert.vind,3),50,'filled','k');
        scatter3(data.vertices(cageVerts,1),data.vertices(cageVerts,2),data.vertices(cageVerts,3),10,'filled','r');
%         for tet = adjTets'
%             for eind = data.tetsToEdges(tet,:);
%                 vs = data.vertices(data.edges(eind,:),:);
%                 plot3(vs(:,1),vs(:,2),vs(:,3),'k');
%             end
%         end
        for eind = cageEdges;
            vs = data.vertices(data.edges(eind,:),:);
            plot3(vs(:,1),vs(:,2),vs(:,3),'k','LineWidth',.1);
        end
        for seind = seinds'
            vs = data.vertices(data.edges(seind,:),:);
            plot3(vs(:,1),vs(:,2),vs(:,3),'b','LineWidth',2);
        end
%         for face = intersect(adjFaces,find(~isTrivialTriangle))'
%             polyPtch = data.vertices(data.triangles(face,:), :);
%             ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .2);
%         end
        for c = 1:numel(mvert.corners)
            corner = mvert.corners{c};
            adddir = (sum(data.vertices(corner.threeEdges2V(:,2),:)-data.vertices(corner.threeEdges2V(:,1),:))/3)*.1;
            for j=1:3
                vs = data.vertices(corner.threeEdges2V(j,:),:);
                vs = vs + [adddir;adddir];

                if(corner.threeEdgesDeg(j)==3)
                    plot3(vs(:,1),vs(:,2),vs(:,3),'b','LineWidth',2);
                else
                    plot3(vs(:,1),vs(:,2),vs(:,3),'g','LineWidth',2);
                end
            end
        end
        
        for i = 1:numel(Paths)
            p = Paths{i};
            ptets = p(2:end-1);
            
            ptet2tet = [ptets;circshift(ptets,-1)]; ptet2tet = ptet2tet(:,1:end-1)';
            for ppart = 1:size(ptet2tet,1)
                tbc = data.tetBarycenters(ptet2tet(ppart,:),:);
                plot3(tbc(:,1),tbc(:,2),tbc(:,3),'r','LineWidth',2); hold on;
            end
        end
    end
end

figure; rotate3d on;
%% Spanning tree
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
%     if(loopsVtriangs~=0)
%         % this is made to count how many dual loops where closed by adding
%         % triangles to the tree. Each dual loop represents constraints that
%         % had to be satisfied. So having more than 1 constrain on the value
%         % of one triangle means it's possible that no value of that
%         % triangle satisfies all constraints.
%         'fewer triangles added to tree than loops closed!! possible contradictions in the system?'
%         %pause;
%     end
end

triangleToTransition = cell(data.numTriangles,1);
triangleHasTransition = zeros(data.numTriangles,1);
triangleToTransition(find(inTreeIndicator))={[]};
triangleHasTransition(find(inTreeIndicator))=1;

surfaceTriInds = find(~inTreeIndicator);
surfaceTriInds = ARemoveB(surfaceTriInds',find(data.isBoundaryTriangle)');

nonMEdges = find((totalTrianglesPerEdge - closedTrianglesPerEdge)>2);

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
%     polyPtch = data.vertices(data.triangles(tri,:), :); ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .8);
        
    StartingTri{genpos}=tri;
    triangleHasTransition(tri)=1;
    triangleToTransition{tri}=genpos;
    
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
    trisAdded = -1;
    todel = {}; todelpos = 1;
    
    while(numel(trianglesToClose)~=0)
%         for tdel = 1:numel(todel); delete(todel{tdel}); end; todel = {}; todelpos = 1; if(trisAdded ~= -1) for t=trisAdded;             polyPtch = data.vertices(data.triangles(t,:), :); ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'black'); alpha(ptc, .3);        end; end;
        
        closedTrianglesPerEdge = e2tIndicator*inTreeIndicator;
        DualLoopsToClose = find((totalTrianglesPerEdge - closedTrianglesPerEdge)==1);
        DualLoopsToClose = ARemoveB(DualLoopsToClose', SEdges')';
        trianglesToClose = find(sum(e2tIndicator(DualLoopsToClose, :),1)~=0);
        triangsintree = sum(inTreeIndicator);
        newTrisAddedInds = find(~inTreeIndicator(trianglesToClose));
        inTreeIndicator(trianglesToClose(newTrisAddedInds)) = 1;
        MetaSurfaceClosed{genpos} = [MetaSurfaceClosed{genpos} trianglesToClose(newTrisAddedInds)];
        
        % indexes into data.triangles
        trisAdded = trianglesToClose(newTrisAddedInds);
        for tris = trisAdded(:)'
%             polyPtch = data.vertices(data.triangles(tris,:), :); ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'red'); alpha(ptc, .3); todel{todelpos}=ptc; todelpos = todelpos+1;
            
            edges = data.trianglesToEdges(tris,:);
            for edge = edges(:)'
                if(data.isBoundaryEdge(edge) || data.isSingularEdge(edge))
                    continue; 
                end
                triCycle = data.edgeTriCycles{edge}(1:end-1);
                if(sum(~triangleHasTransition(triCycle))==1)
                    % only one triangle is left unknown.
                    missingTriInd = find(~triangleHasTransition(triCycle));
                    missingTri = triCycle(missingTriInd);
                    triCycle = circshift(triCycle, numel(triCycle) - missingTriInd);
                    
                    accumTransitions = [];
                    for accumTris = triCycle(1:end-1)
                        % get transitions of a triangle.
                        transition = triangleToTransition{accumTris};
                        
                        side1 = min(find(data.edgeCycles{edge}==data.trianglesToTets{accumTris}(1)));
                        forward = data.edgeCycles{edge}(side1+1)==data.trianglesToTets{accumTris}(2);
                        
                        if(forward)
                            accumTransitions = [accumTransitions transition];
                        else
                            accumTransitions = [accumTransitions fliplr(-1*transition)];
                        end
                    end
                    
                    % orient solution into transition store.
                    side1 = min(find(data.edgeCycles{edge}==data.trianglesToTets{triCycle(end)}(1)));
                    forward = data.edgeCycles{edge}(side1+1)==data.trianglesToTets{triCycle(end)}(2);
                    if(~forward); accumTransitions = -fliplr(accumTransitions); end;
                    
                    triangleToTransition{triCycle(end)} = -fliplr(accumTransitions);
                    triangleHasTransition(triCycle(end))=1;
                    
                    assert(numel(cancelAntiPairs(transitionsPerEdge(edge,data,triangleToTransition)))==0)
                    
%                     if(~(numel(transitions)==1 || ~any(transitions == circshift(transitions,1))))
%                         'paused!'
%                         pause;
%                          alpha(ptc, .9);
%                         vs = data.vertices(data.edges(edge,:),:); hold on; plot3(vs(:,1),vs(:,2),vs(:,3),'g','LineWidth',2)
%                     end
                    
                    break;
                end
            end
        end
    end
    MetaSurfaceClosed{genpos} = unique(MetaSurfaceClosed{genpos});
    surfaceToPuncture = find(~inTreeIndicator);
    surfaceToPuncture = ARemoveB(surfaceToPuncture',find(data.isBoundaryTriangle)');
    genpos = genpos + 1; MetaSurfaceClosed{genpos} = [];
end
MetaSurfaceClosed=MetaSurfaceClosed(1:end-1);
fprintf(['Num Generators ' num2str(numel(H1Generators)) '\n']);

%% prune generators if they are not real d.o.f.
longstr = '';
for i = 1:numel(MetaSurfaceClosed)
    tris = MetaSurfaceClosed{i};
    longstr = [longstr '[' num2str(triangleToTransition{tris(1)})];
    for t = tris(2:end)
        longstr = [longstr ' ; ' num2str(triangleToTransition{t})];
    end
    longstr = [longstr '] \n'];
end
fprintf(['Showing transitions \n' longstr '\n']);

% simplify transitions and sanity checks
pruneTriangles = {};
assert(all(triangleHasTransition(find(~data.isBoundaryTriangle))));
for i = 1:numel(MetaSurfaceClosed)
    tris = MetaSurfaceClosed{i};
    pruneTriangles{i} = [];
    for tind = 1:numel(tris)
        t = tris(tind);
        transitions = triangleToTransition{t};
        
        % find consecutive flipped values and remove
        transitions = cancelAntiPairs(transitions);
        triangleToTransition{t} = transitions;
        
        if(numel(transitions)==0)
            pruneTriangles{i} = [pruneTriangles{i} tind];
        end
    end
end

for i = 1:numel(MetaSurfaceClosed)
    MetaSurfaceClosed{i}(pruneTriangles{i})=[];
end

isTrivialTriangle = zeros(data.numTriangles,1);
for i = 1:numel(triangleHasTransition)
    trans = triangleToTransition{i};
    isTrivialTriangle(i) = numel(trans)==0;
end


%% TODO:  prune generators based on regular holonomy?? currently fails. lots of false positives.
% regular non manifold edges
regularNMEdges = ARemoveB(ARemoveB(nonMEdges, SEdges),find(data.isBoundaryEdge));
generatorsToRemove = [];
for eiter = 1:numel(regularNMEdges)
    edgeind = regularNMEdges(eiter);
    transitions = transitionsPerEdge(edgeind,data,triangleToTransition);
    transitions = cancelAntiPairs(transitions);
    if(numel(transitions)~=0)
        generatorsToRemove = unique([generatorsToRemove transitions]);
    end
end

save('cache/SEdges.mat','SEdges');
save('cache/data.mat','data');
save('cache/surfaceTriInds.mat','surfaceTriInds');
save('cache/MetaVertices.mat','MetaVertices');
save('cache/SEdgeValences.mat','SEdgeValences');
save('cache/H1Generators.mat','H1Generators');
save('cache/MetaSurfaceClosed.mat','MetaSurfaceClosed');
save('cache/StartingTri.mat','StartingTri');
save('cache/H1DualEdgeGenerators.mat','H1DualEdgeGenerators');
save('cache/triangleToTransition.mat','triangleToTransition');
save('cache/PathsIndex.mat','PathsIndex');
save('cache/Paths.mat','Paths');

if(Visualize)
    %% Visualize stuff
    % plot singular edges
    if(numel(SEdges)~=0)
        f = VisualizeEdges(SEdges, data, '-', 0, [1 0 0]);
    else; f = figure; hold on; end;
    % boundary verts
    scatter3(data.vertices(find(data.isBoundaryVertex),1),data.vertices(find(data.isBoundaryVertex),2),data.vertices(find(data.isBoundaryVertex),3),1,'blue');
    % plot nonmanifold edges
    % not needed anymore. moved up. nonMEdges = find((totalTrianglesPerEdge - closedTrianglesPerEdge)>2);
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

