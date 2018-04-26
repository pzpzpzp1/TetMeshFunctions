clear;
close all;

Visualize = true;
VerifyLoops = true;
VisualizeLoops = false;

%function z = SpanningTreeApproach
load('cache/SEdges.mat');
load('cache/data.mat');
load('cache/surfaceTriInds.mat');
load('cache/MetaVertices.mat');
load('cache/SEdgeValences.mat');
load('cache/H1Generators.mat');
load('cache/MetaSurfaceClosed.mat');
load('cache/StartingTri.mat');
load('cache/H1DualEdgeGenerators.mat');
load('cache/triangleToTransition.mat');
load('cache/Paths.mat');
load('cache/PathsIndex.mat');

Surf2Edges = data.trianglesToEdges(surfaceTriInds,:);
E2Fcount = accumarray(Surf2Edges(:),1); 
if(numel(E2Fcount) < data.numEdges); E2Fcount(data.numEdges)=0; end;

SurfaceEdge = find(data.isBoundaryEdge & (E2Fcount ~= 0));
InteriorEdge = find(~data.isBoundaryEdge & (E2Fcount ~= 0) & (E2Fcount ~= 2));
SingularEdges = SEdges;
SingularValences = SEdgeValences;
MetaVertices = MetaVertices;

allverts = data.edges(unique([SurfaceEdge; InteriorEdge; SingularEdges]),:);
vcount = accumarray(allverts(:),1,[data.numVertices 1]);

%% time to reduce the generator set. Some are redundant because of the shape of the msp.
% Q: how can there be more than 12 loops? 
% A: imagine a surface contacts itself, creating an overlap triangle. If
% the msp chooses to peirce that triangle, it creates a useless loop,
% because that triangle has to have identity holonomy. Peircing its source
% surface would have removed this extra triangle too, but in reverse that
% doesn't happen. We get an extra generator, with identity holonomy.
%% TODO: remove extra generators.

% NOTICE: given a dual loop, how do you express it in terms of generators?
% also notice, not necessary for algorithm to work. can use symbolic
% propagation
interiorSingularEdges = intersect(find(~data.isBoundaryEdge),SEdges);
triloop = data.edgesToTrianglesUnoriented{interiorSingularEdges(1)};
%gens = DualLoopToGenerators(data, triloop, H1DualEdgeGenerators, DualFacesToExclude);

%% list corner constraints in terms of symbolic values 
tetdoubleloopscornertype = {}; dloopind = 1;
tettripleloopsedgetype = {}; t1loopind = 1;
tettripleloopsfacetype = {}; t2loopind = 1;
for i = 1:numel(MetaVertices)
    mvert = MetaVertices{i};
    MetaVertices{i}.loops = {}; mvertloopind = 1;
    vind = mvert.vind;
    for c = 1:numel(mvert.corners)
        corner = mvert.corners{c};
        [found, cornerEdgeInds] = FindEdgeIndexOfVertPairs(corner.threeEdges2V, data.edges);
        shiftedInds = circshift(cornerEdgeInds,1);
        
        for cei = 1:3
            estart = cornerEdgeInds(cei);
            eend = shiftedInds(cei);
            P = Paths{PathsIndex(estart,eend)};
            %construct full double loops!
            if(all(~data.isBoundaryEdge([estart, eend])))
                tns = data.edgeCycles{estart}(1:end-1);
                tne = data.edgeCycles{eend}(1:end-1);
                
                eseeOriented = data.edges([estart eend],1)==vind;
                if(eseeOriented(1))
                    tns = fliplr(tns);
                end
                if(eseeOriented(2))
                    tne = fliplr(tne);
                end
                
                assert(any(P(2) == tns));
                assert(any(P(end-1) == tne));
                
                %construct full loops!
                tns = circshift(tns,-find(P(2) == tns,1)+1);
                tne = circshift(tne,-find(P(end-1) == tne,1)+1);
                assert(tns(1)==P(2));
                assert(tne(1)==P(end-1));
                
                Pes = P([1 end]);
                Pmid = P(2:end-1);
                
                Pmid = [tns Pmid tne(2:end) fliplr(Pmid)];
                P = [Pes(1) Pmid Pes(2)];
                assert(P(2)==P(end-1));
                tetdoubleloopscornertype{dloopind} = Pmid; dloopind = dloopind + 1;
                MetaVertices{i}.loops{mvertloopind} = Pmid; mvertloopind = mvertloopind + 1;
            end
        end
        
        %% construct triple loops
        if(all(~data.isBoundaryEdge(cornerEdgeInds)))
            e1 = cornerEdgeInds(1);
            e2 = cornerEdgeInds(2);
            e3 = cornerEdgeInds(3);
            
            tn1 = data.edgeCycles{e1}(1:end-1);
            tn2 = data.edgeCycles{e2}(1:end-1);
            tn3 = data.edgeCycles{e3}(1:end-1);

            eseeOriented = data.edges(cornerEdgeInds,1)==vind;
            if(eseeOriented(1))
                tn1 = fliplr(tn1);
            end
            if(eseeOriented(2))
                tn2 = fliplr(tn2);
            end
            if(eseeOriented(3))
                tn3 = fliplr(tn3);
            end
            
            %construct full triple loop!
            P12 = Paths{PathsIndex(e1,e2)}(2:end-1);
            P23 = Paths{PathsIndex(e2,e3)}(2:end-1);
            P31 = Paths{PathsIndex(e3,e1)}(2:end-1);
            
            tn2 = circshift(tn2,-find(P12(end) == tn2,1));
            tn2 = tn2(1:(find(P23(1)==tn2,1)-1));
            tn3 = circshift(tn3,-find(P23(end) == tn3,1));
            tn3 = tn3(1:(find(P31(1)==tn3,1)-1));
            tn1 = circshift(tn1,-find(P31(end) == tn1,1));
            tn1 = tn1(1:(find(P12(1)==tn1,1)));
            fullPath = [P12 tn2 P23 tn3 P31 tn1];
            
            % insert into the appropriate cell 
            % tettripleloopsedgetype t1loopind tettripleloopsfacetype t2loopind
            if(mod(sum(corner.threeEdgesDeg==5),2)==0)
                % edgetype
                tettripleloopsedgetype{t1loopind}=fullPath;
                t1loopind = t1loopind + 1;
                MetaVertices{i}.loops{mvertloopind} = fullPath; mvertloopind = mvertloopind + 1;
            else
                % facetype
                tettripleloopsfacetype{t2loopind}=fullPath;
                t2loopind = t2loopind + 1;
                MetaVertices{i}.loops{mvertloopind} = fullPath; mvertloopind = mvertloopind + 1;
            end
        end
    end
end

% verify loops
if(VerifyLoops)
    for i = 1:numel(tettripleloopsedgetype)
        tetloop = tettripleloopsedgetype{i};
        tetloopExt = [tetloop tetloop(1)];
        %figure; hold on; axis equal; rotate3d on;
        for j = 1:numel(tetloop)
            pairtets = tetloopExt(j:j+1);
            sharedTri = intersect(data.tetsToTriangles(pairtets(1),:),data.tetsToTriangles(pairtets(2),:));
            assert(numel(sharedTri)==1 || j == numel(tetloop))
            assert(numel(sharedTri)==4 || j ~= numel(tetloop))
%             plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'r');
%             if(exist('hndl'))
%                 delete(hndl);
%             end
%             hndl = plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'b','LineWidth',2);
%             pause(.5);
        end
    end
    for i = 1:numel(tettripleloopsfacetype)
        tetloop = tettripleloopsfacetype{i};
        tetloopExt = [tetloop tetloop(1)];
%         figure; hold on; axis equal; rotate3d on;
        for j = 1:numel(tetloop)
            pairtets = tetloopExt(j:j+1);
            sharedTri = intersect(data.tetsToTriangles(pairtets(1),:),data.tetsToTriangles(pairtets(2),:));
            assert(numel(sharedTri)==1 || j == numel(tetloop))
            assert(numel(sharedTri)==4 || j ~= numel(tetloop))
%             plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'r');
%             if(exist('hndl'))
%                 delete(hndl);
%             end
%             hndl = plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'b','LineWidth',2);
%             pause(.5);
        end
    end
    for i = 1:numel(tetdoubleloopscornertype)
        tetloop = tetdoubleloopscornertype{i};
        tetloopExt = [tetloop tetloop(1)];
%         figure; hold on; axis equal; rotate3d on;
        for j = 1:numel(tetloop)
            pairtets = tetloopExt(j:j+1);
            sharedTri = intersect(data.tetsToTriangles(pairtets(1),:),data.tetsToTriangles(pairtets(2),:));
            assert(numel(sharedTri)==1 || j == numel(tetloop))
            assert(numel(sharedTri)==4 || j ~= numel(tetloop))
%             plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'r');
%             if(exist('hndl'))
%                 delete(hndl);
%             end
%             hndl = plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'b','LineWidth',2);
%             pause(.5);
        end
    end
end

%% todo: do something ??? about the boundary

%% compute constrained equations based on loops
cornerHolonomyLoops = {};
faceHolonomyLoops = {};
edgeHolonomyLoops = {};

for i = 1:numel(tetdoubleloopscornertype)
    cornerHolonomyLoops{i} = transitionsPerTetLoop(tetdoubleloopscornertype{i},data,triangleToTransition);
end
for i = 1:numel(tettripleloopsfacetype)
    faceHolonomyLoops{i} = transitionsPerTetLoop(tettripleloopsfacetype{i},data,triangleToTransition);
end
for i = 1:numel(tettripleloopsedgetype)
    edgeHolonomyLoops{i} = transitionsPerTetLoop(tettripleloopsedgetype{i},data,triangleToTransition);
end
%% todo: compute constrained equations based on edgevalence and add them to faceHolonomyLoops

%% sort and prune holonomy constraints
[cornerHolonomyLoops,faceHolonomyLoops,edgeHolonomyLoops] = reduceHolonomyConstraints(cornerHolonomyLoops,faceHolonomyLoops,edgeHolonomyLoops);


%% verify metavertices and corners by visualizing
if(VisualizeLoops)
    for vertind = 1:numel(MetaVertices)
        mvert = MetaVertices{vertind}; vind = mvert.vind;
        adjTets = find(sum(data.tetrahedra==mvert.vind,2));
        adjVerts = ARemoveB(unique(data.tetrahedra(adjTets,:)),vind);
        figure; hold on; axis equal; rotate3d on;
        scatter3(data.vertices(adjVerts,1),data.vertices(adjVerts,2),data.vertices(adjVerts,3),'r.');
        plot3(data.vertices(vind,1),data.vertices(vind,2),data.vertices(vind,3),'k');
        
        for i = 1:numel(mvert.loops)
            tetloop = mvert.loops{i};
            tetloopExt = [tetloop tetloop(1)];
            hold on; axis equal; rotate3d on;
            for j = 1:numel(tetloop)
                pairtets = tetloopExt(j:j+1);
                plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'k','LineWidth',2);
                if(exist('hndl'))
                    delete(hndl);
                end
                hndl = plot3(data.tetBarycenters(pairtets,1),data.tetBarycenters(pairtets,2),data.tetBarycenters(pairtets,3),'b','LineWidth',4);
                pause(.1);
            end
        end
    end
end


if(Visualize)
    for h = 1:numel(H1Generators)
        f = VisualizeEdges(SEdges, data, '-', 0, [1 0 0]); title(['generator ' num2str(h)]);
        root = data.tetBarycenters(H1Generators{1}(1),:);
        scatter3(root(1),root(2),root(3),300,'r.');

        plot3(data.tetBarycenters(H1Generators{h},1),data.tetBarycenters(H1Generators{h},2),data.tetBarycenters(H1Generators{h},3));
        for i = 1:numel(MetaSurfaceClosed{h})
            if(StartingTri{h} == MetaSurfaceClosed{h}(i)); continue; end;
            vinds = data.triangles(MetaSurfaceClosed{h}(i),:);
            polyPtch = data.vertices(vinds, :);
            
            if(numel(triangleToTransition{MetaSurfaceClosed{h}(i)})==0)
                color = 'black'
            else
                color = 'green'
            end
            
            ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),color); alpha(ptc, .6);
            %pause(.01)
        end
        
        % Starting Triangle: where the loop was closed.
        vinds = data.triangles(StartingTri{h},:);
        polyPtch = data.vertices(vinds, :);
        ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'red'); alpha(ptc, .8);
        
        surfEdges = data.trianglesToEdges(MetaSurfaceClosed{h},:);
        edgeCount = accumarray(surfEdges(:),1);
        surfEdges = find(edgeCount ~= 0 & edgeCount ~= 2);
        VisualizeEdges(surfEdges, data, '-', f, [0 0 1]);
    end
    
    f = VisualizeEdges(SurfaceEdge, data, '-', 0, [0 0 1]); title('Verify Surface');
    VisualizeEdges(InteriorEdge, data, '-', f, [0 1 0]);
    VisualizeEdges(SEdges, data, '-', f, [1 0 0]);
    for i = 1:numel(MetaVertices)
        hold on;
        scatter3(data.vertices(MetaVertices{i}.vind,1),data.vertices(MetaVertices{i}.vind,2),data.vertices(MetaVertices{i}.vind,3),'r');
    end
    
    f = VisualizeEdges(SEdges, data, '-', 0, [1 0 0]); axis equal;
    title('Verify metavertices and edges')
    for i = 1:numel(MetaVertices)
        hold on;
        scatter3(data.vertices(MetaVertices{i}.vind,1),data.vertices(MetaVertices{i}.vind,2),data.vertices(MetaVertices{i}.vind,3),'r');
    end

    figure; hold on; axis equal; title('Show Corners');
    for i = 1:numel(MetaVertices)
        mvert = MetaVertices{i};
        scatter3(data.vertices(MetaVertices{i}.vind,1),data.vertices(MetaVertices{i}.vind,2),data.vertices(MetaVertices{i}.vind,3),'r');

        for j = 1:size(mvert.adjE2V,1)
            vs = data.vertices(mvert.adjE2V(j,:),:);
            plot3(vs(:,1),vs(:,2),vs(:,3),'r');
        end

        for c = 1:numel(mvert.corners)
            corner = mvert.corners{c};
            adddir = sum(data.vertices(corner.threeEdges2V(:,2),:)-data.vertices(corner.threeEdges2V(:,1),:))/3;
            for j=1:3
                vs = data.vertices(corner.threeEdges2V(j,:),:);
                vs = vs + [adddir;adddir];

                if(corner.threeEdgesDeg(j)==3)
                    plot3(vs(:,1),vs(:,2),vs(:,3),'b','LineWidth',2);
                else
                    plot3(vs(:,1),vs(:,2),vs(:,3),'g','LineWidth',2);
                end

                %corner.threeEdgesIsBoundary(j)

                %corners{pos}.threeEdges2V = threeEdges2V;
                %corners{pos}.threeEdgesDeg = threeEdgesDeg;
                %corners{pos}.threeEdgesIsBoundary = threeEdgesIsBoundary;
            end
        end
    end
end

%%






%end