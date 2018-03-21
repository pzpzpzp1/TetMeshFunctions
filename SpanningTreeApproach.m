Visualize = true;

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

%% NOTICE: given a dual loop, how do you express it in terms of generators?
interiorSingularEdges = intersect(find(~data.isBoundaryEdge),SEdges);
triloop = data.edgesToTrianglesUnoriented{interiorSingularEdges(1)};
%gens = DualLoopToGenerators(data, triloop, H1DualEdgeGenerators, DualFacesToExclude);



%% verify metavertices and corners
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
            ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .6);
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







%end