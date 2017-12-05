% edges need to index into data.edges.
function f = VisualizeDualSurface(data, edges, type, f)
    if(nargin == 3)
        f = figure(); hold on; axis equal;
    end
    
    figure(f);
    
    if type == '.'
        boundaryverts = data.vertices(find(data.isBoundaryVertex),:);
        scatter3(boundaryverts(:,1),boundaryverts(:,2),boundaryverts(:,3),1);
        
        EdgesRemaining = edges;
        vertinds = data.edges(EdgesRemaining,:);
        vlocs = (data.vertices(vertinds(:,1),:) + data.vertices(vertinds(:,2),:))/2;
        scatter3(vlocs(:,1),vlocs(:,2),vlocs(:,3),1);
        
    elseif type == 'p'
        boundaryverts = data.vertices(find(data.isBoundaryVertex),:);
        scatter3(boundaryverts(:,1),boundaryverts(:,2),boundaryverts(:,3),1);
        
        EdgesRemaining = edges;
        
        for e = 1:numel(EdgesRemaining)
            edgeInd = EdgesRemaining(e);
            tets = data.edgeCycles{edgeInd};
            polyPtch = data.tetBarycenters(tets',:);
            ptc = patch(polyPtch(:,1), polyPtch(:,2), polyPtch(:,3),'green'); alpha(ptc, .05);
            pause(.001);
        end
    else
        error('type not supported');
    end
    
    


end

