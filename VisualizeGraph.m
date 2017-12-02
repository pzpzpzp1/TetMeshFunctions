% very elementary viewing. Doesn't do much if any safegaurding against bad
% input. so use carefully.
function VisualizeGraph(edges, vertices, color, pauseTim, edgeSubset)

    if(nargin == 4)
        edgeSubset = 1:size(edges,1);
    end
    
    hold on; axis equal;
    startpoint = vertices(edges(edgeSubset(1),1)',:);
    scatter3(startpoint(:,1), startpoint(:,2), startpoint(:,3), 50);
    for i = 1:numel(edgeSubset)
        v = vertices(edges(edgeSubset(i),:)',:);
        plot3(v(:,1), v(:,2), v(:,3),color);
        pause(pauseTim);
    end
end