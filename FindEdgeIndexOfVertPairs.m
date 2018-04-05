function [found, index] = FindEdgeIndexOfVertPairs(vertPairs, edges)

    flipinds = find(edges(:,1) > edges(:,2));
    edges(flipinds,:) = [edges(flipinds,2) edges(flipinds,1)];
    
    flipinds = find(vertPairs(:,1) > vertPairs(:,2));
    vertPairs(flipinds,:) = [vertPairs(flipinds,2) vertPairs(flipinds,1)];
    
    [found,index] = ismember(vertPairs, edges, 'rows');
end