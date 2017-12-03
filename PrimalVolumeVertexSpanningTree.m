function T = PrimalVolumeVertexSpanningTree(edges) %, exclude)
    %if(argin == 1)
    %    exclude = [];
    %end

    verts = sort(unique(edges(:)));
    visitedV = verts(1);
    completedV = [];
    Adjacency = sparse(edges(:,1), edges(:,2), 1:size(edges,1));
    Adjacency(numel(verts),numel(verts))=0; 
    T = []; 
    
    while numel(completedV) ~= numel(verts)
        nleft = numel(verts) - numel(completedV);
        if(mod(nleft,1000)==0)
            fprintf('%d verts left ...\n',nleft);
        end
        
        if numel(visitedV) == 0
            error('Graph is not connected!');
        end
        v = visitedV(1);
        visitedV = visitedV(2:end);
        completedV = [completedV v];
        
        % get neighbors of v
        neighboringV1 = find(Adjacency(v,:));
        neighboringV2 = find(Adjacency(:,v));
        candidateEdges1 = Adjacency(v,neighboringV1);
        candidateEdges2 = Adjacency(neighboringV2,v);
        neighboringV = [neighboringV1 neighboringV2'];
        candidateEdges = [candidateEdges1 candidateEdges2'];
        
        % get neighboringV that aren't in visitedV or completedV
        RemoveV = [visitedV completedV];
        removeIndicator = (sum(RemoveV == neighboringV',2) > 0)';
        neighboringV(removeIndicator ) = [];
        candidateEdges(removeIndicator ) = [];
        
        % add candidate edges to tree, and new Vs to visited.
        visitedV = [visitedV neighboringV];
        T = [T candidateEdges];
    end
    T = full(T);
end