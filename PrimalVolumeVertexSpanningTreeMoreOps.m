function T = PrimalVolumeVertexSpanningTreeMoreOps(edges, exclude, include, startVert)
    assert(numel(include)==0); % include isn't supported yet. 
    edgesWithExclusion = edges; edgesWithExclusion(exclude,:)=[];
    edgesWithExclusionInds = 1:size(edges,1); edgesWithExclusionInds(exclude)=[];
    
    
    verts = sort(unique(edgesWithExclusion));
    visitedV = startVert;
    completedV = [];
    Adjacency = sparse(edgesWithExclusion(:,1), edgesWithExclusion(:,2), 1:size(edgesWithExclusion,1));
    Adjacency(numel(verts),numel(verts))=0; 
    T = []; 
    
    while numel(completedV) ~= numel(verts)
        nleft = numel(verts) - numel(completedV);
        if(mod(nleft,1000)==0)
            fprintf('%d verts left ...\n',nleft);
        end
        
        if numel(visitedV) == 0
            'Graph is not connected! Propagation finished'
            T = full(edgesWithExclusionInds(T));
            return;
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
    T = full(edgesWithExclusionInds(T));
end