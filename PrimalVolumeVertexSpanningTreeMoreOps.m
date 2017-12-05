% startVert allows specifying of which connected component the tree is on.
% but in current use, there's only one component.
function T = PrimalVolumeVertexSpanningTreeMoreOps(edges, exclude, include, startVert)
    assert(numel(intersect(include, exclude))==0); % include needs to have no intersect with exclude. or it doesn't make sense.
    assert(numel(include) * numel(startVert) == 0); % specifying include means a startVert is easily obtained already. don't allow overspecification. that leads to bugs and confusing usage.
    
    edgesWithExclusion = edges; edgesWithExclusion(exclude,:)=[];
    edgesWithExclusionInds = 1:size(edges,1); edgesWithExclusionInds(exclude)=[];
    
    isIncluded = sparse(include,ones(numel(include),1),ones(numel(include),1)); isIncluded(size(edges,1)+1)=0;
    isIncluded (exclude)=[]; includeWithExclusion = find(isIncluded);
    
    verts = sort(unique(edgesWithExclusion));
    visitedV = [unique(edges(include,:))' startVert];
    completedV = [];
    Adjacency = sparse(edgesWithExclusion(:,1), edgesWithExclusion(:,2), 1:size(edgesWithExclusion,1));
    Adjacency(numel(verts),numel(verts))=0; 
    T = includeWithExclusion'; 
    
    while numel(completedV) ~= numel(verts)
        nleft = numel(verts) - numel(completedV);
        if(mod(nleft,1000)==0)
            fprintf('%d verts left ...\n',nleft);
        end
        
        if numel(visitedV) == 0
            % this part isn't actually reachable in the way this method is currently used. 
            % no point engineering a general purpose method when I only have a specific use. 
            % exclude won't disconnect vertices. it effectively removes some but the remainder is still connected.
            assert(0==1); 
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