function [X2, T2, BrokenSingularEdges] = SubdivideTMeshAtVertex(X, T, vind, data)
    if(strcmp(class(data),'double') && data == 0)
        data = paul_getTetData(T,X,1);
    end
    
    adjEdges = find(sum(data.edges==vind,2)~=0);
    adjFaces = find(sum(data.triangles==vind,2)~=0);
    adjTets = find(sum(data.tetrahedra==vind,2)~=0);
    cageEdges = unique(ARemoveB(data.tetsToEdges(adjTets,:),adjEdges));
    cageVerts = unique(ARemoveB(data.tetrahedra(adjTets,:),vind));
    nV = size(X,1);
    ncMV = numel(cageVerts);
    ncMVMV = numel(cageEdges);
    
    X2 = X;
    T2 = data.tetrahedra;
    T2(adjTets,:)=[];
    
    % indexed by cageVerts
    cageMidwayVerts = (data.vertices(vind,:)+data.vertices(cageVerts,:))/2;
    [~,ind] = ismember(data.edges(cageEdges,:),cageVerts);
    % indexed by cageEdges
    cageMidwayVertsMidwayVerts = (cageMidwayVerts(ind(:,1),:) + cageMidwayVerts(ind(:,2),:))/2;
    % add to new vertex list. V + cageVerts + cageEdges
    X2 = [X2;cageMidwayVerts;cageMidwayVertsMidwayVerts];
    
    % put vind at the beginning, while maintaining orientation.
    adjT = T(adjTets,:);
    reorientInds = find(T(adjTets,4)==vind);
    cycle1Inds = find(T(adjTets,2)==vind);
    cycle2Inds = find(T(adjTets,3)==vind);
    adjT(reorientInds,:) = fliplr(adjT(reorientInds,:));
    adjT(cycle1Inds,:) = adjT(cycle1Inds,[2 3 1 4]);
    adjT(cycle2Inds,:) = adjT(cycle2Inds,[3 1 2 4]);
    assert(all(adjT(:,1)==vind));
    
    [~, iv2] = ismember(adjT(:,2), cageVerts);
    [~, iv3] = ismember(adjT(:,3), cageVerts);
    [~, iv4] = ismember(adjT(:,4), cageVerts);
    
    [~, ie23] = FindEdgeIndexOfVertPairs(adjT(:,[2 3]), data.edges(cageEdges,:));
    [~, ie34] = FindEdgeIndexOfVertPairs(adjT(:,[3 4]), data.edges(cageEdges,:));
    [~, ie42] = FindEdgeIndexOfVertPairs(adjT(:,[4 2]), data.edges(cageEdges,:));
    
    iv2=iv2+nV;
    iv3=iv3+nV;
    iv4=iv4+nV;
    ie23=ie23+nV+ncMV;
    ie34=ie34+nV+ncMV;
    ie42=ie42+nV+ncMV;
    
    % 4 inner
    T2Add = [...
    [adjT(:,1) ie42 iv2 ie23];...
    [adjT(:,1) ie23 iv3 ie34];...
    [adjT(:,1) ie34 iv4 ie42];...
    [adjT(:,1) ie23 ie34 ie42];...
    % 3 side peices
    [ie42 iv2 ie23 adjT(:,2)];...
    [ie23 iv3 ie34 adjT(:,3)];...
    [ie34 iv4 ie42 adjT(:,4)];...
    % inner 4 peices
    [ie42 adjT(:,3) adjT(:,4) adjT(:,2)];...
    [ie42 adjT(:,3) adjT(:,2) ie23];...
    [ie42 adjT(:,3) ie23 ie34];...
    [ie42 adjT(:,3) ie34 adjT(:,4)]];

    BrokenSingularEdges = [[adjT(:,1) iv2 adjT(:,2)];...
    [adjT(:,1) iv3 adjT(:,3)];...
    [adjT(:,1) iv4 adjT(:,4)]];
    BrokenSingularEdges = unique(BrokenSingularEdges,'rows');

    T2=[T2; T2Add];
end