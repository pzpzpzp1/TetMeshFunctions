function accumTransitions = transitionsPerTetLoop(tetloop,data,triangleToTransition)
    assert(max(tetloop)<data.numTetrahedra);
    tetloopExt = [tetloop tetloop(1)];
    
    accumTransitions = [];
    for j = 1:numel(tetloop)-1
        pairtets = tetloopExt(j:j+1);
        sharedTri = intersect(data.tetsToTriangles(pairtets(1),:),data.tetsToTriangles(pairtets(2),:));
        assert(numel(sharedTri)==1 || j == numel(tetloop));
        
        transition = triangleToTransition{sharedTri};
        forward = data.trianglesToTets{sharedTri}(1) == pairtets(1) && data.trianglesToTets{sharedTri}(2) == pairtets(2);
        if(~forward)
            assert(data.trianglesToTets{sharedTri}(2) == pairtets(1) && data.trianglesToTets{sharedTri}(1) == pairtets(2));
        end
        
        if(forward)
            accumTransitions = [accumTransitions transition];
        else
            accumTransitions = [accumTransitions fliplr(-1*transition)];
        end
    end
    
    accumTransitions = cancelAntiPairs(accumTransitions);
end