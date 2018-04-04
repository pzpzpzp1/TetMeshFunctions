function accumTransitions = transitionsPerEdge(edgeind,data,triangleToTransition)
    assert(~data.isBoundaryEdge(edgeind));

    triCycle = data.edgeTriCycles{edgeind}(1:end-1);
    accumTransitions = [];
    for accumTris = triCycle
        % get transitions of a triangle.
        transition = triangleToTransition{accumTris};

        side1 = min(find(data.edgeCycles{edgeind}==data.trianglesToTets{accumTris}(1)));
        forward = data.edgeCycles{edgeind}(side1+1)==data.trianglesToTets{accumTris}(2);

        if(forward)
            accumTransitions = [accumTransitions transition];
        else
            accumTransitions = [accumTransitions fliplr(-1*transition)];
        end
    end
end