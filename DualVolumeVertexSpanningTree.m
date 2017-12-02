function DT = DualVolumeVertexSpanningTree(data)
    m = data.nonBoundaryTrianglesToTets;

    DT = PrimalVolumeVertexSpanningTree(m);
end