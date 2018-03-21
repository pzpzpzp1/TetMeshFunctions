% given tet mesh data, and tetloop, a dual edge loop made of tets, and
% generators also expressed as tet loops, DualFacesToExclude as in dual edges
% in 1-ring of these edges will not be. Output is the tetloop expressed as
% ordered generators.
function gens = DualLoopToGenerators(data, dEloop, generators, DualFacesToExclude)
    numDualEdges = sum(~data.isBoundaryTriangle);
    G = sparse(numDualEdges,numel(generators));
    for i = 1:numel(generators)
        G(:,i) = sparse(generators{i},ones(numel(generators{i}),1),ones(numel(generators{i}),1),size(data.nonBoundaryTrianglesToTets,1),1);
    end
    nbE2nbTri = data.edgesToTrianglesIndicator(~data.isBoundaryEdge,~data.isBoundaryTriangle);
    G = [G nbE2nbTri'];
    
    % looking for x s.t. G*x=dEloop.
    



end