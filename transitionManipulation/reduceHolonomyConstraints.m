function [cornerHolonomyLoops,faceHolonomyLoops,edgeHolonomyLoops] = reduceHolonomyConstraints(cornerHolonomyLoops,faceHolonomyLoops,edgeHolonomyLoops)
    numels = cellfun(@numel,cornerHolonomyLoops);
    [~, reorder] = sort(numels);
    cornerHolonomyLoops = cornerHolonomyLoops(reorder);
    cornerHolonomyLoops = cornerHolonomyLoops(numels(reorder)~=0);

    numels = cellfun(@numel,faceHolonomyLoops);
    [~, reorder] = sort(numels);
    faceHolonomyLoops = faceHolonomyLoops(reorder);
    faceHolonomyLoops = faceHolonomyLoops(numels(reorder)~=0);

    numels = cellfun(@numel,edgeHolonomyLoops);
    [~, reorder] = sort(numels);
    edgeHolonomyLoops = edgeHolonomyLoops(reorder);
    edgeHolonomyLoops = edgeHolonomyLoops(numels(reorder)~=0);
    
    loops = cornerHolonomyLoops;
    numels = cellfun(@numel,loops);
    fullreducedLoops = {};
    for i = unique(numels);
        subloops = loops(numels==i);
        reducedLoops = {}; rloopind = 1;
        for j = 1:numel(subloops)
            subloop = subloops{j};
            % check if subloop is in reducedLoops
            anymatch = false;
            for k = 1:numel(reducedLoops)
                rloop = reducedLoops{k};
                % check if rloop equals subloop up to cyclic or invert.
                if areHolonomyConstraintsEquivalent(rloop, subloop)
                    anymatch = true;
                    break;
                end
            end
            if(~anymatch)
                reducedLoops{rloopind} = subloop; rloopind = rloopind + 1;
            end
        end
        fullreducedLoops = [fullreducedLoops reducedLoops];
    end
    cornerHolonomyLoops = fullreducedLoops;
    
    loops = faceHolonomyLoops;
    numels = cellfun(@numel,loops);
    fullreducedLoops = {};
    for i = unique(numels);
        subloops = loops(numels==i);
        reducedLoops = {}; rloopind = 1;
        for j = 1:numel(subloops)
            subloop = subloops{j};
            % check if subloop is in reducedLoops
            anymatch = false;
            for k = 1:numel(reducedLoops)
                rloop = reducedLoops{k};
                % check if rloop equals subloop up to cyclic or invert.
                if areHolonomyConstraintsEquivalent(rloop, subloop)
                    anymatch = true;
                    break;
                end
            end
            if(~anymatch)
                reducedLoops{rloopind} = subloop; rloopind = rloopind + 1;
            end
        end
        fullreducedLoops = [fullreducedLoops reducedLoops];
    end
    faceHolonomyLoops = fullreducedLoops;
    
    loops = edgeHolonomyLoops;
    numels = cellfun(@numel,loops);
    fullreducedLoops = {};
    for i = unique(numels);
        subloops = loops(numels==i);
        reducedLoops = {}; rloopind = 1;
        for j = 1:numel(subloops)
            subloop = subloops{j};
            % check if subloop is in reducedLoops
            anymatch = false;
            for k = 1:numel(reducedLoops)
                rloop = reducedLoops{k};
                % check if rloop equals subloop up to cyclic or invert.
                if areHolonomyConstraintsEquivalent(rloop, subloop)
                    anymatch = true;
                    break;
                end
            end
            if(~anymatch)
                reducedLoops{rloopind} = subloop; rloopind = rloopind + 1;
            end
        end
        fullreducedLoops = [fullreducedLoops reducedLoops];
    end
    edgeHolonomyLoops = fullreducedLoops;
end


