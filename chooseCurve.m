%% create test singular primal curve
function [curveEdges, curveV, f] = chooseCurve(data, f, color, otherCurveEdges)
    if color == 0
        color = 'b'
    end

    if f==0
        f = figure; hold on; axis equal;
    else
        figure(f);
    end

    donechoosing  = false;
    while ~donechoosing 
        hold off; scatter3(data.vertices(:,1),data.vertices(:,2),data.vertices(:,3),.01,'b'); hold on;
        for i = 1:numel(otherCurveEdges)
            VisualizeEdges(otherCurveEdges{i}, data, '.', f)
        end
        
        nonboundaryedges = data.edges(find(~data.isBoundaryEdge),:);
        Adj = sparse(nonboundaryedges(:,1),nonboundaryedges(:,2),1:(size(nonboundaryedges,1)));
        startedge = randi(size(nonboundaryedges,1));
        v = data.vertices(nonboundaryedges(startedge,:)',:);
        endpoints = nonboundaryedges(startedge,:);
        plot3(v(:,1),v(:,2),v(:,3),color); desdir = v(2,:)-v(1,:); desdir=desdir/norm(desdir);
        notdone = any(~data.isBoundaryVertex(endpoints));
        curve = startedge;
        curveV = nonboundaryedges(startedge,:);
        while notdone
            for endpointInd = 1:2    
                if(~data.isBoundaryVertex(endpoints(endpointInd)))
                    e1 = endpoints(endpointInd);

                    candidateNv1 = find(Adj(e1,:));
                    candidateNv2 = find(Adj(:,e1));

                    candidateEdges1 = Adj(e1,candidateNv1);
                    candidateEdges2 = Adj(candidateNv2,e1);

                    candidateNv = [candidateNv1 candidateNv2'];
                    candidateEdges = [candidateEdges1 candidateEdges2'];

                    RemoveV = [curveV];
                    removeIndicator = (sum(RemoveV == candidateNv',2) > 0)';
                    candidateNv(removeIndicator ) = [];
                    candidateEdges(removeIndicator ) = [];

                    choice = randi(numel(candidateNv));
                    dists = data.vertices(candidateNv',:) - repmat(data.vertices(e1,:),numel(candidateNv),1);
                    dists = dists./sqrt(sum(dists.*dists,2));
                    factor = (endpointInd-1)*2-1;
                    [a b] = max(factor * (sum(dists.*repmat(desdir, numel(candidateNv),1),2)));
                    choice = b;

                    curveV = [curveV candidateNv(choice)];
                    curve = [curve candidateEdges(choice)];

                    endpoints(endpointInd)=candidateNv(choice);
                    v = data.vertices([e1 candidateNv(choice)]',:);
                    plot3(v(:,1),v(:,2),v(:,3),color);
                    %pause(1);
                end
            end

           notdone = any(~data.isBoundaryVertex(endpoints));
        end
        scatter3(data.vertices(endpoints',1),data.vertices(endpoints',2),data.vertices(endpoints',3),1);
        assert(all(data.isBoundaryVertex(endpoints)));
        donechoosing = input('Is this curve good enough? Answer 0 or 1. ');
    end
    inds = find(~data.isBoundaryEdge); curveEdges = inds(curve); % reindex curve to data.edges.
    fprintf("Curve accepted \n");

end