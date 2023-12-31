function f = VisualizeEdges(edgeInds, data, type, f, color)
    if nargin == 3 || f == 0
        f = figure; hold on; axis equal;
    end
    if nargin == 4
        color = [1 0 1];
    end
    figure(f);

    if type == '.'
        vertInds = data.edges(edgeInds,:);
        p1=data.vertices(vertInds(:,1),:);
        p2=data.vertices(vertInds(:,2),:);
        interp = 0:.04:1;
        final = [];
        for i = 1:numel(interp)
            final = [final; p1*interp(i) + p2*(1-interp(i))];
            
        end
        scatter3(final(:,1),final(:,2),final(:,3),.1,'green');
        
    elseif type == '-'
        for i = 1:numel(edgeInds)
            e = data.edges(edgeInds(i),:);
            vs = data.vertices(e',:);
            plot3(vs(:,1),vs(:,2),vs(:,3),'Color',color,'LineWidth',1.5);
        end
    else
        error('type not supported');
    end
    
end