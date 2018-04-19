
% subdivide boundary tets with more than 2 boundary triangles
% permute tets so boundary tets are at the end
function [X2,T3,data3,tetPermutation] = preprocess_data(X,T)

    data = paul_getTetData(T,X,1,1);
    tets2TriSprse = sparse(repmat(1:data.numTetrahedra,1,4)', data.tetsToTriangles(:), ones(size(data.tetsToTriangles(:))), data.numTetrahedra, data.numTriangles);
    badTets = find(sum(tets2TriSprse(:,find(data.isBoundaryTriangle)),2)>1);
    
    X2 = X; T2 = T;
    if(numel(badTets)>0)
    
        % compute new splitpoints.
        xs = X(T(badTets,:)',:);
        c1 = sum(reshape(xs(:,1),4,[]));
        c2 = sum(reshape(xs(:,2),4,[]));
        c3 = sum(reshape(xs(:,3),4,[]));
        X2 = [X; [c1' c2' c3']/4];

        % 123 214 134 324
        newVinds = [size(X,1)+1:size(X2,1)]';
        T2 = T;
        badT = T2(badTets,:);
        T2(badTets,:)=[];

        Bt2V = data.triangles(find(data.isBoundaryTriangle),:);
        BT2t = data.tetsToTriangles(badTets,:);
        BT2tv = reshape(data.triangles(BT2t',:)',12,[])';
        newT = [[BT2tv(:,1:3) newVinds]; [BT2tv(:,4:6)  newVinds]; [BT2tv(:,7:9)  newVinds]; [BT2tv(:,10:12) newVinds]];
        T2 = [T2; newT];
    end
    
    data2 = paul_getTetData(T2,X2,1,1);
    T3 = [T2(~data2.isBoundaryTet,:); T2(data2.isBoundaryTet,:)];
    tetPermutation = [find(~data2.isBoundaryTet); find(data2.isBoundaryTet)];
    
    if(nargout > 2)
        data3 = paul_getTetData(T3,X2,0,1);
    end
    
    %{
        hold on;
        scatter3(X(:,1),X(:,2),X(:,3))
        scatter3(c1/4,c2/4,c3/4,'r')
        for i = 1:size(data2.edges,1)
            vs = data2.vertices(data2.edges(i,:),:);
            plot3(vs(:,1),vs(:,2),vs(:,3),'k');
        end
    %}
end