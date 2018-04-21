function [X2, T2] = SubdivideTMeshAtTets(X, T, Tinds)
    
    X2 = X;
    T2 = T;
    
    newX = permute(sum(reshape(X(T(Tinds,:)',:)',[3,4,numel(Tinds)]),2)/4,[1 3 2])';
    X2 = [X2;newX];
    Told = T2(Tinds,:);
    T2(Tinds,:)=[];
    nV = size(X,1); nnX = size(newX,1); nX2 = size(X2,1);
    
    nVInds = [nV+1:nX2]';
    T2 = [T2; [Told(:,1:3) nVInds]];
    T2 = [T2; [Told(:,[1 4 2]) nVInds]];
    T2 = [T2; [Told(:,[2 4 3]) nVInds]];
    T2 = [T2; [Told(:,[3 4 1]) nVInds]];
    
    assert(size(T2,1) == size(T,1)+3*numel(Tinds));
end