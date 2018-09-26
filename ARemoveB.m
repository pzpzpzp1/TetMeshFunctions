function A = ARemoveB(A,B)
    A = A(:)';
    B = B(:)';

    inds = ismember(A,B);
    A(inds)=[];
    
%     removeIndicator = (sum(B == A',2) > 0)';
%     A(removeIndicator) = [];
%     C = A;
%    assert(all(sort(C)==sort(A1)))
end