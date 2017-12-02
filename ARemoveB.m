function C = ARemoveB(A,B)
    assert(size(A,1)==1)
    assert(size(B,1)==1)

    removeIndicator = (sum(B == A',2) > 0)';
    A(removeIndicator) = [];
    C = A;
end