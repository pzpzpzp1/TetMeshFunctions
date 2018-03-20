function C = ARemoveB(A,B)
    A = A(:)';
    B = B(:)';

    removeIndicator = (sum(B == A',2) > 0)';
    A(removeIndicator) = [];
    C = A;
end