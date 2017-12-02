% used by paul_getTetData
function z = findOrderOfAInB(A,B)
    assert(numel(B) >= numel(A));
    
    if(numel(A)==2)
        if(numel(B)==2)
            z = A(1)==B(1) && A(2)==B(2);
            return;
        end
        
        ind1 = find(B==A(1), 1);
        ind2 = find(B==A(2), 1);
        z = (mod(ind2,  numel(B)) == mod(ind1 + 1, numel(B)));
    elseif numel(A)==1
        z = B(1)==A;
    else
        assert(false);
    end
end
