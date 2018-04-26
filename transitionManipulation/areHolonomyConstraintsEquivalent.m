function z = areHolonomyConstraintsEquivalent(constr1, constr2)
    z = false;
    constr1 = cancelAntiPairs(constr1);
    constr2 = cancelAntiPairs(constr2);
    if(numel(constr2)~=numel(constr1))
        z = false; return;
    end
    
    for i = 0:numel(constr1)
        if(all(constr1==circshift(constr2,i)))
            z = true; 
        end
    end

    constr1 = fliplr(-constr1);
    for i = 0:numel(constr1)
        if(all(constr1==circshift(constr2,i)))
            z = true; 
        end
    end
    
end


