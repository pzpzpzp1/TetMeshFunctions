function transitions = cancelAntiPairs(transitions)

    inds2rmv = find(transitions == -circshift(transitions,-1));
    while(numel(inds2rmv)~=0 && numel(transitions)>1)
        cind = inds2rmv(1);
        nind = cind + 1;
        if(cind == numel(transitions)); nind = 1; end;
        transitions([cind nind])=[];

        inds2rmv = find(transitions == -circshift(transitions,-1));
    end 
end