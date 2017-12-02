%% takes matrix M and duplicates each row N times.
function Mdup = repRows(M, N)

    Mdup = repmat(M,1,N);
    Mdup = reshape(Mdup', size(M,2), size(M,1)*N)';

end