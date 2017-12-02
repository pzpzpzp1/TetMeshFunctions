%% takes a list of quats q1 and q2. returns list of quaternions t12i and
% list of rotation matrics such that Ri*t12i(q1(F0))==q2(F0)

function [t12i, Ri] = computeTransitions(q1, q2, quats, R, rev)
    if(nargin == 0)
        N = randi(100,1);
        q1 = randQuat(N);
        q2 = randQuat(N);
    end
    if(nargin == 2 || nargin == 0)
        quats = octahedralGroup('quaternion');
        R = octahedralGroup('rotationCells');
        rev = octahedralGroup('inverseIndex');
    end
    assert(size(q1,2)==4 && size(q2,2)==4 && size(q1,1)==size(q2,1));
    N = size(q1,1);
    
    t12iter_stacked = zeros(N*24,4);
    for iter = 1:24
        t12iter_stacked((iter-1)*N+1 : (iter-1)*N+N  ,:) = quatmultiply(quatinv(q1), quatmultiply(repmat(quats(rev(iter),:), N,1), q2));
    end
    axangs = quat2axang(t12iter_stacked);
    angs = reshape(abs(axangs(:,4)), N, 24);
    [minangles,index] = min(angs');
    fullindex = [0:24:24*(N-1)];
    
    % final results
    t12i = t12iter_stacked(((index-1)*N+1)+[1:N]-1,:);
    Ri = R(index);
    
    %% checks that the method worked t12i(q1(F0))==R*q2(F0)
    if nargin==0
        for checki = 1:N
            right = R{index(checki)} * quatrotate(repRows(q2(checki,:),3),eye(3));
            left = quatrotate(t12i(checki,:), quatrotate(q1(checki,:), eye(3)));
            assert(norm(left-right)<.0001);
        end
    end
    
    assert(all(size(t12i)==size(q1)));
end

%% Sanity checks
function sanityChecks()
    
    quats = octahedralGroup('quaternion');
    R = octahedralGroup('rotationCells');
    rev = octahedralGroup('inverseIndex');

    F0 = eye(3);
    q1 = randQuat;
    q2 = randQuat;
    
    % find t12 s.t. t12(q1(F0))=q2(F0)
    t12 = quatmultiply(quatinv(q1), q2);
    assert(norm(quatrotate(t12, quatrotate(q1,F0))-quatrotate(q2,F0))<.0001);

    for iter=1:24
        invInd = rev(iter);
        
        % for iter, compute desired frame
        F2iter = R{iter} * quatrotate(q2, F0);

        % find q2iter s.t. q2iter(F0) = F2iter = R{iter} * q2(F0)        
        q2iter = quatmultiply(quats(invInd,:), q2);
        assert(norm(quatrotate(q2, quatrotate(quats(invInd,:), F0))-F2iter)<.00001);
        assert(norm(quatrotate(quatmultiply(quats(invInd,:), q2), F0)-F2iter)<.00001);
        assert(norm(quatrotate(q2iter,F0)-F2iter)<.0001);

        % find t12iter s.t. t12iter(q1(F0)) = q2iter(F0)
        t12iter = quatmultiply(quatinv(q1), quatmultiply(quats(invInd,:), q2));
        t12iter_associative = quatmultiply(quatmultiply(quatinv(q1), quats(invInd,:)), q2);
        assert(norm(t12iter - quatmultiply(quatinv(q1), q2iter))<.0001);
        assert(norm(quatrotate(t12iter, quatrotate(q1,F0))-quatrotate(q2iter,F0))<.0001);
        assert(norm(t12iter-t12iter_associative) <.00001);
    end
end