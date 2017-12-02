function rq = randQuat(N)

    if(nargin == 0)
        N = 1;
    end

    ypr = rand(N,3)*2*pi;
    rq = angle2quat(ypr(:,1), ypr(:,2), ypr(:,3));
end