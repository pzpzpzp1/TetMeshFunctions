% returns the octahedral group in various representations
function z = octahedralGroup(type)

    if(isa(type,'double') && all(size(type)==[3,3]))
        %input is matrix. find closest octahedral rotation, and type
        %(I,C,E,F) ident corner edge face
        O = octahedralGroup('rotationCells');
        faces = O(2:10);
        edges = O(11:16);
        corners = O(17:24);
        matching = zeros(24,1);
        for i = 1:24
            matching(i)=norm(O{i}-type);
        end
        [m,i] = min(matching);
        z = i;
        return;
    elseif(isa(type,'double') && numel(type)==1)
        if(type == 1) % ident
            z=0;
        elseif(ismember(type,[2:10])) % faces
            z = 1;
        elseif(ismember(type,[11:16])) % edges
            z = 2;
        elseif(ismember(type,[17:24])) % corners
            z = 3;
        end
        return;
    end

    if strcmp(type, 'symmetryPermutation')
        % front face of cube has corners labeled 1 2 3 4 ccw. opposite
        % corners are labeled the same.
        z = [1 2 3 4; 3 1 4 2 ;4 3 2 1   ; 2 4 1 3   ; 4 1 2 3   ; 3 4 1 2   ; 2 3 4 1   ; 3 4 2 1   ; 2 1 4 3   ; 4 3 1 2   ; 1 2 4 3   ; 4 2 3 1   ; 1 4 3 2   ; 2 1 3 4   ; 1 3 2 4   ; 3 2 1 4   ; 2 3 1 4   ; 4 2 1 3   ; 1 3 4 2   ; 4 1 3 2   ; 3 1 2 4   ; 3 2 4 1   ; 1 4 2 3   ; 2 4 3 1 ];
        %sum(sum(sortrows(perms([1:4]))==sortrows(z)))==numel(z); sanity check.
    elseif strcmp(type, 'quaternion')
        axang = [1 0 0 0;  % identity
            
            % face rotations
            1 0 0 pi/2; 1 0 0 pi; 1 0 0 3*pi/2; 
            0 1 0 pi/2; 0 1 0 pi; 0 1 0 3*pi/2; 
            0 0 1 pi/2; 0 0 1 pi; 0 0 1 3*pi/2; 
            
            % edge rotations
            1 1 0 pi;
            0 1 1 pi;
            1 0 1 pi;
            -1 1 0 pi;
            0 -1 1 pi;
            -1 0 1 pi;
            
            % corner rotations
            1 1 1 2*pi/3;
            1 1 -1 2*pi/3;
            1 -1 1 2*pi/3;
            1 -1 -1 2*pi/3;
            1 1 1 4*pi/3;
            1 1 -1 4*pi/3;
            1 -1 1 4*pi/3;
            1 -1 -1 4*pi/3;
            ];
        
        z = axang2quat(axang);
    elseif strcmp(type, 'axang')
        axang = [1 0 0 0;  % identity
            
            % face rotations
            1 0 0 pi/2; 1 0 0 pi; 1 0 0 3*pi/2; 
            0 1 0 pi/2; 0 1 0 pi; 0 1 0 3*pi/2; 
            0 0 1 pi/2; 0 0 1 pi; 0 0 1 3*pi/2; 
            
            % edge rotations
            1 1 0 pi;
            0 1 1 pi;
            1 0 1 pi;
            -1 1 0 pi;
            0 -1 1 pi;
            -1 0 1 pi;
            
            % corner rotations
            1 1 1 2*pi/3;
            1 1 -1 2*pi/3;
            1 -1 1 2*pi/3;
            1 -1 -1 2*pi/3;
            1 1 1 4*pi/3;
            1 1 -1 4*pi/3;
            1 -1 1 4*pi/3;
            1 -1 -1 4*pi/3;
            ];
        
        z = axang;
    elseif strcmp(type, 'rotationCells')
        quats = octahedralGroup('quaternion');
        quatsdup = repRows(quats,3);
        %quatsdup = repmat(quats,1,3);
        %quatsdup = reshape(quatsdup', 4, 24*3)';
        mats = quatrotate(quatsdup, repmat(eye(3),24,1));
        
        mats(find(mats>.5))=1;
        mats(find(mats<-.5))=-1;
        mats(find(mats<.5 & mats>-.5))=0;
        
        z = mat2cell(mats,ones(24,1)*3,3);
        z = cellfun(@transpose,z,'un',0);
        
        % sanity checks
        %{
        randv = rand(100,3); i = randi(24,100);
        for iter = 1:100
            left = (rotationCells{i(iter)}*randv(iter,:)')';
            right = quatrotate(quats(i(iter),:), randv(iter,:));
            assert(norm(left-right)<.0001);
        end
        %}
        
        % further sanity checks. also useful to see the relationships
        %{
        ypr = rand(1,3)*2*pi;
        randq = angle2quat(ypr(:,1), ypr(:,2), ypr(:,3));
        randv = rand(1,3);
        i = randi(24,1); 
        
        % NOTE quatrot(quatmult(q1,q2),v) IS NOT quatrot(q1, quatrot(q2, v))
        left = (rotationCells{i}*quatrotate(randq, randv)')';
        right = quatrotate(quatmultiply(randq, quats(i,:)), randv);
        right1 = quatrotate(quats(i,:), quatrotate(randq, randv));
        
        assert(norm(left-right) + norm(right-right1) < .000001);
        %}
        
    elseif strcmp(type, 'rotationMats')
        z = octahedralGroup('rotationCells');
        z = cell2mat(z);
        
        z(find(z>.5))=1;
        z(find(z<-.5))=-1;
        z(find(z<.5 & z>-.5))=0;
    elseif strcmp(type, 'inverseIndex')
        z = [1 4 3 2 7 6 5 10 9 8 11 12 13 14 15 16 21 22 23 24 17 18 19 20];
        
        % sanity check
        %{
        rcs = octahedralGroup('rotationCells');
        for i = 1:24
            assert(norm(rcs{i}*rcs{z(i)}-eye(3))<.00001);
        end
        %}
    elseif strcmp(type, 'name')
        z = {'I'
            ,'fx1/4'
            ,'fx1/2'
            ,'fx3/4'
            ,'fy1/4'
            ,'fy1/2'
            ,'fy3/4'
            ,'fz1/4'
            ,'fz1/2'
            ,'fz3/4'
            ,'exy'
            ,'eyz'
            ,'exz'
            ,'e-xy'
            ,'e-yz'
            ,'e-xz'
            ,'c++1/3'
            ,'c+-1/3'
            ,'c-+1/3'
            ,'c--1/3'
            ,'c++2/3'
            ,'c+-2/3'
            ,'c-+2/3'
            ,'c--2/3'};
    else
        error('type not supported');
    end
end