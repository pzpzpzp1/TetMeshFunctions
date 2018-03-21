%% Test loop corner constraints
O = octahedralGroup('rotationCells');
rx = O{2};
ry = O{5};
rz = O{8};

rots = {rx, rx', ry, ry', rz, rz'};

resultstrings = [];
pos = 1;
i1 = []; i2 = []; i3 = [];
isCxy = []; isCyz = []; isCzx = [];
xyzType = [];
numVal5 = [];
validCorner = [];
for x=1:6
    for y=1:6
        for z = 1:6
            r1 = rots{x};
            r2 = rots{y};
            r3 = rots{z};
            
            pxy = r1*r2;
            pyz = r2*r3;
            pzx = r3*r1;
            
            pxyz = r1*r2*r3;
            
            isCxyt = octahedralGroup(octahedralGroup(pxy))==3;
            isCyzt = octahedralGroup(octahedralGroup(pyz))==3;
            isCzxt = octahedralGroup(octahedralGroup(pzx))==3;
            
            type = octahedralGroup(octahedralGroup(pxyz));
            
            numVal5t = sum(ismember([x,y,z],[2,4,6]));
            
            c1 = ismember(x,[1 2]) & ismember(y,[3 4]) & ismember(z,[5 6]);
            c2 = ismember(x,[3 4]) & ismember(y,[5 6]) & ismember(z,[1 2]);
            c3 = ismember(x,[5 6]) & ismember(y,[1 2]) & ismember(z,[3 4]);
            validcornert = c1 || c2 || c3;
            
            % save results
            i1(pos) = x;
            i2(pos) = y;
            i3(pos) = z;
            
            isCxy(pos) = isCxyt;
            isCyz(pos) = isCyzt;
            isCzx(pos) = isCzxt;
            
            xyzType(pos) = type;
            
            numVal5(pos) = numVal5t;
            
            validCorner(pos) = validcornert;
            
            pos = pos + 1;
        end
    end
end

T = table(i1',i2',i3',numVal5',validCorner',isCxy',isCyz',isCzx',xyzType');
T.Properties.VariableNames = {'i1','i2','i3','numVal5','validCorner','isCxy','isCyz','isCzx','xyzType'};

isEven = mod(numVal5,2);
isEdge = xyzType == 2;
isFace = xyzType == 1;
cornercond = (isEven & isEdge) | (~isEven & isFace);
pairtest = isCxy&isCyz&isCzx;
fulltest = pairtest&cornercond;

redT = table(validCorner', pairtest', cornercond',fulltest');
redT.Properties.VariableNames = {'validCorner','pairtest','tripletest','fulltest'};

failind = find(fulltest ~= validCorner)



