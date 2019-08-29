
%% CAN'T LOAD MIXED ELEMENT MESHES!!!
function writeOVM(filename,X,tets,properties)
    error('not written!!');
    %fid = fopen(filename,'w');
    
    
end

function [verts, edges, faces2verts, cells2verts, faces, cells, properties] = LoadOVM(filename)
    ovmfid = fopen(filename,'r');
    fscanf(ovmfid,'OVM %s\n');
    fscanf(ovmfid,'Vertices\n');
    nv = fscanf(ovmfid,'%d',1);
    verts = fscanf(ovmfid,'%f',3*nv);
    verts = reshape(verts,3,[])';
    
    fscanf(ovmfid,'\nEdges\n');
    nE = fscanf(ovmfid,'%d',1);
    edges = fscanf(ovmfid,'%d',2*nE);
    edges = reshape(edges,[],nE)'+1;
    
    fscanf(ovmfid,'\nFaces\n');
    nF = fscanf(ovmfid,'%d',1);
    faces = fscanf(ovmfid,'%d');
    faces = reshape(faces,[],nF)';
    faces = faces(:,2:end)+1;
    
    fscanf(ovmfid,'\nPolyhedra\n');
    nC = fscanf(ovmfid,'%d',1);
    cells = fscanf(ovmfid,'%d');
    cells = reshape(cells,[],nC)';
    cells = cells(:,2:end)+1;
    
    % convert half edge indices to edge indices 
    % convert half face indices to face indices 
    faces = ceil(faces/2);
    cells = ceil(cells/2);
    
    % extract face to vertices and cell to vertices
    faces2verts = [];
    for i=1:size(faces,2)
        faces2verts = [faces2verts edges(faces(:,i),:)];
    end
    faces2verts = sort(faces2verts,2);
    faces2verts = faces2verts(:,1:2:end);
    
    cells2verts = [];
    for i=1:size(cells,2)
        cells2verts = [cells2verts faces2verts(cells(:,i),:)];
    end
    cells2verts = sort(cells2verts,2);
    cells2verts = cells2verts(:,1:size(faces,2):end);
    
    % load properties if there are any
    
    properties={};
    propertiesLeft = true;
    while propertiesLeft
        proptype = fscanf(ovmfid,'\n%sProp');
        typename  = fscanf(ovmfid,' %s+ "');
        fieldname = fscanf(ovmfid,'%s+"\n'); fieldname = fieldname(2:end-1);
        
        if numel(proptype)==0
            break;
        end
        
        if strcmp(typename,'double')
            nC = fscanf(ovmfid,'%f');
        else
            nC = fscanf(ovmfid,'%d');
        end
        
        newstruct = struct(fieldname,{nC});
        newstruct.fieldname = fieldname;
        newstruct.proptype = proptype;
        properties{end+1} = newstruct;
    end
    
    
    fclose(ovmfid);
end

%% example ovm
%{
OVM ASCII
Vertices
59
0.320145 0.676773 -0.210422
-0.151702 -0.0359468 -0.140349
-0.409798 -0.410891 0.0382738
-0.256426 -0.411409 0.319033
-0.165833 -0.41095 -0.259362
0.276023 -0.410771 -0.366416
-0.164997 -0.410861 -0.186358
0.330274 0.705439 -0.340405
0.420801 0.697739 -0.281433
0.406717 0.678116 -0.112494
0.341807 0.70536 -0.265922
-0.00354654 -0.142634 0.266957
-0.383886 0.105565 0.488158
-0.0487128 -0.138246 0.325967
0.242065 -0.374765 -0.534241
-0.188865 0.70526 0.32093
0.120587 0.0211298 0.0133852
-0.41021 0.168494 -0.552018
0.0181195 0.67456 0.124116
-0.373386 0.431248 0.488146
-0.410326 0.672896 0.265783
0.0166591 -0.411409 0.402849
-0.410022 0.693697 -0.296046
0.216736 -0.410859 -0.536694
-0.148944 0.682631 0.12442
-0.0853676 0.705194 -0.41258
-0.409914 0.704529 -0.54432
-0.113918 -0.151414 0.2549
0.070973 -0.332476 0.488059
-0.410155 0.00304086 -0.557231
0.00854929 -0.411409 0.327735
-0.102242 -0.411409 0.299418
-0.195082 0.200819 0.488107
-0.136341 -0.251541 0.488458
-0.142665 -0.0618451 0.488009
-0.410408 0.641555 -0.528118
-0.367277 0.636676 -0.57734
-0.410092 0.683312 -0.0531555
-0.0834923 0.705292 0.0952132
-0.177708 0.55833 0.487878
-0.392532 0.307793 0.48805
0.0759695 -0.40516 0.186912
-0.376957 0.666584 0.488137
0.100712 0.705139 -0.57148
0.41593 -0.400705 -0.14379
0.265121 -0.401678 0.0676253
-0.409933 0.133715 0.485957
0.0760656 -0.353435 0.242789
-0.266406 0.705447 0.451527
0.065622 -0.147368 0.353397
0.0715572 -0.229658 0.406645
0.0661202 -0.41101 0.420107
0.312208 -0.0855332 -0.485127
0.407768 -0.127089 -0.335691
0.24832 -0.366713 0.0783011
-0.191166 -0.411409 0.305674
-0.410178 -0.10668 0.467231
-0.410075 -0.169177 0.218835
-0.15153 0.378612 0.0518961
Edges
0
Faces
0
Polyhedra
0
%}