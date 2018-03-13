% load ray sokolov mesh file
function [X,T]=LoadMeshFile(filename)
    if nargin==0
        filename = '..\..\jsolomon\octahedral_frames\comparison_data\FF_ray_sokolov\FFinit\elk_18k.mesh';
    end
    fid = fopen(filename,'r');
    metadata = fscanf(fid,'%s',5);
    n = fscanf(fid,'%d',1);
    X = fscanf(fid,'%f',4*n);
    X = reshape(X,4,[])';
    X = X(:,1:3);

    metadata = fscanf(fid,'%s',1);
    t = fscanf(fid,'%d',1);
    tri = fscanf(fid,'%f',4*t);
    metadata = fscanf(fid,'%s',1);
    E = fscanf(fid,'%d',1);
    edges = fscanf(fid,'%f',3*E);
    metadata = fscanf(fid,'%s',1);
    nT = fscanf(fid,'%d',1);
    T = fscanf(fid,'%d',5*nT);
    T = reshape(T,5,[])';
    T = T(:,1:4);
end