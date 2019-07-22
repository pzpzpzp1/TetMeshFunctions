% load ray sokolov mesh file
function [X,T]=LoadTetFile(filename)
    if nargin==0
        filename = 'D:\Documents\hexdom_geogram_1.6.1\example.tet';
    end
    fid = fopen(filename,'r');
    n = fscanf(fid,'%d vertices\n',1);
    m = fscanf(fid,'%d tets\n',1);
    
    X = fscanf(fid,'%f',3*n);
    X = reshape(X,3,[])';
    
    T = fscanf(fid,'%f',5*m);
    T = reshape(T,5,[])';
    T=T(:,2:end)+1;
    fclose(fid);
end