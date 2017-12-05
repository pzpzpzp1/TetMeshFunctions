function [X,T] = loadTetGenMesh(baseFile)

% Read vertices
nodefile = [baseFile '.node'];
fid = fopen(nodefile,'r');
metadata = fscanf(fid,'%g',4);
n = metadata(1);
X = fscanf(fid,'%g',4*n);
X = reshape(X,4,[])';
X(:,1) = [];
fclose(fid);

% Read tets
tetfile = [baseFile '.ele'];
fid = fopen(tetfile,'r');
metadata = fscanf(fid,'%d',3);
n = metadata(1);
T = fscanf(fid,'%d',5*n);
T = reshape(T,5,[])';
T(:,1) = [];
T = T + 1;
fclose(fid);