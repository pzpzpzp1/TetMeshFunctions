% load ray sokolov frames file
% n is number of vertices of the mesh.
function frames=LoadFramesFile(filename,n)
    if nargin==0
        filename = '..\..\jsolomon\octahedral_frames\comparison_data\FF_ray_sokolov\FFinit\elk_18k_frame.txt';
        n = 4748;
    end
    fid = fopen(filename,'r');
    R = fscanf(fid,'%f',9*n);
    frames = reshape(R,3,3,[]);
end