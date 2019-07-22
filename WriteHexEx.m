
function WriteHexEx(filename, X, Tets, uvw)
    
    % uvw = uvwnew(Tets(:,1),:) uvwnew(Tets(:,2),:) uvwnew(Tets(:,3),:) uvwnew(Tets(:,4),:);

    outfid = fopen(filename,'w');
    fprintf(outfid,'%d\n',size(X,1));
    fprintf(outfid,'%f %f %f\n',X');
    fprintf(outfid,'%d\n',size(Tets,1));
    outmat = [Tets-1 uvw];
    fprintf(outfid,'%d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f\n',outmat');
    fclose(outfid);

end

%% readme hexex
%{

Input File Format:

n
v1_x v1_y v1_z
v2_x v2_y v2_z
...
vn_x vn_y vn_z
m
c1_i1 c1_i2 c1_i3 c1_i4 p_c1_i1_u p_c1_i1_v p_c1_i1_w p_c1_i2_u p_c1_i2_v p_c1_i2_w p_c1_i3_u p_c1_i3_v p_c1_i3_w p_c1_i4_u p_c1_i4_v p_c1_i4_w
c2_i1 c2_i2 c2_i3 c2_i4 p_c2_i1_u p_c2_i1_v p_c2_i1_w p_c2_i2_u p_c2_i2_v p_c2_i2_w p_c2_i3_u p_c2_i3_v p_c2_i3_w p_c2_i4_u p_c2_i4_v p_c2_i4_w
...
cm_i1 cm_i2 cm_i3 cm_i4 p_cm_i1_u p_cm_i1_v p_cm_i1_w p_cm_i2_u p_cm_i2_v p_cm_i2_w p_cm_i3_u p_cm_i3_v p_cm_i3_w p_cm_i4_u p_cm_i4_v p_cm_i4_w 

The input file starts with an integer n representing the number of vertices of the input mesh. It follow the n positions of the vertices, each as three floating point values in its own line.  After that there is again an integer m representing the number of tets which are to be read, followed by m tet definitions, each in its own line.  A tet definition consists of a line starting with 4 integers, which reference the vertices.  The vertices of a tet should be ordered such that if p1, p2, p3 and p4 are the corresponding positions, the determinant det(p2-p1,p3-p1,p4-p1) is positive.  The four vertex indices are followed by the four corresponding parameter locations, each specified as three floating point numbers. 
An example file should have been shipped with this readme.  

%}