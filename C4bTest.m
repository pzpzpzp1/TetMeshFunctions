O = octahedralGroup('rotationCells');


count = 0;
x = [1;0;0];
for i = 1:24;
    for j = 1:24;
        for k = 1:24;
            R1 = O{i}*x;
            R2 = O{j}*x;
            R3 = O{k}*x;
            R4 = O{i}*O{j}*O{k};
            %if(abs(R1'*R2)<.000001 && abs(R3'*R1)<.000001 && abs(R2'*R3)<.000001 && abs(det([R1 R2 R3])-1)<.000001 && norm(R4-eye(3))<.0000001)
            if(abs(R1'*R2)<.000001 && abs(R3'*R1)<.000001 && abs(R2'*R3)<.000001 && abs(det([R1 R2 R3])-1)<.000001)
                [i j k]
                count = count + 1;
            end
        end
    end
end
count

count = zeros(1,24);
x = [1;0;0];
y = [0;1;0];
z = [0;0;1];
for j = 1:24;
    for k = 1:24;
        R2 = O{j}*x;
        R3 = O{k}*x;
        R4 = eye(3)*O{j}*O{k};
        %if(abs(x'*R2)<.000001 && abs(R3'*x)<.000001 && abs(R2'*R3)<.000001 && abs(det([x R2 R3])-1)<.000001 && norm(O{j}*y-y)<.00001 && norm(O{k}*z-z)<.00001)
        if(abs(x'*R2)<.000001 && abs(R3'*x)<.000001 && abs(R2'*R3)<.000001 && abs(det([x R2 R3])-1)<.000001)
            [j k]
            count(j) = count(j) + 1;
        end
    end
end
count


