N = 1000; rad = 1.9; da = .1;
x = rand(N,3)*rad*pi-rad*pi/2;
rminds = find(sqrt(sum(x.*x,2)) > pi);
x(rminds,:)=[];

V = x;
angle = sqrt(sum(x.*x,2));
axes = x./angle;

axang = [axes, angle];
q=axang2quat([axes,angle]);
n = size(q,1);
qxdir = repmat(axang2quat([1,0,0,da]),n,1);
qydir = repmat(axang2quat([0,1,0,da]),n,1);
qzdir = repmat(axang2quat([0,0,1,da]),n,1);

qx = quatmultiply(q,qxdir);
qy = quatmultiply(q,qydir);
qz = quatmultiply(q,qzdir);

axang = quat2axang(q);
xaxang = quat2axang(qx);
yaxang = quat2axang(qy);
zaxang = quat2axang(qz);

V1 = axang(:,1:3).*axang(:,4);
V2 = axang(:,1:3).*axang(:,4);
V3 = axang(:,1:3).*axang(:,4);
Vx = xaxang(:,1:3).*xaxang(:,4);
Vy = yaxang(:,1:3).*yaxang(:,4);
Vz = zaxang(:,1:3).*zaxang(:,4);

inds1 = find(sqrt(sum((V1-Vx).*(V1-Vx),2)) < pi);
inds2 = find(sqrt(sum((V2-Vy).*(V2-Vy),2)) < pi);
inds3 = find(sqrt(sum((V3-Vz).*(V3-Vz),2)) < pi);

interp = 0:.01:1;
xscatter = [];
for i = 1:numel(interp)
    xscatter = [xscatter;V(inds1,:)*interp(i)+Vx(inds1,:)*(1-interp(i))];
end
yscatter = [];
for i = 1:numel(interp)
    yscatter = [yscatter;V(inds2,:)*interp(i)+Vy(inds2,:)*(1-interp(i))];
end
zscatter = [];
for i = 1:numel(interp)
    zscatter = [zscatter;V(inds3,:)*interp(i)+Vz(inds3,:)*(1-interp(i))];
end

hold on
scatter3(V(:,1),V(:,2),V(:,3))
scatter3(xscatter(:,1),xscatter(:,2),xscatter(:,3),.4,'red')
scatter3(yscatter(:,1),yscatter(:,2),yscatter(:,3),.4,'green')
scatter3(zscatter(:,1),zscatter(:,2),zscatter(:,3),.4,'blue')

