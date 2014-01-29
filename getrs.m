nx = 40;
ny = 40;
[xx yy] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));
[xx1 yy1] = meshgrid(linspace(0,1,nx+1),linspace(0,1,ny+1));
[xx2 yy2] = meshgrid(linspace(0,1,nx+2),linspace(0,1,ny+2));
u = frdr('u.dat',nx,ny);
uraw = frdr('uraw.dat',nx,ny+1);
us = frdr('us.dat',nx,ny+1);
v = frdr('v.dat',nx,ny);
vraw = frdr('vraw.dat',nx,ny);
vs = frdr('vs.dat',nx+1,ny);
omega = frdr('omega.dat',nx,ny);
phi = frdr('phi.dat',nx+1,ny+1);
strm = frdr('stream.dat',nx+2,ny+2);
getvo
figure(3)
quiver(xx,yy,u',v')
figure(4)
plot(xx(end/2,:),u(end/2,:))
xlabel('Y');
ylabel('U Velocity');
title('Central Channel Profile');
figure(5)
plot(xx(end/2,:),v(end/2,:))
xlabel('Y');
ylabel('V-Velocity');
title('Central Channel Profile');


