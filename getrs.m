nx = 40
ny = 40
[xx yy] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));
u = frdr('u.dat',nx,ny);
uraw = frdr('uraw.dat',nx,ny+1);
us = frdr('us.dat',nx,ny+1);
v = frdr('v.dat',nx,ny);
vraw = frdr('vraw.dat',nx,ny);
vs = frdr('vs.dat',nx+1,ny);
omega = frdr('omega.dat',nx,ny);
phi = frdr('phi.dat',nx+1,ny+1);
strm = frdr('stream.dat',nx+2,ny+2);
quiver(xx,yy,u',v')
