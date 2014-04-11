function getrs(nx,ny,nprefix)

    nghost = 3;
    dx = 1/(nx-1);
    dy = 1/(ny-1);
    close all
    load('./ghia/Re100-two.mat');
    [xx, yy] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));
    [xx1, yy1] = meshgrid(linspace(0,1,nx+1),linspace(0,1,ny+1));
    [xx2, yy2] = meshgrid(linspace(0,1,nx+2),linspace(0,1,ny+2));
    [xxg, yyg] = meshgrid(linspace(dx/2-(nghost*dx),1-dx/2+(nghost*dx),nx+2*nghost-1),...
        linspace(dy/2-(nghost*dy),1-dy/2+(nghost*dy),ny+2*nghost-1));
    u = frdr('u.dat',nx+2*nghost-1,ny+2*nghost-1);
    uraw = frdr('uraw.dat',nx+2*nghost-1,ny+2*nghost);
    us = frdr('us.dat',nx,ny+1);
    v = frdr('v.dat',nx+2*nghost-1,ny+2*nghost-1);
    vraw = frdr('vraw.dat',nx+2*nghost+1,ny+2*nghost);
    vs = frdr('vs.dat',nx+1,ny);
    omega = frdr('omega.dat',nx,ny);
    phi = frdr('phi.dat',nx+1,ny+1);
    strm = frdr('stream.dat',nx+2,ny+2);
    G = frdr('G.dat',nx+2*nghost-1, ny+2*nghost-1);
    alpha = frdr('alpha.dat',nx+2*nghost-1, ny+2*nghost-1);
    d = frdr('d.dat',nx+2*nghost-1, ny+2*nghost-1);
    [lsegx lsegy] = frdrsegs('lsegs.dat',nx+2*nghost-1, ny+2*nghost-1);
%     G1 = frdr('G1.dat',nx+2*nghost-1, ny+2*nghost-1);
%     Ge1 = frdr('G.1',nx+2*nghost-1, ny+2*nghost-1);
%     Ge2 = frdr('G.2',nx+2*nghost-1, ny+2*nghost-1);
%     Ge3 = frdr('G.3',nx+2*nghost-1, ny+2*nghost-1);
%     Ge0 = frdr('G.0',nx+2*nghost-1, ny+2*nghost-1);
%     
%     dGdt = frdr('dGdt.dat',nx+2*nghost-1, ny+2*nghost-1);
     Markers = frdrint('Markers.dat',nx+2*nghost-1, ny+2*nghost-1);
%     Markers0 = frdrint('Markers.0',nx+2*nghost-1, ny+2*nghost-1);
%     Markers1 = frdrint('Markers.1',nx+2*nghost-1, ny+2*nghost-1);
%     Markers2 = frdrint('Markers.2',nx+2*nghost-1, ny+2*nghost-1);
%     Markers30 = frdrint('Markers.30',nx+2*nghost-1, ny+2*nghost-1);
%     Markers40 = frdrint('Markers.40',nx+2*nghost-1, ny+2*nghost-1);
%     Markers50 = frdrint('Markers.50',nx+2*nghost-1, ny+2*nghost-1);
%     Markers100 = frdrint('Markers.100',nx+2*nghost-1, ny+2*nghost-1);
    
%     %% Fig 1
%     [xx, yy] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));
%     figure(1);
%     contour(xx,yy,strm(2:nx+1,2:ny+1)',strmlevels)
%     axis equal
%     axis([0 1 0 1])
%     xlabel('x');
%     ylabel('y');
%     title('Stream Function Contour Re=100');
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 4 4])
%     print( gcf, '-dpng', [nprefix 'Stream'])
%     
%     %% Fig 2
%     figure(2)
%     contour(xx,yy,-omega',omegalevels)
%     axis equal
%     axis([0 1 0 1])
%     xlabel('x');
%     ylabel('y');
%     title('Vorticity Contour Re=100');
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 4 4])
%     print( gcf, '-dpng', [nprefix 'Vort'])
%     
    %% Fig 3
%     figure(3)
%     quiver(xxg,yyg,u',v',2)
%     grid on
%     axis equal
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 4 4])
%     print( gcf, '-dpng', [nprefix 'Quiv'])
%  
%     
%     umax = max(u(ceil(end/2),:));
%     xsol = linspace(0,1,100);
%     usol = -4*umax.*(xsol.^2-xsol);
%     %%
%     figure(4)
%     plot(xx(ceil(end/2),:),u(ceil(end/2),:), xsol, usol);
%     legend('Numerical','Analytical')
%     xlabel('Y');
%     ylabel('U Velocity');
%     title('Vertical Centerline Profile');
%     grid on
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 4 4])
%     print( gcf, '-dpng', [nprefix 'UVel'])
%     
%     figure(5)
%     plot(yy(:,ceil(end/2)),v(:,ceil(end/2)))
%     xlabel('Y');
%     ylabel('V-Velocity');
%     title('Horizontal Centerline Profile');
%     grid on
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 4 4])
%     print( gcf, '-dpng', [nprefix 'VVel'])
    
    
    figure(6)
    [c,h] = contour(xxg,yyg,G',[0 0],'k','LineWidth',1);
    %clabel(c,h)1
%     hold on
%     [c,h] = contour(xxg,yyg,G1',[0 0],'k','LineWidth',2);
    %clabel(c,h)
    hold off
    xlabel('X');
    ylabel('Y');
    axis equal
    grid on
%     axis([0.2 0.8 0.5 1])
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 6 4])
    print( gcf, '-dpng', [nprefix 'Cont'])
%     
%     figure(7)
%     pcolor(xxg-dx/2,yyg-dy/2,alpha');
%     colorbar
%     colormap(flipud(bone))
    figure(8)
    pcolor(xxg-dx/2,yyg-dy/2,alpha');
    colorbar
    colormap(flipud(bone))
    axis equal
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 6 4])
    print( gcf, '-dpng', [nprefix 'alpha'])
    
%     figure(6)
%     [c,h] = contour(xxg,yyg,G',[-1:.1:1],'k','LineWidth',1);
%     %clabel(c,h)1
%     hold on
%     [c,h] = contour(xxg,yyg,G1',[0 0],'k','LineWidth',2);
%     %clabel(c,h)
%     hold off
%     xlabel('X');
%     ylabel('Y');
%     axis equal
%     grid on
% %     axis([0.2 0.8 0.5 1])
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 6 4])
%     print( gcf, '-dpng', [nprefix 'Cont'])
    
    
%     
%     figure(8)
%     [c,h] = contour(xxg,yyg,Ge2',[-1:.05:.5],'k');
%     xlabel('X');
%     ylabel('Y');
%     set(gcf, 'PaperPositionMode', 'manual');
%     set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 6 4])
%     print( gcf, '-dpng', [nprefix 'DCont'])

% figure(9)
% surf(xxg,yyg,G')

% figure(10)
% hold on
% [c,h] = contour(xxg,yyg,G1',[0 0],'k');
% [c,h] = contour(xxg,yyg,Ge0',[0 0],'k');
% % clabel(c,h);
% % subplot 222
% [c,h] = contour(xxg,yyg,Ge1',[0 0],'k');
% % clabel(c,h);
% % subplot 223
% [c,h] = contour(xxg,yyg,Ge2',[0 0],'k');
% % clabel(c,h);
% % subplot 224
% [c,h] = contour(xxg,yyg,Ge3',[0 0],'k');
% % clabel(c,h);
% hold off
% axis equal
% grid on
% % axis([0.2 0.8 0.5 1])
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 6 4])
% print( gcf, '-dpng', [nprefix 'Contf'])
% figure(10)
% image([xxg(1,1), xxg(end,end)],[yyg(1,1) yyg(end,end)],Markers'*64/3);
% axis ij;
% %image(xxg-dx/2,yyg-dy/2,Markers')
% hold on;
%  plot(min(max(-.05,lsegx),1.05)',min(max(-0.05,lsegy),1.05)','-*','LineWidth',2)
%  set(gca,'LineStyleOrder', '-*|:|o')
%  set(gca,'xtick',linspace(0,1,nx))
%  set(gca,'ytick',linspace(0,1,ny))
%  grid on
%  axis equal
%  hold off;
%  colorbar
% contour(xxg,yyg,Markers',[0 0])
figure(7)
 hold on
 plot(min(max(0,lsegx),1)',min(max(0,lsegy),1)','-k','LineWidth',2);
 set(0,'DefaultAxesColorOrder',[0,0,0])
 grid on
 axis equal
 set(gcf, 'PaperPositionMode', 'manual');
 set(gcf,'PaperUnits','inches','PaperPositionMode','manual','PaperPosition',[0 0 6 4])
 print( gcf, '-dpng', [nprefix 'plic'])
 hold off;

    
! mv ./*.png ~/Documents/ASU/MAE598Interfaces/HW4/


