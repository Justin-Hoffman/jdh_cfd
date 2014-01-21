
    [nxo nyo] = size(do);
    [nxv nyv] = size(dv);
    ho = [];
    figure(1)
    close
    figure(2)
    close
    [xx yy] = meshgrid(linspace(0,1,nx),linspace(0,1,ny));
    
    figure(1)
    for j=1:2:nyo
       for i=1:nxo 
            if(isnan(do(i,j))==0)
               o(i,1:2) = do(i,j:j+1);
            end
       end
       seenan = 0;
       hold on
       plot(o(:,1),o(:,2),'.');
       clear o
    end
    contour(xx,yy,strm(2:nx+1,2:ny+1)',strmlevels)
    axis equal
    axis([0 1 0 1])
    xlabel('x');
    ylabel('y');
    title('Stream Function Contour Implicit Re=100');
    
    
    figure(2)
    for j=1:2:nyv
       for i=1:nxv 
            if(isnan(dv(i,j))==0)
               v(i,1:2) = dv(i,j:j+1);
            end
       end
       hold on
       plot(v(:,1),v(:,2),'.');
       clear v
    end
    contour(xx,yy,-omega',omegalevels)
        
    axis equal
    axis([0 1 0 1])
    xlabel('x');
    ylabel('y');
    title('Vorticity Contour Implicit Re=100');
    
