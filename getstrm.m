function strm = getstrm(omega)
    [nx ny] = size(omega);
    dx = 1/nx;
    dy = 1/ny;
    strm = zeros(nx+2,ny+2);
    for t = 1:50000
       for j=2:ny+1 
        for i=2:nx+1
           strm(i,j) = 1/4*(strm(i+1,j)+strm(i-1,j)+strm(i,j+1)+strm(i,j-1))+1/4*dx*dy*omega(i-1,j-1);
        end
       end
    end
end