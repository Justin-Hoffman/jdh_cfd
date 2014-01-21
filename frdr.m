function dat2 = frdr(name,nx,ny)
    fid = fopen(name);
    dat1 = fread(fid,nx*ny, 'double');
    size(dat1);
    for i = 1:nx
       dat2(i,:) = dat1(int32((i-1)*ny+1):int32((i-1)*ny+ny));
    end
end