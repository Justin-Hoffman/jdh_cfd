function dat2 = frdrint(name,nx,ny)
    fid = fopen(name);
    dat1 = fread(fid,nx*ny, 'int32');
    size(dat1);
    for i = 1:nx
       dat2(i,:) = dat1(int32((i-1)*ny+1):int32((i-1)*ny+ny));
    end
end