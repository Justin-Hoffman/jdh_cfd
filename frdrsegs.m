function [lsegx lsegy] = frdrsegs(name,nx,ny)
    fid = fopen(name);
    dat1 = fread(fid,(nx*ny*4), 'double');
    size(dat1);
    lsegx = [];
    lsegy = [
        ];
    for i = 1:nx
        for j = 1:ny
           lsegx = [lsegx; dat1(4*(nx*(i-1)+(j-1))+1)   dat1(4*(nx*(i-1)+(j-1))+3)];
           lsegy = [lsegy; dat1(4*(nx*(i-1)+(j-1))+2) dat1(4*(nx*(i-1)+(j-1))+4)];

        end
    end
end