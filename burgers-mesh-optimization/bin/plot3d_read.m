function [x,y,z] = plot3d_read(filename)


data = importdata(filename);

nblocks = data(1);
imax = data(2);
jmax = data(3);
kmax = data(4);

nodes = imax*jmax*kmax;

if (length(data)~= 3*nblocks*nodes+4) 
    fprintf('Error! Incorrect plot3d file.')
    x = 0.0;
    y = 0.0;
    z = 0.0;
    return
end


xvec = data(5:5+nodes-1);
yvec = data(5+nodes:5+2*nodes-1);
zvec = data(5+2*nodes:5+3*nodes-1) ;

x = reshape(xvec,imax,jmax,kmax);
y = reshape(yvec,imax,jmax,kmax);
z = reshape(zvec,imax,jmax,kmax);

end