function plot3d_write(filename,x,y,z)

imax = size(x,1);
jmax = size(x,2);
kmax = size(x,3);

nodes = imax*jmax*kmax;
nblocks = size(x,4);

fid = fopen(filename,'w');

fprintf(fid,'%3.0f\n',nblocks);
fprintf(fid,'%4.0f %4.0f %4.0f \n',imax,jmax,kmax);


xvec = reshape(x,nodes,1);
yvec = reshape(y,nodes,1);
zvec = reshape(z,nodes,1);


%write x
for i = 1:nodes
    newline=1;%flag to skip a line
    
    fprintf(fid,'%23.15e ',xvec(i));
    if mod(i,4)==0
      fprintf(fid,'\n');
      newline=0;
   end
end

if newline==1
  fprintf(fid,'\n');
end

%write y
for i = 1:nodes
    newline=1;%flag to skip a line
    
    fprintf(fid,'%23.15e ',yvec(i));
    if mod(i,4)==0
      fprintf(fid,'\n');
      newline=0;
   end
end

if newline==1
  fprintf(fid,'\n');
end


%write z
for i = 1:nodes
    newline=1;%flag to skip a line
    
    fprintf(fid,'%23.15e ',zvec(i));
    if mod(i,4)==0
      fprintf(fid,'\n');
      newline=0;
   end
end

if newline==1
  fprintf(fid,'\n');
end


end