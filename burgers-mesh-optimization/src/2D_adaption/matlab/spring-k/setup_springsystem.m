function [east west north south ixrowsleft ixrowsright ixrowsbottom ixrowstop ...
       Ixelemleft Ixelemright Ixelembottom Ixelemtop] = setup_springsystem(imax,jmax)
   
% imax = 5;
% jmax = 5;

zero = imax*jmax;

ecnt = 1;
wcnt = 1;
ncnt = 1;
scnt = 1;
cnt =  1;
for jj = 1:jmax
    for ii = 1:imax

        if ii<imax
            east(cnt) = ecnt;
            ecnt = ecnt + 1;
        else
            east(cnt) = zero;
        end
        
        if ii>1
            west(cnt) = wcnt;
            wcnt = wcnt + 1;
        else
            west(cnt) = zero;
        end
        
        if jj<jmax
            north(cnt) = ncnt;
            ncnt = ncnt + 1;
            
        else
            north(cnt) = zero;
        end
        
        if jj>1
            south(cnt) = scnt;
            scnt = scnt + 1;
        else
            south(cnt) = zero;
        end
        
        cnt = cnt + 1;
        
    end
end


%get indicies for left and right boundary conditions
ixrowsleft = [1:imax:imax*jmax];
temp = zeros(imax*jmax);
for jj = 1:length(ixrowsleft)
  ii = ixrowsleft(jj);
  temp(ii,ii) = 1;
end
Ixelemleft = find(temp==1);

ixrowsright = [imax:imax:imax*jmax];
temp = zeros(imax*jmax);
for jj = 1:length(ixrowsright)
  ii = ixrowsright(jj);
  temp(ii,ii) = 1;
end
Ixelemright = find(temp==1);


%get indicies for bottom and top boundary conditions
ixrowsbottom = [1:imax];
temp = zeros(imax*jmax);
for jj = 1:length(ixrowsbottom)
  ii = ixrowsbottom(jj);
  temp(ii,ii) = 1;
end
Ixelembottom = find(temp==1);

ixrowstop = [(imax-1)*jmax+1:imax*jmax];
temp = zeros(imax*jmax);
for jj = 1:length(ixrowstop)
  ii = ixrowstop(jj);
  temp(ii,ii) = 1;
end
Ixelemtop = find(temp==1);
