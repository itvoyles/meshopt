function k = reversespringsystem2d_coupled(x,y)
global imax jmax

% imax = 5;
% jmax = 5;
% x = repmat(linspace(0,1,imax)',1,jmax);
% y = repmat(linspace(0,1,jmax),imax,1);


nk = (imax-1)*jmax;
nj = (jmax-1)*imax;
nn = nk+nj;

dxi = zeros(imax+1,jmax);
dxj = zeros(imax,jmax+1);
dyi = dxi;
dyj = dxj;
% YY = XX;

% XX(2:imax+1,2:jmax+1) = x;
% YY(2:imax+1,2:jmax+1) = y;

dxi(2:imax,1:jmax) = x(2:end,1:jmax)-x(1:end-1,1:jmax);
dxj(1:imax,2:jmax) = x(1:imax,2:jmax)-x(1:imax,1:jmax-1);

dyi(2:imax,1:jmax) = y(2:end,1:jmax)-y(1:end-1,1:jmax);
dyj(1:imax,2:jmax) = y(1:imax,2:jmax)-y(1:imax,1:jmax-1);


% xi = XX(2:imax+1,2:jmax+1);
dxim1 = dxi(1:imax,1:jmax);
dxip1 = dxi(2:imax+1,1:jmax);
dxjm1 = dxj(1:imax,1:jmax);
dxjp1 = dxj(1:imax,2:jmax+1);

% yi = YY(2:imax+1,2:jmax+1);
dyim1 = dyi(1:imax,1:jmax);
dyip1 = dyi(2:imax+1,1:jmax);
dyjm1 = dyj(1:imax,1:jmax);
dyjp1 = dyj(1:imax,2:jmax+1);



Axtemp = zeros(nk,nn);
Aytemp = zeros(nj,nn);
bx = zeros(nk,1);
by = zeros(nj,1);


nkeq = (imax-2)*jmax;
njeq = (jmax-2)*imax;

sx = 1;
sy = nk+2;
%x force balance
kind = reshape([1:(imax-1)*(jmax)],imax-1,jmax);
jind = reshape([1:(imax)*(jmax-1)]+nk,imax,jmax-1);

kim1 = kind(1:imax-2,1:jmax);
kip1 = kind(2:imax-1,1:jmax);
jjm1 = jind(2:imax-1,1:jmax-1);
jjp1 = jind(2:imax-1,1:jmax-1);

sx = 1;
for j = 1:jmax
    
    Axtemp(sx,kind(1,j)) = 1;
    bx(sx) = 1;
    
    sx = sx + 1;
    for i = 2:imax-1
       Axtemp(sx,kim1(i-1,j)) = dxim1(i,j);
       Axtemp(sx,kip1(i-1,j)) = -dxip1(i,j);
       
      
       if j>1
        Axtemp(sx,jjm1(i-1,j-1)) = dxjm1(i,j);
       end
       
       if j<jmax
        Axtemp(sx,jjp1(i-1,j)) = -dxjp1(i,j);
       end

       sx = sx + 1;
       sy = sy + 1;
    end
end


%x force balance
kind = reshape([1:(imax-1)*(jmax)],imax-1,jmax);
jind = reshape([1:(imax)*(jmax-1)]+nk,imax,jmax-1);

kim1 = kind(1:imax-1,2:jmax-1);
kip1 = kind(1:imax-1,2:jmax-1);
jjm1 = jind(1:imax,1:jmax-2);
jjp1 = jind(1:imax,2:jmax-1);

sx = 1;


for j = 1:jmax-1
    for i = 1:imax
       if j == 1
            Aytemp(sx,jind(i,1)) = 1;
            by(sx) = 1;
            sx = sx + 1;
            
       else
       
           Aytemp(sx,jjp1(i,j-1)) = -dyjp1(i,j);
           Aytemp(sx,jjm1(i,j-1)) = dyjm1(i,j);

           if i>1
            Aytemp(sx,kim1(i-1,j-1)) = dyim1(i,j);
           end

           if i<imax
            Aytemp(sx,kip1(i,j-1)) = -dyip1(i,j);
           end

           sx = sx + 1;
           sy = sy + 1;
       
       end
    end
end



A = sparse([Axtemp;Aytemp]);
b = [bx;by];

    



k = A\b;

    
end
