 function [xmat, ymat] = springsystem2d_coupled(kin)
global imax jmax Lref east west north south ixrowsleft ixrowsright ixrowsbottom...
       ixrowstop ixelemleft ixelemright ixelembottom ixelemtop xgs ygs solver

%EXAMPLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imax = 5;
% jmax = 5;
% Lref = 8;
% xx1 = repmat(linspace(0,pi,imax-1)',1,jmax);
% yy1 = repmat(linspace(0,pi,jmax),imax-1,1);
% xx2 = repmat(linspace(0,pi/2,imax)',1,jmax-1);
% yy2 = repmat(linspace(0,pi,jmax-1),imax,1);
% kk = 10*sin(xx1)+10*sin(yy1)+1;
% jj = 10*sin(yy2)+10*sin(xx2)+1;
% kin = [reshape(kk,(imax-1)*jmax,1);reshape(jj,(jmax-1)*imax,1)];
% % kin = ones(length(kin),1);
% 
% 
% [east west north south ixrowsleft ixrowsright ixrowsbottom ixrowstop ...
%        ixelemleft ixelemright ixelembottom ixelemtop] = setup_springsystem(imax,jmax);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nk = jmax*(imax-1);
nj = imax*(jmax-1);
n = imax*jmax;

jmat = kin(nk+1:nk+nj);
kmat = kin(1:nk);
kmat(imax*jmax) = 0;
jmat(imax*jmax) = 0;


A1 = -diag(kmat(east(1:end-1)),1);
A2 = -diag(kmat(west(2:end)),-1);
A3 = -diag(jmat(north(1:end-imax)),imax);
A4 = -diag(jmat(south(imax+1:end)),-imax);
AA = diag(kmat(east)+kmat(west)+jmat(north)+jmat(south));

A = A1+A2+A3+A4+AA;


%% setup coefficient matrix for x coordinates
Ax = A;
bx = zeros(imax*jmax,1);

%left bc
Ax(ixrowsleft,:) = 0;
Ax(ixelemleft) = 1;
bx(ixrowsleft,1) = -Lref/2;    

%right bc
Ax(ixrowsright,:) = 0;
Ax(ixelemright) = 1;
bx(ixrowsright,1) = Lref/2;




%% setup coefficient matrix for x coordinates
Ay = A;
by = zeros(imax*jmax,1);

%bottom bc
Ay(ixrowsbottom,:) = 0;
Ay(ixelembottom) = 1;
by(ixrowsbottom,1) = -Lref/2;    

%top bc
Ay(ixrowstop,:) = 0;
Ay(ixelemtop) = 1;
by(ixrowstop,1) = Lref/2;
%bottom bc









if solver == 1
    x = Ax\bx;
    y = Ay\by;
elseif solver == 2
    x = pentasolve(Ax,bx,imax); %faster for imax = jmax > 41
    y = pentasolve(Ay,by,imax); %faster for imax = jmax > 41
elseif solver == 3
    x=gauss_seidel(Ax,bx,xgs,imax);
    y=gauss_seidel(Ay,by,ygs,imax);
end

xmat = reshape(x,imax,jmax);
ymat = reshape(y,imax,jmax);

xgs = reshape(x,imax*jmax,1);
ygs = reshape(y,imax*jmax,1);


% 
% %number of elements
% nk = jmax*(imax-1);
% nj = imax*(jmax-1);
% n = imax*jmax;
% 
% %initialize elements
% kmat = zeros(imax+1,jmax);
% jmat = zeros(imax,jmax+1);
% Atemp = zeros((imax+2)*(jmax+2));
% 
% %setup matrices and relative spring locations
% kmat(2:imax,1:jmax) = reshape(k(1:nk),imax-1,jmax);
% jmat(1:imax,2:jmax) = reshape(k(nk+1:nk+nj),imax,jmax-1);
% kim1 = kmat(1:imax,1:jmax);
% kip1 = kmat(2:imax+1,1:jmax);
% kjm1 = jmat(1:imax,1:jmax);
% kjp1 = jmat(1:imax,2:jmax+1);
% 
% 
% 
% %setup primary A matrix (treats every point as an interior point)
% s = (imax+1)+2;
% for j = 1:jmax
%     for i = 1:imax
%         Atemp(s,s-(imax)) = -kjm1(i,j);
%         Atemp(s,s-1) = -kim1(i,j);
%         Atemp(s,s) = kim1(i,j)+kip1(i,j)+kjm1(i,j)+kjp1(i,j);
%         Atemp(s,s+1) = -kip1(i,j);
%         Atemp(s,s+(imax)) = -kjp1(i,j);
%         s = s+1;
%     end
% end
% 
% %extracts actual A matrix from the Atemp
% sx = (imax+1)+1;
% A = Atemp(sx+1:sx+n,sx+1:sx+n);
% 
% %setup coefficient matrix for x coordinates
% Ax = A;
% ileft = 1:imax:n;
% iright = imax:imax:n;
% Ax(ileft,:) = 0;
% Ax(iright,:) = 0;
% 
% Adiag = zeros(size(Ax,1),1);
% Adiag(ileft) = 1;
% Adiag(iright) = 1;
% Axtemp = diag(Adiag);
% Ax = Ax+Axtemp;
% 
% bx = zeros(size(Ax,1),1);
% bx(ileft) = -Lref/2;
% bx(iright) = Lref/2;
% 
% x = Ax\bx;
% xmat = reshape(x,imax,jmax);
% 
% %setup coefficient matrix for y coordinates
% Ay = A;
% jbottom = 1:imax;
% jtop = (imax)*(jmax-1)+1:n;
% Ay(jbottom,:) = 0;
% Ay(jtop,:) = 0;
% 
% Adiag = zeros(size(Ay,1),1);
% Adiag(jbottom) = 1;
% Adiag(jtop) = 1;
% Aytemp = diag(Adiag);
% Ay = Ay+Aytemp;
% 
% by = zeros(size(Ax,1),1);
% by(jbottom) = -Lref/2;
% by(jtop) = Lref/2;
% 
% y = Ay\by;
% ymat = reshape(y,imax,jmax);
% 
% 
% % Z = zeros(size(Ax));
% % AA = [Ax,Z; Z,Ay];
% % B = [bx;by];
% % XX = AA\B;
% % xmat = reshape(XX(1:n),imax,jmax);
% % ymat = reshape(XX(n+1:2*n),imax,jmax);
% 
% % mesh(xmat,ymat,xmat*0)
% % view([0,0,90])
% % axis equal
% 
% 
% 
% 
% % 