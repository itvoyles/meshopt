function springsystem2d_force_setup(imax,jmax)
global SSAxinv SSAx SSAyinv SSAy AA AAinv dxdfx dydfy

% slow for loop to setup A matrix
% cnt = 1;
% for j = 1:jmax
%     for i = 1:imax
%         if i > 1
%           A(cnt,cnt-1) = -1;
%         end
%         
%         if j > 1
%             A(cnt,cnt-jmax) = -1;
%         end
%         
%         if i < imax
%             A(cnt,cnt+1) = -1;
%         end
%         
%         if j < jmax
%             A(cnt,cnt+imax) = -1;
%         end
%         
%         A(cnt,cnt) = -sum(A(cnt,:));
%         cnt = cnt + 1;
%     end
% end


% quick A matrix setup compared to for loop (above)
D = -ones(imax*jmax,1);
D1=D; D1(1:imax) = 0; A = diag(D1(imax+1:end),-jmax);
D1=D; D1(imax*(jmax-1)+1:imax*jmax)=0; A = A + diag(D1(1:end-imax),jmax);
D1=D; D1(1:imax:imax*jmax)=0; A = A + diag(D1(2:end),-1);
D1=D; D1(imax:imax:imax*jmax)=0; A = A + diag(D1(1:end-1),1);
A = A - diag(sum(A,2));


Ax = A;
Ay = A;


%setup x boundary conditions
Ax(1:imax:imax*jmax,:) = 0;
for i = 1:imax:imax*jmax
    Ax(i,i) = 1;
end

Ax(imax:imax:imax*jmax,:) = 0;
for i = imax:imax:imax*jmax
    Ax(i,i) = 1;
end


%setup y boundary conditions
Ay(1:imax,:) = 0;
for i = 1:imax
    Ay(i,i) = 1;
end

Ay(imax*(jmax-1)+1:imax*jmax,:) = 0;
for i = imax*(jmax-1)+1:imax*jmax
    Ay(i,i) = 1;
end
Ax(imax,:) = 0;
Ax(imax,imax) = 1;

%full matrix
% SSAfull = [Ax,zeros(size(Ax)); zeros(size(Ax)),Ay];
% SSAfullinv = SSAfull^-1;
% SSAfull = sparse(SSAfull);
% SSAfullinv = sparse(SSAfullinv);


SSAx = sparse(Ax);
SSAy = sparse(Ay);


AA = [Ax zeros(size(Ax,1),size(Ay,2))
      zeros(size(Ax,2),size(Ay,1)) Ay];
AAinv = inv(AA);

dxdfx = AAinv(1:imax*jmax,1:imax*jmax);
dydfy = AAinv(imax*jmax+1:end,imax*jmax+1:end);

dxdfx(1:imax:imax*jmax,:)=0;
dxdfx(imax:imax:imax*jmax,:)=0;
dxdfx(:,1:imax:imax*jmax)=0;
dxdfx(:,imax:imax:imax*jmax)=0;

dydfy(1:imax,:)=0;
dydfy(imax*(jmax-1)+1:imax*jmax,:)=0;
dydfy(:,1:imax)=0;
dydfy(:,imax*(jmax-1)+1:imax*jmax)=0;




%invert
% SSAxinv = invertmatrix(SSAx);
% SSAyinv = SSAxinv;
SSAxinv = SSAx^-1;
SSAyinv = SSAy^-1;