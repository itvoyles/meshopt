function TE_firstorder = Burgers2d_TE1storder_terms(x,y,argin)

% imax = 33;
% jmax = 33;
% x = repmat(linspace(0,2,imax)',1,jmax);
% y = repmat(linspace(0,2,jmax),imax,1);
% Re = 8;
% Lref = 8;
% 
% theta = 45;
% Rot = [cosd(theta), -sind(theta); sind(theta),cosd(theta)];
% xnew = x.*Rot(1,1)+y.*Rot(1,2);
% ynew = x.*Rot(2,1)+y.*Rot(2,2);
% 
% x = xnew;
% y = ynew;

Re = argin(1);
Lref = argin(2);

nu = (2*Lref/Re);

[imax,jmax] = size(x);
% TE = Burgers2d_TE3(x,y,Re,Lref); %continuous residual
% uij = uxn(x,y,Re,Lref,0);

% i = 2:imax-1;
% j = 2:jmax-1;
% TE = uij(i,j).*((uij(i+1,j)-uij(i-1,j))./(x(i+1,j)-x(i-1,j))+(uij(i,j+1)-uij(i,j-1))./(y(i,j+1)-y(i,j-1)))...
%      -nu.*((uij(i+1,j)-2.*uij(i,j)+uij(i-1,j))./( ((x(i+1,j)-x(i-1,j))/2).^2)+(uij(i,j+1)-2.*uij(i,j)+uij(i,j-1))./(((y(i,j+1)-y(i,j-1)))/2).^2);




deta = 1;
dxsi = 1;


xxsi = (x(3:end,2:end-1)-x(1:end-2,2:end-1))/(2.*dxsi);
yxsi = (y(3:end,2:end-1)-y(1:end-2,2:end-1))/(2.*dxsi);
xeta = (x(2:end-1,3:end)-x(2:end-1,1:end-2))/(2.*deta);
yeta = (y(2:end-1,3:end)-y(2:end-1,1:end-2))/(2.*deta);
Ji = (xxsi.*yeta-xeta.*yxsi);

xsix = yeta./Ji;
etax = -yxsi./Ji;
xsiy = -xeta./Ji;
etay = xxsi./Ji;


x = x(2:end-1,2:end-1);
y = y(2:end-1,2:end-1);


uij = uxn(x,y,Re,Lref,0);
ux2y0 = uxn(x,y,Re,Lref,2);
ux3y0 = uxn(x,y,Re,Lref,3);
ux4y0 = uxn(x,y,Re,Lref,4);

ux0y2 = uxn(x,y,Re,Lref,2);
ux0y3 = uxn(x,y,Re,Lref,3);
ux0y4 = uxn(x,y,Re,Lref,4);

ux1y1 = uxn(x,y,Re,Lref,2);
ux2y1 = uxn(x,y,Re,Lref,3);
ux3y1 = uxn(x,y,Re,Lref,4);


ux1y2 = uxn(x,y,Re,Lref,3);
ux2y2 = uxn(x,y,Re,Lref,4);

ux1y3 = uxn(x,y,Re,Lref,4);



TE_firstorder = zeros(imax,jmax);
TE_firstorder(2:end-1,2:end-1) = deta.^2.*((uij.*ux3y0.*xeta.^3.*xxsi)./(6.*Ji) + (nu.*ux4y0.*xeta.^4.*xxsi.^2)./(4.*Ji.^2) + (uij.*ux2y1.*xeta.^2.*xxsi.*yeta)./(2.*Ji) + ...
           (2.*nu.*ux3y1.*xeta.^3.*xxsi.^2.*yeta)./(3.*Ji.^2) + (uij.*ux1y2.*xeta.*xxsi.*yeta.^2)./(2.*Ji) + (nu.*ux2y2.*xeta.^2.*xxsi.^2.*yeta.^2)./(2.*Ji.^2) + ...
           (uij.*ux0y3.*xxsi.*yeta.^3)./(6.*Ji) - (nu.*ux0y4.*xxsi.^2.*yeta.^4)./(12.*Ji.^2) - (uij.*ux3y0.*xeta.^3.*yxsi)./(6.*Ji) + (nu.*ux3y1.*xeta.^4.*xxsi.*yxsi)./(3.*Ji.^2) - ...
           (uij.*ux2y1.*xeta.^2.*yeta.*yxsi)./(2.*Ji) + (nu.*ux2y2.*xeta.^3.*xxsi.*yeta.*yxsi)./Ji.^2 + (nu.*ux4y0.*xeta.^3.*xxsi.*yeta.*yxsi)./(3.*Ji.^2) - ...
           (uij.*ux1y2.*xeta.*yeta.^2.*yxsi)./(2.*Ji) + (nu.*ux1y3.*xeta.^2.*xxsi.*yeta.^2.*yxsi)./Ji.^2 + (nu.*ux3y1.*xeta.^2.*xxsi.*yeta.^2.*yxsi)./Ji.^2 - ...
           (uij.*ux0y3.*yeta.^3.*yxsi)./(6.*Ji) + (nu.*ux0y4.*xeta.*xxsi.*yeta.^3.*yxsi)./(3.*Ji.^2) + (nu.*ux2y2.*xeta.*xxsi.*yeta.^3.*yxsi)./Ji.^2 + ...
           (nu.*ux1y3.*xxsi.*yeta.^4.*yxsi)./(3.*Ji.^2) - (nu.*ux4y0.*xeta.^4.*yxsi.^2)./(12.*Ji.^2) + (nu.*ux2y2.*xeta.^2.*yeta.^2.*yxsi.^2)./(2.*Ji.^2) + ...
           (2.*nu.*ux1y3.*xeta.*yeta.^3.*yxsi.^2)./(3.*Ji.^2) + (nu.*ux0y4.*yeta.^4.*yxsi.^2)./(4.*Ji.^2)) + ...
        dxsi.^2.*(-(uij.*ux3y0.*xeta.*xxsi.^3)./(6.*Ji) + (nu.*ux4y0.*xeta.^2.*xxsi.^4)./(4.*Ji.^2) + (uij.*ux3y0.*xxsi.^3.*yeta)./(6.*Ji) + (nu.*ux3y1.*xeta.*xxsi.^4.*yeta)./(3.*Ji.^2) - ...
           (nu.*ux4y0.*xxsi.^4.*yeta.^2)./(12.*Ji.^2) - (uij.*ux2y1.*xeta.*xxsi.^2.*yxsi)./(2.*Ji) + (2.*nu.*ux3y1.*xeta.^2.*xxsi.^3.*yxsi)./(3.*Ji.^2) + ...
           (uij.*ux2y1.*xxsi.^2.*yeta.*yxsi)./(2.*Ji) + (nu.*ux2y2.*xeta.*xxsi.^3.*yeta.*yxsi)./Ji.^2 + (nu.*ux4y0.*xeta.*xxsi.^3.*yeta.*yxsi)./(3.*Ji.^2) - ...
           (uij.*ux1y2.*xeta.*xxsi.*yxsi.^2)./(2.*Ji) + (nu.*ux2y2.*xeta.^2.*xxsi.^2.*yxsi.^2)./(2.*Ji.^2) + (uij.*ux1y2.*xxsi.*yeta.*yxsi.^2)./(2.*Ji) + ...
           (nu.*ux1y3.*xeta.*xxsi.^2.*yeta.*yxsi.^2)./Ji.^2 + (nu.*ux3y1.*xeta.*xxsi.^2.*yeta.*yxsi.^2)./Ji.^2 + (nu.*ux2y2.*xxsi.^2.*yeta.^2.*yxsi.^2)./(2.*Ji.^2) - ...
           (uij.*ux0y3.*xeta.*yxsi.^3)./(6.*Ji) + (uij.*ux0y3.*yeta.*yxsi.^3)./(6.*Ji) + (nu.*ux0y4.*xeta.*xxsi.*yeta.*yxsi.^3)./(3.*Ji.^2) + (nu.*ux2y2.*xeta.*xxsi.*yeta.*yxsi.^3)./Ji.^2 + ...
           (2.*nu.*ux1y3.*xxsi.*yeta.^2.*yxsi.^3)./(3.*Ji.^2) - (nu.*ux0y4.*xeta.^2.*yxsi.^4)./(12.*Ji.^2) + (nu.*ux1y3.*xeta.*yeta.*yxsi.^4)./(3.*Ji.^2) + ...
           (nu.*ux0y4.*yeta.^2.*yxsi.^4)./(4.*Ji.^2));
       
       



% max(max(abs(TE_firstorder-TE)))