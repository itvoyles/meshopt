function [x y] = springsystem2d_force(F)
global SSAxinv SSAyinv imax jmax SSAfullinv
%spring system with force as the design variable
 
% imax = 5;
% jmax = 5;
% springsystem2d_force_setup(imax,jmax)
% 
% Fx = zeros(imax,jmax);
% Fx(1,:) = 0;
% Fx(imax,:) = 1;
% Fy = zeros(imax,jmax);
% Fy(:,1) = 0;
% Fy(:,jmax) = 1;
% F = [reshape(Fx,imax*jmax,1);reshape(Fy,imax*jmax,1)];

n = length(F);


% F(1) = 1;
% F(imax) =  1;

xvec = SSAxinv*F(1:n/2);
yvec = SSAyinv*F(n/2+1:n);


x = reshape(xvec,imax,jmax);
y = reshape(yvec,imax,jmax);


% dx = x(2:end)-x(1:end-1);
% xc = (x(2:end)+x(1:end-1))/2;
% plot(xc,dx,'-o')


