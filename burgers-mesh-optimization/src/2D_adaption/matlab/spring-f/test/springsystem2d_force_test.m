%% springsystem2d force test
global Lref imax jmax

clc
% clear all

imax = 33;
jmax = imax;
Lref = 8;
x = repmat(linspace(-4,4,imax)',1,jmax);
y = repmat(linspace(-4,4,jmax),jmax,1);

% x = repmat(linspace(atan(-4),atan(4),imax)',1,jmax);
% x = tan(x);
% y = x';

springsystem2d_force_setup(imax,jmax)

F = reversespringsystem2d_force(x,y);

xtest = springsystem2d_force(F);

disp(max(max(xtest-x)));
