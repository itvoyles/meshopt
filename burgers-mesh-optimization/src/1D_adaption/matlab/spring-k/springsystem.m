function x = springsystem(k)
global Lref iter leftbndry rightbndry
%spring system
%Ax=b
%A=f(k)
% clear all
% clc
% imax = 25;
% 
% xx = linspace(0,pi,imax-1);
% k = sin(xx)+0.1;

imax = length(k)+1;

A = zeros(imax,imax);
b = zeros(imax,1);
for ii = 2:imax-1
   A(ii,ii-1) = -k(ii-1);
   A(ii,ii) = ( k(ii)+k(ii-1));
   A(ii,ii+1) = -k(ii);
   
   b(ii) = 0;
end
A(1,1) = 1;
b(1) = leftbndry;
A(imax,imax) = 1;
b(imax) =  rightbndry;


x = A\b;


