
clc
% clear all


imax = 25;
xx = linspace(0,pi,imax-1);
k = sin(xx)+2;

X = springsystem(k);


% dx = (X(2:end)-X(1:end-1))';
% E0 = sum(1/2.*k.*dx.^2);
E0 = sum(k);



K = reversespringsystem(X,E0);

% K = K*(min(k(1))/min(K(1)));

plot(xx,k,'-ok',xx,K,'-*k')


% function [E] = springenergy(x,k)


