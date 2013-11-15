%% springsystem force test
global Lref

clc
% clear all

n = 33;
Lref = 8;
% x = linspace(-4,4,n)';
x = linspace(atan(-4),atan(4),n)';
x = tan(x);

springsystem_force_setup(n)

F = reversespringsystem_force(x);

xtest = springsystem_force(F);

disp(max(xtest-x));
