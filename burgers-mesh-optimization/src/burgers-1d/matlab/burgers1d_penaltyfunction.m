function [phi] = burgers1d_penaltyfunction(x)
global nu alpha Lref Reynoldsnumber
% 
% Lref = 8;
% Reynoldsnumber = 16;
% alpha = Reynoldsnumber/2;
% nu = 2*Lref/Reynoldsnumber;
% 
% 
% x = linspace(-Lref/2,Lref/2,17);

dudx = abs(uxn(x,nu,alpha,Lref,1));

phi = (max(dudx)-dudx)/max(dudx);

