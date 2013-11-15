function [ TE ] = Euler2D_TE( x,y,x_xsi,y_xsi,x_eta,y_eta,J,exact1,exact2,exact3,exact4,parameters,F1,G1,TEMatrix )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
exact1=exact1*0;
exact2=exact2*0;
exact3=exact3*0;
exact4=exact4*0;
totalTE = truncErrorU1(x,y,x_xsi,y_xsi,x_eta,y_eta,J,exact1,exact2,exact3,exact4,parameters,F1,G1,TEMatrix);
TE = TEMatrix;

end

