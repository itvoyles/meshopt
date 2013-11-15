function [ M2, T2, P2] = PM_exp_fan( M1, T1, P1, angle )
%Takes mach number, temperature and pressure of incoming flow and uses turn angle (degrees) to 
%   calculate the outgoing mach, temp and pressure

%convert angle to radians
angle = pi*angle/180;

%calculate prandtl-meyer function v(M1)
V1 = sqrt( (1.4+1)/(1.4-1) )*atan( sqrt( (1.4-1)*(M1^2-1)/(1.4+1) ) ) - atan( sqrt( M1^2-1 ) );

%Note angle = v(M2)-v(M1)
V2 = angle + V1;

%Solve for second mach number M2--assuming 1 < M2 < 30 
M2 = fzero(@(x)( (sqrt( (1.4+1)/(1.4-1) )*atan( sqrt( (1.4-1)*(x^2-1)/(1.4+1) ) ) - atan( sqrt(x^2-1) )) - V2 ),[1,30] );

%Compute second temperature
T2 = T1*( (1+(1.4-1)*M1^2/2)/(1+(1.4-1)*M2^2/2) );

%Compute second pressure
P2 = P1*( (1+(1.4-1)*M2^2/2)/(1+(1.4-1)*M1^2/2) )^(-1.4/(1.4-1));
end

