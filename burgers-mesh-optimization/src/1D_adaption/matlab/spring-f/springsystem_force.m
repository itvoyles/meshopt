function x = springsystem_force(F)
global SSAinv
%spring system with force as the design variable
 

x = SSAinv*F;




