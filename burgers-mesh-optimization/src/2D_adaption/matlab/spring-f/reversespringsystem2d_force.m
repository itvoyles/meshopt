function F = reversespringsystem2d_force(x,y)
global SSAx SSAy imax jmax

Fx = SSAx*reshape(x,imax*jmax,1);
Fy = SSAy*reshape(y,imax*jmax,1);
F = [Fx;Fy];
