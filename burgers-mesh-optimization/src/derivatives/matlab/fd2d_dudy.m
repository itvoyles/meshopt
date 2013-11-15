function uy = fd2d_dudy(f,dx)

  uy = (f(2:end-1,3:end)-f(2:end-1,1:end-2))/(2*dx);
