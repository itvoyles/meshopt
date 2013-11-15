function uyy = fd2d_d2udy2(f,dx)
  uyy = (f(2:end-1,3:end)-2*f(2:end-1,2:end-1)+f(2:end-1,1:end-2))/(dx^2);
