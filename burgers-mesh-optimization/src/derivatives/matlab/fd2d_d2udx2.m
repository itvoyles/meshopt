function uxx = fd2d_d2udx2(f,dx)

  uxx = (f(3:end,2:end-1)-2*f(2:end-1,2:end-1)+f(1:end-2,2:end-1))/(dx^2);
