function ux = fd2d_dudx(f,dx)

ux = (f(3:end,2:end-1)-f(1:end-2,2:end-1))./(2*dx);
