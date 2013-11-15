function uxy = fd2d_d2udxy(u,dxi,deta)

  uxy = (u(3:end,3:end)-u(3:end,1:end-2)-u(1:end-2,3:end)+u(1:end-2,1:end-2))/(4*dxi*deta);
