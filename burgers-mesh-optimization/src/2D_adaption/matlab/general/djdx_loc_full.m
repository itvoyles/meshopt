function djdx = djdx_loc_full(x,y,ii,jj,dx,ex,J)
%given a mesh, this function calculates the mesh sensitivity of the functional
%J by perturbing the mesh by dx. Richardson extrapolation is then used to 
%estimate the error. If the error is below the defined tolerance (err_lim) then
%djdx is returned, if not the required refinement factor is calculated to reach 
%this limit and the derivative is calculated again.

  xip1 = x;
  xim1 = x;

  %calc 1
  xip1(ii,jj) = x(ii,jj)+dx;
  Jip1 = calc_j(xip1,y,ex);
  
  xim1(ii,jj) = x(ii,jj)-dx;
  Jim1 = calc_j(xim1,y,ex);
  
%   djdx = (Jip1-J)/(dx); % forward difference
  
  djdx = (Jip1-Jim1)/(2.*dx);
  
 
  

end
