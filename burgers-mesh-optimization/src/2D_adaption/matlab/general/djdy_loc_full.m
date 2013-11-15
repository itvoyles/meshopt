function djdy = djdy_loc_full(x,y,ii,jj,dy,ex,J)
%given a mesh, this function calculates the mesh sensitivity of the functional
%J by perturbing the mesh by dx. Richardson extrapolation is then used to 
%estimate the error. If the error is below the defined tolerance (err_lim) then
%djdy is returned, if not the required refinement factor is calculated to reach 
%this limit and the derivative is calculated again.

  yip1 = y;
  yim1 = y;


  %calc 1
  yip1(ii,jj) = y(ii,jj)+dy;
  Jip1 = calc_j(x,yip1,ex);

  yim1(ii,jj) = y(ii,jj)-dy;
  Jim1 = calc_j(x,yim1,ex);

%   djdy = (Jip1-J)/(dy); % forward differnece
  djdy = (Jip1-Jim1)/(2.*dy);



end
