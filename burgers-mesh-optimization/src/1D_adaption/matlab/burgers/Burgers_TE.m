function [TE]=Burgers_TE(x,x_xi,x_xixi,dxi,v,a,L,f)


  if nargin==8
    flag=f;
  else
    flag = 0;
  end

  u = uxn(x,v,a,L,0);
  uxx = uxn(x,v,a,L,2);
  uxxx = uxn(x,v,a,L,3);
  uxxxx = uxn(x,v,a,L,4);
  
%   !standard term
  TE_stand = dxi.^2.*x_xi.^2.*(...
                              1./6.*u.*uxxx    ...
                             -1./12.*v.*uxxxx  ...
                              );

  TE_stretch = dxi.^2.*x_xixi.*(                                             ...
                               -1./3.*v.*uxxx                     ...
                               +0.5.*u.*uxx                             ...
                               );
 
%   !mixed term
  TE_mixed = dxi.^2.*v.*(x_xixi./x_xi).^2.*uxx/4.;


  if (flag == 1)
    TE = TE_stand;
  elseif (flag == 2)
    TE = TE_stretch;
  elseif (flag == 3)
    TE = TE_mixed;
  else
    TE = TE_stand+TE_stretch+TE_mixed;
  end


  end
