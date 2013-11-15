function [ lim_os ] = limiter_ospre( neq,r )
%Function "limiter_ospre"
%   Adapted for MATLAB from Joe Derlaga's "limiter_ospre.f90", 071213.
%   Ospre limiter.

lim_os = zeros(neq,1);

 for eq = 1:neq,
    lim_os(eq) = 1.5*(r(eq).*r(eq) + r(eq)) ./ (r(eq).*r(eq)+ r(eq) + 1);
    
    if (r(eq) < 0),
       lim_os(eq) = 0;
    end
    
 end


end

