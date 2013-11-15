function [ lim_va ] = limiter_vanalbada( neq,r )
%Function "limiter_vanalbada"
%   Adapted for MATLAB from Joe Derlaga's "limiter_vanalbada.f90", 071213.
%   Van Albada limiter.

lim_va = zeros(neq,1);

for eq = 1:neq,
    lim_va(eq) = (r(eq).*r(eq) + r(eq)) ./ (r(eq).*r(eq) + 1);
    if (r(eq) < 0),
        lim_va(eq) = 0;
    end
    
end

end

