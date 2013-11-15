function [ lim_v ] = limiter_vanleer( neq,r  )
%Function "limiter_vanleer"
%   Adapted for MATLAB from Jow Derlaga's "limiter_vanleer.f90", 071213.
%   Van Leer limiter.

lim_v = zeros(neq,1);

for eq = 1:neq,
    lim_v(eq) = (r(eq) + abs(r(eq))) ./ (1 + abs(r(eq)));
end


end

