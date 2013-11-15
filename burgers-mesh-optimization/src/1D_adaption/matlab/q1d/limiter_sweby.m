function [ lim_sw ] = limiter_sweby( neq,r,beta )
%Function "limiter_sweby"
%   Adapted for MATLAB from Joe Derlaga's "limiter_sweby.f90" 071213.
%   Sweby limiter.

lim_sw = zeros(neq,1);

for i = 1:neq,
    
     lim_sw(i) = max( 0.0, max( min(beta.*r(i), 1), min(r(i), beta) ) );
    
end


end

