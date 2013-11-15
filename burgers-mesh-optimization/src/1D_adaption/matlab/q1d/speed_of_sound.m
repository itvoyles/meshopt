function [ a ] = speed_of_sound( p, rho )
%Function "speed_of_sound"
%   MATLAB version of Joe Derlaga's "speed_of_sound.f90",
%   created 070313.
%   Takes a pressure and a density, returns the perfect gas
%   speed of sound.

global gmma 

a = sqrt(gmma*p/rho);

end

