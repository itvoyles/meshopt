function [ar] = geom(x)
%Function "geom"
%   Contains the geometry of the CD nozzle, for use with the 
%   Q1D Euler equations.

n = length(x);
ar = zeros(n,1);

for i = 1:n;
    
    if abs(x(i)) <= 1,
        
        ar(i) = 0.2 + 0.4*(1 + sin(pi*(x(i) - 0.5))); %Consider function handle call
        
    elseif abs(x(i)) > 1,
        
        ar(i) = 1.0;
        
    else
%         disp('Outside of the domain')
        ar(i) = NaN;
    end
    
end

end

