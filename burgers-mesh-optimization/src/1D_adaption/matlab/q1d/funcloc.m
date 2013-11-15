function [ out ] = funcloc( Ms,pb,p01 )
%Function "funcloc"
%   Function to be driven to zero in order to solve for the 
%     shock location. Compares the pre-shock area ratio ( As/A1* )
%     calculated in two ways.
%          (As/A1*)_guess - (As/A1*)_calc = 0
%     Uses the pre-shock Mach number (Ms) as the independent variable.

global astar gm1 gp1 gmma rightbndry

% Ms = abs(Ms);
AsAe = astar(1)*( (1./Ms).*(2*(1+(Ms.^2)*gm1/2)/gp1).^(-1/2 + gmma/gm1) )/geom(rightbndry); % As/Ae 
As2As1 = (1 + 2*gmma*(Ms.^2 - 1)/gp1).^(1/gm1) * ( (2+gm1*Ms.^2)/(gp1*Ms.^2) ).^(gmma/gm1); % A2*/A1*

p02 = p01 * 1./As2As1; % Downstream stagnation pressure

Me = sqrt( (pb./p02).^(-gm1./gmma) - 1) * sqrt(2./gm1); % Exit Mach number
Me = abs(Me);
AeAs2 = (1./Me) .* (2*(1 + Me.^2 * gm1/2)/gp1).^(-1/2 + gmma/gm1); % Ae/A2*

out = (1./Ms).*(2*(1+(Ms.^2)*gm1/2)/gp1).^(-1/2 + gmma/gm1) ... % (As/A1*)_guess
    - AsAe .* As2As1 .* AeAs2; % (As/A1*)_calc
 

end

