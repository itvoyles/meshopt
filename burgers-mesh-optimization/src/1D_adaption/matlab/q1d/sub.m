function [ rho,u,p,M ] = sub( x,~,pb,p0,f )
%Function "sub"
%   Computes an exact solution at an x-location for a subsonic, 
%   isentropic flow through a Q1D nozzle.
%   -------
%   Inputs: sub(x,~,pb,p0,f)
%            x    :: Vector of x-locations where solution is desired
%            pb   :: Nozzle back pressure
%            p0   :: Local stagnation pressure
%            f    :: Flagging variable, for selecting the A* value
%   Outputs: [rho,u,p,M], the solution
%            rho  :: Density (kg/m3)
%            u    :: Velocity (m2/s)
%            p    :: Pressure (Pa)
%            M    :: Mach number
global toler gmma gm1 gp1 astar tzero

n = length(x);
pr = pb/p0; % Pressure ratio

rho = zeros(n,1);
u = zeros(n,1);
p = zeros(n,1);
M = zeros(n,1);
M(:) = 0.2; % Initial guess for subsonic solution
phi = zeros(n,1);
fcn = ones(n,1);
delfcn = zeros(n,1);
aratio = zeros(n,1);

h=1;
if f~=0,
    h = 2;
end

aratio = geom(x) / astar(h) ;

count = 0;
phi = (2/gp1)*(1 + (gm1/2)* M .* M);
fcn = (aratio.^2) .* (M.^2) - phi.^(gp1/gm1);

while max(abs(fcn)) > toler && count < 500,
    
    delfcn = 2*M .* (phi.^(2/gm1) - aratio.^2);
    delmach = fcn ./ delfcn;
    M = M + 0.25*delmach;
    
    phi = (2/gp1)*(1 + (gm1/2)* M .* M);
    fcn = (aratio.^2) .* (M.^2) - phi.^(gp1/gm1);
    
    count = count + 1;
end

% for j=1:500,
%     
%     phi = (2/gp1)*(1 + (gm1/2)* M .* M);
%     fcn = (aratio.^2) .* (M.^2) - phi.^(gp1/gm1);
%     
%     % Check Tolerance
%     if max(( abs(fcn) ) ) < toler,
%         break
%     end
%     
%     delfcn = 2*M .* (phi.^(2/gm1) - aratio.^2);
%     delmach = fcn ./ delfcn;
%     M = M + 0.25*delmach;
%     
% end

phi = 1 + (gm1/2)*M.^2;
p = p0 ./ ( phi.^(gmma/gm1) );
rho = p ./ (286.9*tzero./phi);
u = M .* sqrt(gmma*286.9*tzero./phi);

end

