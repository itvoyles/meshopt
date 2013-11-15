function [ rho,u,p,M ] = shock( x,shloc,pb,p01,f )
%Function "shock"
%   Computes an exact solution at an x-location for a Q1D nozzle
%   with an internal shock. The solution is decomposed into two
%   isentropic solutions, on either side of the shock.
%   -------
%   Inputs: sub(x,shloc,pb,p01,f)
%            x     :: Vector of x-locations where solution is desired
%            shloc :: Shock location (m)
%            pb    :: Nozzle back pressure
%            p01   :: Upstream stagnation pressure
%            f     :: Flagging variable, for selecting the A* value
%   Outputs: [rho,u,p,M], the solution
%            rho   :: Density (kg/m3)
%            u     :: Velocity (m2/s)
%            p     :: Pressure (Pa)
%            M     :: Mach number
global astar %toler gmma gm1 gp1

p02 = p01 * astar(1)/astar(2);

% Divide chosen domain:
g = find(x >= shloc,1); % Locates shock within vector of x-values

if numel(g) == 0, % Entire domain is pre-shock
    [rho,u,p,M] = super(x,shloc,pb,p01,f);
elseif g == 1, % Entire domain is post-shock
    [rho,u,p,M] = sub(x,shloc,pb,p02,f);
else           % Shock is within domain
    
    x1 = x(1:g-1); % Pre-shock
    x2 = x(g:end); % Post-shock
    
    
    % Call isentropic solutions:
    [rho1, u1, p1, M1] = super(x1,shloc,pb,p01,f);
    [rho2, u2, p2, M2] = sub(x2,shloc,pb,p02,f);
    
    % Concatenate solutions:
    rho = [rho1;rho2];
    u = [u1;u2];
    p = [p1;p2];
    M = [M1;M2];
end

end

