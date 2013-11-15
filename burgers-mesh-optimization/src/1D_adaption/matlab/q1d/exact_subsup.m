% Q1D exact solutions

sol = 'sub'; % Specify exact solution type: 'sub' -- Isentropic Subsonic-subsonic
             %                              'super' -- Isentropic Subsonic-supersonic
             %                              'shock' -- Internal shock
             
% Isentropic Sub-Sup
k = ceil( (imax+1)/2) ; % index of throat (symmetric nozzle)
toler = 1e-12;
tzero = 300.0; % Stagnation temperature. K 
pzero = 300000.0; % Stagnation pressure, Pa

pb = 297330; % For subsonic case: Back pressure, Pa

% Initialize
mach = zeros(imax,1);
aratio = zeros(imax,1);
phi = zeros(imax,1);
fcn = zeros(imax,1);
delfcn = zeros(imax,1);
delmach = zeros(imax,1);

p = zeros(imax,1);
rho = zeros(imax,1);
u = zeros(imax,1);

gmma = 1.4;
gp1 = gmma+1; gm1 = gmma-1;

if strcmp(sol,'super'),
    mach(1:k) = 0.2;
    mach(k+1:imax) = 10.0;
    astar = geom(0); % A_t = A*
    
elseif strcmp(sol,'sub'),
    mach(:) = 0.2;
    prat = pb/pzero;
    fcn2 = 1000.0; % Initialize large
        
    % Check exit pressure ratio
    if prat < 0.937 || prat > 1.0,
        disp(sprintf('Exit pressure ratio is outside of the isentropic subsonic range.'))
    end
        
    for j=1:500, % Iteratively obtain exit Mach number from pressure ratio
        
        fcn2 = ( 1 + (mach(end)^2) * (gm1/2))^(-gmma/gm1) - prat;     
        
        if abs(fcn2) < toler,
            break
        end
        
        delfcn2 = (-gmma/gm1) * 2*mach(end)*(gm1/2);
        delmach2 = -fcn2/delfcn2;
        mach(end) = mach(end) + delmach2;
        
    end
    hold = ( (1/mach(end) ) * ( ((2/gp1) * (1+(gm1/2)*mach(end)^2))^(gp1/(2*gm1))) );
    astar = geom(x(end)) / hold;
    
    keyboard % debug
else
    disp('Enter a correct exact solution type')
    return
end


% Area ratio
for i = 1:imax,
    aratio(i) = geom(x(i)) / astar;
end
keyboard % debug
% Solve for Mach numbers using Newton iteration
ii = 1:imax;

for j=1:500, % Sufficient for toler >= 1E-12. Does not converge for smaller tolerances.
    
    phi(ii) = (2/gp1)*(1 + (gm1/2)* mach(ii) .* mach(ii));
    fcn(ii) = (aratio(ii).^2) .* (mach(ii).^2) - phi(ii).^(gp1/gm1);
    
    % Check Tolerance
    if max(( abs(fcn) ) ) < toler,
        break
    end
    
    delfcn(ii) = 2*mach(ii) .* (phi(ii).^(2/gm1) - aratio(ii).^2);
    delmach(ii) = fcn(ii) ./ delfcn(ii);
    mach(ii) = mach(ii) + delmach(ii)/4;
    
end


% Primitive variables

phi = 1 + (gm1/2)*mach.^2;
p = pzero./ ( phi.^(gmma/gm1) );
rho = p ./ (286.9*tzero./phi);
u(ii) = mach .* sqrt(gmma*286.9*tzero./phi);

% Plotting

figure(1)
subplot(2,3,1), plot(mach,aratio)
subplot(2,3,2), plot(x,rho)
subplot(2,3,3), plot(x,u)
subplot(2,3,4), plot(x,p)
subplot(2,3,5), plot(x,mach)
