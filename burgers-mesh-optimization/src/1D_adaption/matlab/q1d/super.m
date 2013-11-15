function [ rho,u,p,M ] = super( x,~,~,p0,~ )
%Function "super"
%   Computes an exact solution at an x-location for a subsonic-supersonic 
%   isentropic flow through a Q1D nozzle. Defined by an A* value.
%   -------
%   Inputs: sub(x,~,~,p0,~)
%            x    :: Vector of x-locations where solution is desired
%            p0   :: Local stagnation pressure
%              The number of input arguments is retained generality. 
%   Outputs: [rho,u,p,M], the solution
%            rho  :: Density (kg/m3)
%            u    :: Velocity (m2/s)
%            p    :: Pressure (Pa)
%            M    :: Mach number
global k xthroat toler gmma gm1 gp1 astar tzero

n = length(x);

rho = zeros(n,1);
u = zeros(n,1);
p = zeros(n,1);
M = zeros(n,1);
% M(1:k) = 0.2;      % ] Initial guess for sub-sup solution
% M(k+1:end) = 10.0; % ]
psi = ones(n,1);
phi = ones(n,1);
fcn = ones(n,1);
delfcn = zeros(n,1);
aratio = zeros(n,1);

% Determine if domain contains the throat; initialize appropriately
h=find(x>=xthroat,1);

if numel(h)==0,
    M(:) = 0.2; % All x's are before throat-- subsonic initialization
%     h = n;
elseif h==1,
    M(:) = 5.0; % All x's are after the throat-- supersonic initialization
else 
    M(1:h-1) = 0.2;  % throat is within domain-- initialize subsonic
    M(h:end) = 5.0;  %      before the throat; supersonic, after.
end

aratio = geom(x) / astar(1) ;

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

count = 0;
psi = (2/gp1)*(1 + (gm1/2)* M .* M);
fcn = (aratio.^2) .* (M.^2) - psi.^(gp1/gm1);
% 
% while max( abs(fcn(1)) ) > toler && count < 1000;
%     
%     delfcn(1) = 2*M(1) * (psi(1)^(2/gm1) - aratio(1).^2);
%     M(1) = abs( M(1) + 0.6*(fcn(1)/delfcn(1)) );
%     psi(1) = (2/gp1)*(1 + (gm1/2)* M(1) * M(1));
%     fcn(1) = (aratio(1)^2) * (M(1)^2) - psi(1)^(gp1/gm1);
%     disp(count)
%     count = count + 1;
% end
% count = 0;


%~~~~~ Unrolling loops
% for i=2:h-1,
%~~~~~~~~~
    while max( abs(fcn) ) > toler && count < 1000;
        
        delfcn = 2*M .* (psi.^(2/gm1) - aratio.^2);
        M = M + 0.6*(fcn./delfcn);
        psi = (2/gp1)*(1 + (gm1/2)* M .* M);
        fcn = (aratio.^2) .* (M.^2) - psi.^(gp1/gm1);
        
        count = count + 1;
    end
%     M(i) = M(i-1);
%     psi(i) = (2/gp1)*(1 + (gm1/2)* M(i) * M(i));
%     fcn(i) = (aratio(i)^2) * (M(i)^2) - psi(i)^(gp1/gm1);  
%     
%     while max( abs(fcn(i)) ) > toler && count < 1000;
%         
%         delfcn(i) = 2*M(i) * (psi(i)^(2/gm1) - aratio(i)^2);
%         M(i) = abs( M(i) + 0.6*(fcn(i)/delfcn(i)) );
%         psi(i) = (2/gp1)*(1 + (gm1/2)* M(i) * M(i));
%         fcn(i) = (aratio(i)^2) * (M(i)^2) - psi(i)^(gp1/gm1);
%         disp(count)
%         count = count + 1;
%     end
%     count = 0;
%     %~~~~~~~~~
% end
% 
% if n > h-1,
% 
% %     M(h) = M(h-1) + 2.0;
% %     psi(h) = (2/gp1)*(1 + (gm1/2)* M(h) * M(h));
% %     fcn(h) = (aratio(h)^2) * (M(h)^2) - psi(h)^(gp1/gm1);
%     
%     for i=h:n,
%         
%         M(i) = abs(M(i-1));
%         psi(i) = (2/gp1)*(1 + (gm1/2)* M(i) * M(i));
%         fcn(i) = (aratio(i)^2) .* (M(i)^2) - psi(i)^(gp1/gm1);
%         
%         while max( abs(fcn(i)) ) > toler && count < 1000;
%             
%             delfcn(i) = 2*M(i) * (psi(i)^(2/gm1) - aratio(i)^2);
%             M(i) = abs( M(i) + 0.6*(fcn(i)/delfcn(i)) );
%             psi(i) = (2/gp1)*(1 + (gm1/2)* M(i) * M(i));
%             fcn(i) = (aratio(i)^2) * (M(i)^2) - psi(i)^(gp1/gm1);
%             disp(count)
%             count = count + 1;
%         end
%         count = 0;
%     end
%     
% end

%~~~~~~~~~

phi = 1 + (gm1/2)*M.^2;
p = p0 ./ ( phi.^(gmma/gm1) );
rho = p ./ (286.9*tzero./phi);
u = M .* sqrt(gmma*286.9*tzero./phi);

end

