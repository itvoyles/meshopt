function [ solu,pdadxc ] = exactsol(x,param,sol,pb,pzero,order,shloc,flag,As2 )
% Function "exactsol"
% Chooses the type of exact solution, computes one-time argument values,
% and calls individual exact solution functions.

global xthroat toler astar tzero 
global gmma gm1 gp1 xg xgm1 xgp1 xg2m1 gxgm1 gxg2m1 gm1xgp1 gp1xgm1 cv cp


% ii=1:imax-1;
ii=1:length(x)-1;
% xc(ii) = 0.5*(x(ii+1) + x(ii) ); % Cell-center locations 
prat = pb/pzero; % Pressure ratio

% shloc = 10.0; % Shock location-- modified if needed
% flag = 0; % Used to indicate the that there is a shock in the solution,
          % for purposes of selecting the correct A* value.
          % Modified automatically.
          
% astar = zeros(2,1);
rhoc = zeros(param-1,1); 
uc = zeros(param-1,1);
pc = zeros(param-1,1);
machc = zeros(param-1,1);
pdadxc = zeros(param-1,1);
% Select exact solution:

if strcmp(sol,'sub'),
    
    Me = sqrt( prat^(-gm1/gmma) - 1) * sqrt(2/gm1); % Exit Mach
    astar(:) = geom(x(end)) / ( (1/Me ) * ( ((2/gp1) * (1+(gm1/2)*Me^2))^(gp1/(2*gm1))) );
    exsol = @sub;

elseif strcmp(sol,'super')
    
    astar(:) = geom(xthroat);
    exsol = @super;

elseif strcmp(sol,'shock')

%     flag = 1; % Activate use of two A* values
%     astar(:) = geom(xthroat); % A* for upstream, A* = A_t
%     
%     [ shloc, As2 ] = shockloc(pb,pzero); % astar(2) modified in this function
% 
     astar(2) = As2;
    
    exsol = @shock;

end

% Cell-averaged values:
% Scaled quadrature locations and weights:

% 7-point Gauss quadrature 
w =[0.4179591836734694; 0.3818300505051189; 0.3818300505051189; ...
    0.2797053914892766; 0.2797053914892766; 0.1294849661688697; ...
    0.1294849661688697];
t = [0.0000000000000000; 0.4058451513773972; -0.4058451513773972; ...
    -0.7415311855993945; 0.7415311855993945; -0.9491079123427585; ...
    0.9491079123427585];

%% General order Curtis Clenshaw
% [ t,w ] = curtis_clenshaw( order );

% [ t,I ] = sort(t); % Puts x-locations in order
% w = w(I);

% Cell limits:
a(ii) = x(ii); % left limits
b(ii) = x(ii+1); % right limits

for j=1:length(x)-1,
    
%     xq = (b(j) - a(j))*t + a(j); % x locations for Curtis-Clenshaw
    xq = 0.5*(b(j)-a(j))*t + 0.5*(b(j) + a(j)); % x locations for Gauss Quad   

    [rhoq,uq,pq,machq] = exsol(xq,shloc,pb,pzero,flag);
    [pdadxq] = pdadx( xq, pq );
    
    rhoc(j) = 0.5*sum(rhoq .* w);
    uc(j) = 0.5*sum(uq .* w);
    pc(j) = 0.5*sum(pq .* w);
    machc(j) = 0.5*sum(machq .* w);
    pdadxc(j) = 0.5*sum(pdadxq .* w);
end

% figure(1)
% subplot(2,2,1), plot(xc,rhoc)
% subplot(2,2,2), plot(xc,uc)
% subplot(2,2,3), plot(xc,pc)
% subplot(2,2,4), plot(xc,machc)

% To plot solution over the domain:

%[rho,u,p,mach] = exsol(x,shloc,pb,pzero,flag);

% figure(2)
% plot(x,rho,'-*',xc,rhoc,'-o')
% title('Shock Solution: Density, pb/pc = 0.99')
% legend('Exact (Faces)','Average (Centers)','Location','SouthWest')

% subplot(2,2,1), plot(x,rho,xc,rhoc)
% subplot(2,2,2), plot(x,u)
% subplot(2,2,3), plot(x,p)
% subplot(2,2,4), plot(x,mach)

solu = [rhoc, uc, pc];

end
