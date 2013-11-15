function runner 

global reynoldsnumber nu Lref alpha TEmax curiter curfeval fidhistory ...
       leftbndry rightbndry method freq_history forward_mapping TEtolerance
global k xthroat toler astar tzero clls imax eqn cell_vol
global pb pzero sol order flux_type lim_type kappa shloc As2 flag functional
global gmma gm1 gp1 xg xgm1 xgp1 xg2m1 gxgm1 gxg2m1 gm1xgp1 gp1xgm1 cv cp R

imax = 33;
xc = zeros(imax-1,1);
ii = 1:imax-1;
eqn='q1d';

pb = 297000.0; pzero=300000.0; tzero=300.0;sol='super';order=7;kappa=-1;flag=0;
 flux_type='vanleer';lim_type='vanleer'; shloc=-99.9; As2=-99.9;
 xthroat=0.0; readgrid=0; interp=1; Lbnd=-2; Rbnd=2;
 toler = 1E-12;
 astar=zeros(2,1);
 astar(:) = 0.2;
rightbndry = Rbnd;

  gmma = 1.4; % Ratio of specific heats
    R = 286.9; % Metric gas constant for air
    gm1  = gmma - 1;
    gp1  = gmma + 1;
    xg   = 1/gmma;
    xgm1 = 1/(gmma - 1);
    xgp1 = 1/(gmma + 1);
    xg2m1 = 1/(gmma*gmma - 1);
    gxgm1 = gmma / (gmma - 1);
    gxg2m1 = gmma / (gmma*gmma - 1);
    gm1xgp1 = gm1/gp1;
    gp1xgm1 = gp1/gm1;
    
    cv = R*xgm1;
    cp = R*gxgm1;

    [xs As2] = shockloc(pb,pzero);
 keyboard
%  if readgrid == 1
%     fprintf('Reading grid...');
%     x0 = plot3d_read('grid.grd');
%     if length(x0)~=imax && interp==0
%         warning('Grid size does not match specified imax variable. Adjusting imax accordingly.')
%     else
%         fprintf('Done\n');
%     end
%     imax = length(x0);
% else
%     fprintf('Created evenly spaced grid.\n');
%     
%     if strcmp(eqn,'burgers')
%         x0 = linspace(-Lref/2,Lref/2,imax);
%     elseif strcmp(eqn,'q1d')
%         x0 = linspace(Lbnd,Rbnd,imax);
%     end
%     
% end
%     x0 = reshape(x0,length(x0),1);
%     
%     [ rho,u,p,M ] = super( x0,3,3,pzero,3 );
% 
% 
%     dx(ii) = (x0(ii+1) - x0(ii) ); % Cell sizes
%     %     cell_vol = [geom(xc) .* dx',geom(xc) .* dx',geom(xc) .* dx'];
%     [ solu ] = exactsol(x0,imax,sol,pb,pzero,order,shloc,flag,As2 );
%     keyboard
%     save('exact.dat','solu','-ascii','-double');
%     prim_cc = [solu(1,:);solu;solu(end,:)];
%     psi = tzero/(tzero - (gm1*(prim_cc(1,2).^2)./(2*gmma*R)));
%     
%     prim_cc(1,1) = pzero/(R*tzero.*psi.^xgm1);
%     
%     prim_cc(1,3) = pzero/(psi.^gxgm1);
%     
%     save prim.mat prim_cc
%     
%     [ TE ] = q1d_TE( imax,x0,pb,pzero,sol,order,kappa,flag,flux_type,lim_type,shloc,As2 );
%     
%     TE(2:imax,:) = TE(2:imax,:)./cell_vol(2:imax,:);
%     xc(ii) = 0.5*(x0(ii+1) + x0(ii) ); % Cell-center locations
%     TEcomp=[xc, TE(2:imax,:)];
%     save('TEcomp.dat','TEcomp','-ascii','-double')
% %     keyboard
%     
%     % Write restart file in the style of the Q1D code
%     
%     fidrst=fopen('q1d.rst','w');
%     
%     fprintf(fidrst, '%7.6f %7.6f %7.6f\n', 1.0,1.0,1.0);
%     fprintf(fidrst, '%7.6f %7.6f %7.6f\n', 1.0,1.0,1.0);
%     fprintf(fidrst, '%7.6f %7.6f %7.6f\n', 1.0,1.0,1.0);
%     
%     for cell = 1:imax+1,
%         fprintf(fidrst, '%23.15f %23.15f %23.15f\n', prim_cc(cell,1),prim_cc(cell,2),prim_cc(cell,3));
%     end
%     
%     for cell = 1:imax,
%         fprintf(fidrst,'%7.6f %7.6f %7.6f %7.6f %7.6f %7.6f\n',1.0,1.0,1.0,1.0,1.0,1.0);
%     end
%     
%     fclose(fidrst);
        
    
    
end
    