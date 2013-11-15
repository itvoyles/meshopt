function [ TE ] = q1d_TE( parameter,x,pb,pzero,sol,order,kappa,flag,flux_type,lim_type,shloc,As2 )
%Function "q1d_TE"
%   Adapted from Joe Derlaga's "residual.f90", 071013.
%   Takes an FVM exact solution, flux, limiter, and MUSCL designations
%   and calulates the truncation error over the domain.
global cell_vol tzero gmma gm1 xgm1 gxgm1 R

ii = 1:parameter-1;

xc(ii) = 0.5*(x(ii+1) + x(ii) ); % Cell-center locations
xc = [xc(1),xc,xc(end)]; % For extended nozzle, ghost cell areas equal to 
                     % areas at inlet and outlet respectively; so,
                     % repeat inlet and outlet locations 
dx(ii) = (x(ii+1) - x(ii) ); % Cell sizes
dx = [dx(1),dx,dx(end)]; % so that entry (1) corresponds to a ghost cell
cell_vol = [geom(xc) .* dx',geom(xc) .* dx',geom(xc) .* dx'];
% dadx_cc = zeros(imax+1,1);
dadx_cc = dadx( xc ); % Area derivatives

% Exact solution evaluation

[ solu, pdadx_cc ] = exactsol( x,parameter,sol,pb,pzero,order,shloc,flag,As2 );

% Modify matrix of primitive variables to include ghost cells.
% This will allow for consistent indexing within the MATLAB functions.
% For consistency with the Fortran Q1D code, velocities in the first
% interior cells are copied to the ghost cells; then, density and pressure
% are calculated from isentropic relations in the inlet.

prim_cc = [solu(1,:);solu;solu(end,:)];

psi = tzero/(tzero - (gm1*(prim_cc(1,2).^2)./(2*gmma*R)));

prim_cc(1,1) = pzero/(R*tzero.*psi.^xgm1);

prim_cc(1,3) = pzero/(psi.^gxgm1);

TE = zeros(parameter+1,3); % Initialize: # of cells + 2 ghost cells
F = zeros(parameter,3); % # of faces
% S = zeros(imax+1,3); % # of cells + 2 ghost cells
% keyboard

[F] = construct_flux( prim_cc, kappa, flux_type, lim_type ); % Call flux routine
[S] = construct_source( prim_cc, dadx_cc ); % Call area source term routine
% load prim_diff.mat
% prim_diff = abs(prim_diff);
% [F] = construct_flux( prim_diff, kappa, flux_type, lim_type ); % Call flux routine
% [S] = construct_source( prim_diff, dadx_cc ); % Call area source term routine
% keyboard

% [S] = pdadx_cc; % For a quadrature of pdadx over cells

% Because MATLAB does not allow an index of zero:
%    +----+----+--]  [--+----+----+     
%    | G  |    |  ]  [  |    |  G |      G => ghost cell
%    +----+----+--]  [--+----+----+
%     TE(1)                   TE(end)      Cells (imax+1)
%      S(1)                    S(end)
%         ^                  ^
%        F(1)                F(end)        Faces (imax)


% for cell = 2:imax,
% 
%     TE(cell,:) = geom(x(cell)).*F(cell,:) - geom(x(cell-1)).*F(cell-1,:) ...
%         - S(cell,:).*dx(cell)';
%     
% end

for eq = 1:3,
    
    TE(2:parameter,eq) = geom(x(2:parameter)).*F(2:parameter,eq) - geom(x(1:parameter-1)).*F(1:parameter-1,eq); % ...
%         - S(2:imax,eq).*dx(2:imax)';
end

TE(2:parameter,2) = TE(2:parameter,2) - S.*dx(2:parameter)'; 

% Scale by cell volume:

% TE = TE./cell_vol;

save('flux.dat','F','-ascii','-double')
save('source.dat','S','-ascii','-double')

end

