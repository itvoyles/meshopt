function [ flux_roe ] = flux_roe( qL,qR )
%Function "flux_roe"
%   MATLAB version of Joe Derlaga's "flux_roe.f90",
%   created 070313.
%   Takes a left and right state primitive variables [rho,u,p],
%   returns Roe FDS flux.

global gm1

% ----------------

  lambdaeps = 0.1;
  lambdaRoe = zeros(3,1);
  
  rhoL = qL(1);
  uL   = qL(2);
  pL   = qL(3);

  rhoR = qR(1);
  uR   = qR(2);
  pR   = qR(3);

  consL = primitive_to_conserved_1D(qL);
  consR = primitive_to_conserved_1D(qR);

  FL = [rhoL*uL, rhoL*uL^2 + pL, uL*(consL(3) + pL)];
  FR = [rhoR*uR, rhoR*uR^2 + pR, uR*(consR(3) + pR)];

% Roe interface variable
  Rhalf = sqrt(rhoR/rhoL);

% Calculate Roe average variables
  RoeAvgrho = Rhalf*rhoL;
  RoeAvgu   = (Rhalf*uR + uL) / (Rhalf + 1);
  RoeAvght  = (Rhalf*(consR(3) + pR)/rhoR + (consL(3) + pL)/rhoL) / (Rhalf + 1);
  RoeAvga   = sqrt(gm1*(RoeAvght-0.5*RoeAvgu^2));

  lambdaRoe(1) = RoeAvgu;
  lambdaRoe(2) = RoeAvgu + RoeAvga;
  lambdaRoe(3) = RoeAvgu - RoeAvga;

%Entropy fix
  for i = 1:3,
    if ( abs(lambdaRoe(i)) <= 2*lambdaeps*RoeAvga ),
      lambdaRoe(i) = (lambdaRoe(i)^2)/(4*lambdaeps*RoeAvga) ...
                   + lambdaeps*RoeAvga;
    end
  end

  dw(1) = (rhoR-rhoL) - (pR-pL)/RoeAvga^2;
  dw(2) = 0.5*( (pR-pL) + RoeAvgrho*RoeAvga*(uR-uL) ) / (RoeAvga*RoeAvga);
  dw(3) = 0.5*( (pR-pL) - RoeAvgrho*RoeAvga*(uR-uL) ) / (RoeAvga*RoeAvga);

  r1RoeAvg = [1, RoeAvgu, 0.5*RoeAvgu^2];
  r2RoeAvg = [1, RoeAvgu+RoeAvga, RoeAvght+RoeAvgu*RoeAvga];
  r3RoeAvg = [1, RoeAvgu-RoeAvga, RoeAvght-RoeAvgu*RoeAvga];

%Calculate Interface Fluxes
  F = 0.5*( (FL + FR)                       ...
    - (abs(lambdaRoe(1))*dw(1)*r1RoeAvg      ...
    +  abs(lambdaRoe(2))*dw(2)*r2RoeAvg      ...
    +  abs(lambdaRoe(3))*dw(3)*r3RoeAvg) );

  flux_roe = F;
% ----------------


end

