function [ flux_vanleer ] = flux_vanleer( qL,qR )
%Function "flux_vanleer"
%   MATLAB version of Joe Derlaga's "flux_vanleer.f90",
%   created 070313.
%   Takes left and right state primitive variables [rho,u,p] and
%   returns Van Leer FVS flux.

global gm1 xg xgm1 xg2m1

% ----------------

%Calculate Left (+) Flux
  a = speed_of_sound(qL(3), qL(1));
  M = qL(2)/a;

  if ( abs(M) < 1 ), %Left subsonic flux
      fa = 0.25*qL(1)*a*(M+1)^2;
      fb = a*(gm1*M + 2);
      
      FL(1) = fa;
      FL(2) = fa*fb*xg;
      FL(3) = 0.5*fa*fb*fb*xg2m1;
      
  elseif ( M >= 1 ), %Left supersonic flux
      FL(1) = qL(1)*a*M;
      FL(2) = qL(1)*a^2*(M^2+xg);
      FL(3) = qL(1)*a^3*M*(0.5*M^2+xgm1);
  else
      FL = 0;
  end

%Calculate Right (-) Fluxes
  a = speed_of_sound(qR(3), qR(1));
  M = qR(2)/a;

  if ( abs(M) < 1 ), %Right subsonic flux
      fa = -0.25*qR(1)*a*(M-1)^2;
      fb = a*(gm1*M - 2);
      
      FR(1) = fa;
      FR(2) = fa*fb*xg;
      FR(3) = 0.5*fa*fb*fb*xg2m1;
      
  elseif ( M <= -1 ), %Right supersonic flux
      
      FR(1) = qR(1)*a*M;
      FR(2) = qR(1)*a^2*(M^2+xg);
      FR(3) = qR(1)*a^3*M*(0.5*M^2+xgm1);
      
  else
      FR = 0;
  end

  flux_vanleer = FL + FR;
% ----------------


end

