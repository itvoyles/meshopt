function [ p2c ] = primitive_to_conserved_1D( qp )
%Function "primitive_to_conserved_1D:
%   MATLAB version of Joe Derlaga's "primitive_to_conservative_1D.f90:,
%   created 070313.
%   Takes primitive variables [rho,u,p] and returns conserved variables.

global xgm1
p2c = zeros(3,1);
% --------------

p2c(1) = qp(1);
p2c(2) = qp(1) * qp(2);
p2c(3) = qp(3)*xgm1 + 0.5*qp(1)*qp(2)*qp(2);


end

