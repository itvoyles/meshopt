function [ source_out ] = construct_source( prim_cc,dadx_cc )
%Function "construct_source"
%   Adapted for MATLAB from Joe Derlaga's "residual.f90", 071113.
%   Constructs the Q1D area source term for the x-momentum equation.

global param

source = prim_cc(:,3) .* dadx_cc(:);

% Do not return values over ghost cells:

source_out = source(2:param);

end

