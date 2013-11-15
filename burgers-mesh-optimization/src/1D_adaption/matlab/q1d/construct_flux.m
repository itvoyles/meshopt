function [ flux ] = construct_flux( prim_cc,kappa,flux_type,lim_type )
%Function "construct_flux"
%   Adapted for MATLAB from Joe Derlaga's "residual.f90", 071113.
%   Selects and implements chosen flux scheme for the Q1D Euler
%   Euler equations. Calls the MUSCL extrapolation.

global param

[prim_L,prim_R] = muscl_extrap( prim_cc,kappa,lim_type );

flux = zeros(param,3);

switch flux_type
    
    case 'roe'
        for i = 1:param,
            flux(i,:) = flux_roe(prim_L(i,:),prim_R(i,:));
        end
    case 'vanleer'
        for i = 1:param,
            flux(i,:) = flux_vanleer(prim_L(i,:),prim_R(i,:));
        end
        
end


end

