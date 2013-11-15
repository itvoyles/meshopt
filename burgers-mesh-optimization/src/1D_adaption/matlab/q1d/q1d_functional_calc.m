function [ J ] = q1d_functional_calc( DV )
%Function "q1d_functional"
%   Takes the design variable; calls the TE functional for the 
%   Q1D equations and returns it for optimization with fmincon.

global forward_mapping param cell_vol
global pb pzero sol order flux_type lim_type kappa shloc As2 flag

% Map design variable to grid nodes
x = forward_mapping(DV);

% Calculate truncation error
x = sort(x);
[ TE ] = q1d_TE( param,x,pb,pzero,sol,order,kappa,flag,flux_type,lim_type,shloc,As2);

% Calculate functional

J = functional_j_q1d(x,TE);

end

