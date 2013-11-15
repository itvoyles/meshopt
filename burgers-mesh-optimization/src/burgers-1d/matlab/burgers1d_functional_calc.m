function [J dJdx] = burgers1d_functional_calc(DV)
global nu alpha Lref forward_mapping method


% map design variable to grid nodes
x = forward_mapping(DV);

% calculate truncation error
dxi = 1;
x_xi = dudx(x,dxi);
x_xixi = d2udx2(x,dxi);
TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref);

% calculate functional
J = functional_j(x,TE);


% calculate gradient for 'pure' method 
if nargout > 1
    dJdx = calc_dJdx(x,nu,alpha,Lref);
    if ~strcmp(method,'pure')
        error('Gradient can only be calculated for "pure" method!')
    end
end


end