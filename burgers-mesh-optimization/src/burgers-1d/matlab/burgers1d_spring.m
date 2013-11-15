function [J dJdx] = burgers1d_spring(k)
global nu alpha Lref forward mapping


x = springsystem(k);


dxi = 1;
% x_xi = dnfdxn(x,dxi,1,2);
% x_xixi = dnfdxn(x,dxi,2,2);
x_xi = dudx(x,dxi);
x_xixi = d2udx2(x,dxi);

TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref);

J = functional_j(x,TE);
% J = sum(abs(TE.^2));

% J = sqrt( sum(abs(TE).^2)/length(TE) );
% J = max(abs(TE));
end