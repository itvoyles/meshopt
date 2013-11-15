function J = calc_J(x,y,RE,Lref,ex)
global TEfunc argin

    

% TE = Burgers2d_TE3(x,y,RE,Lref);
TE = TEfunc(x,y,argin);
% TE1storder_terms = Burgers2d_TE1storder_terms(x,y,RE,Lref);

if (nargin <5)
    ex = [0,0,0,0];
    J = trap_sum2d(x,y,TE.^2);
    
else
        
    J = trap_sum2d(x(1+ex(1):end-ex(2),1+ex(3):end-ex(4)),y(1+ex(1):end-ex(2),1+ex(3):end-ex(4)),...
        TE(1+ex(1):end-ex(2),1+ex(3):end-ex(4) ).^2 );
end


% J = trap_sum2d(x,y,TE.^2)+trap_sum2d(x,y,TE1storder_terms.^2);
% J = trap_sum2d(x,y,(TE-TE1storder_terms).^2)+trap_sum2d(x,y,TE1storder_terms.^2);
% J = trap_sum2d(x,y,TE1storder_terms.^2);

% J = sum(sum(abs(TE).^2));

end