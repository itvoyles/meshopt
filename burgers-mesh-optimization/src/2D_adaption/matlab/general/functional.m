function [varargout] = functional(DV)
global RE Lref forward_mapping func_method dxdfx dydfy method imax jmax graditer iter grad 


[x,y] = forward_mapping(DV);

J = calc_j(x,y,RE,Lref);

% if func_method = 2, then truncation error and function are evaluated on a
% mesh refined by a factor of 2. For a formal order of accuracy of 2, TE
% should decrease by a factor of 4, and the functional is the integral of
% TE^2 making the expected difference between the two functional values
% (2^2)^2 = 16
if func_method == 2
    x2 = x(1:2:end,1:2:end);
    y2 = y(1:2:end,1:2:end);
    J2 = calc_j(x2,y2,RE,Lref);
    
    J = 16*J+J2;
end

if nargout>1
    
% if iter~=graditer &&  strcmp('coupled-spring-f',method)
%     graditer = iter;
    
        djdxmat = calc_djdx(x,y);
        
        if strcmp('coupled-spring-f',method)
        grad = [(reshape(djdxmat(:,:,1),1,imax*jmax)*dxdfx)'
                (reshape(djdxmat(:,:,2),1,imax*jmax)*dydfy)'];
            
        elseif strcmp('pure',method)
        grad= [reshape(djdxmat(:,:,1),imax*jmax,1)
               reshape(djdxmat(:,:,2),imax*jmax,1)];
        
        else
            error('No user defined gradient for given method.')
        end
            

% end

end
    


varargout{1} = J;
varargout{2} = grad;

end