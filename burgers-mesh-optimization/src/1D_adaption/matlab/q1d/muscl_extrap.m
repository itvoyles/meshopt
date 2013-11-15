function [ qL,qR ] = muscl_extrap( vars_cc,kappa,lim_type )
%Function "muscl_extrap"
%   Adapted for MATLAB from Joe Derlaga's "residual.f90", 071113.
%   Takes cell-centered primitive variables, limiter type, and
%   MUSCL parameter; returns left and right states at faces.

global param

epsil = 0.0000001; % small value
qL = zeros(param,3);
qR = zeros(param,3);

if ~strcmp(lim_type,'off'),
    
    % Left variations
    r_L = zeros(param+2,3);
    
    % for i = 1:imax-1,
    %       r_L(i,:) = max( 0.0, (vars_cc(i+2,:) - vars_cc(i+1,:)) ...
    %           ./ (vars_cc(i+1,:) - vars_cc(i,:) + epsil));
    % end
    
        r_L(1:param-1,:) = max( 0.0, (vars_cc(3:param+1,:) - vars_cc(2:param,:)) ...
            ./ (vars_cc(2:param,:) - vars_cc(1:param-1,:) + epsil));

    
    % Right variations
    r_R = zeros(param+2,3);
    
    
    % for i = 2:imax,
    %      r_R(i,:) = max( 0.0, (vars_cc(i,:) - vars_cc(i-1,:)) ...
    %          ./ (vars_cc(i+1,:) - vars_cc(i,:) + epsil) );
    % end
    
        r_R(2:param,:) = max( 0.0, (vars_cc(2:param,:) - vars_cc(1:param-1,:)) ...
            ./ (vars_cc(3:param+1,:) - vars_cc(2:param,:) + epsil) );
    

    %~~~ Debugging
    save test.mat vars_cc r_L r_R
    %~~~ Debugging
    
    psi_L = zeros(param,3);
    psi_R = zeros(param,3);

%     switch lim_type
%         case 'minmod'
%             for i = 1:imax,
%                 psi_L(i,:) = limiter_sweby(3, r_L(i+1,:),1.0);
%                 psi_R(i,:) = limiter_sweby(3, r_R(i+1,:),1.0);
%             end
%         case 'superbee'
%             for i = 1:imax,
%                 psi_L(i,:) = limiter_sweby(3, r_L(i+1,:),2.0);
%                 psi_R(i,:) = limiter_sweby(3, r_R(i+1,:),2.0);
%             end
%         case 'sweby'
%             for i = 1:imax,
%                 psi_L(i,:) = limiter_sweby(3, r_L(i+1,:),1.5);
%                 psi_R(i,:) = limiter_sweby(3, r_R(i+1,:),1.5);
%             end
%         case 'ospre'
%             for i = 1:imax,
%                 psi_L(i,:) = limiter_ospre(3, r_L(i+1,:));
%                 psi_R(i,:) = limiter_ospre(3, r_R(i+1,:));
%             end
%         case 'vanleer'
%             for i = 1:imax,
%                 psi_L(i,:) = limiter_vanleer(3, r_L(i+1,:));
%                 psi_R(i,:) = limiter_vanleer(3, r_R(i+1,:));
%             end
%         case 'vanalbada'
%             for i = 1:imax,
%                 psi_L(i,:) = limiter_vanalbada(3, r_L(i+1,:));
%                 psi_R(i,:) = limiter_vanalbada(3, r_R(i+1,:));
%             end
%     end
%     
        switch lim_type
        case 'minmod'
            for i = 1:param,
                psi_L(i,:) = limiter_sweby(3, r_L(i,:),1.0);
                psi_R(i,:) = limiter_sweby(3, r_R(i,:),1.0);
            end
        case 'superbee'
            for i = 1:param,
                psi_L(i,:) = limiter_sweby(3, r_L(i,:),2.0);
                psi_R(i,:) = limiter_sweby(3, r_R(i,:),2.0);
            end
        case 'sweby'
            for i = 1:param,
                psi_L(i,:) = limiter_sweby(3, r_L(i,:),1.5);
                psi_R(i,:) = limiter_sweby(3, r_R(i,:),1.5);
            end
        case 'ospre'
            for i = 1:param,
                psi_L(i,:) = limiter_ospre(3, r_L(i,:));
                psi_R(i,:) = limiter_ospre(3, r_R(i,:));
            end
        case 'vanleer'
            for i = 1:param,
                psi_L(i,:) = limiter_vanleer(3, r_L(i,:));
                psi_R(i,:) = limiter_vanleer(3, r_R(i,:));
            end
        case 'vanalbada'
            for i = 1:param,
                psi_L(i,:) = limiter_vanalbada(3, r_L(i,:));
                psi_R(i,:) = limiter_vanalbada(3, r_R(i,:));
            end
        end
    
elseif strcmp(lim_type,'off'),
    psi_L = zeros(param,3);
    psi_R = zeros(param,3);
    psi_L(:,:) = 1.0;
    psi_R(:,:) = 1.0;
end
    


% MUSCL extrapolation

% Inflow face: "zero gradient extrapolation to ghost cell"

% qL(1,:) = (1/3)*( 4*vars_cc(2,:) - vars_cc(3,:) );
qL(1,:) = vars_cc(1,:);
% for i = 2:imax,
%     
%     qL(i,:) = vars_cc(i,:) + 0.25.*( (1+kappa).*psi_R(i,:) ...
%         .* (vars_cc(i+1,:) - vars_cc(i,:)) ...
%         + (1-kappa).*psi_L(i-1,:) .* (vars_cc(i,:) - vars_cc(i-1,:)));
%     
% end

qL(2:param,:) = vars_cc(2:param,:) + 0.25.*( (1+kappa).*psi_R(2:param,:) ...
    .* (vars_cc(3:param+1,:) - vars_cc(2:param,:)) ...
    + (1-kappa).*psi_L(1:param-1,:) .* (vars_cc(2:param,:) - vars_cc(1:param-1,:)));
        
% for i = 1:imax-1,
%     
%     qR(i,:) = vars_cc(i+1,:) - 0.25.*( (1-kappa).*psi_R(i+1,:) ...
%         .* (vars_cc(i+2,:) - vars_cc(i+1,:)) ...
%         + (1+kappa).*psi_L(i,:) .* (vars_cc(i+1,:) - vars_cc(i,:)));
%     
% end

qR(1:param-1,:) = vars_cc(2:param,:) - 0.25.*( (1-kappa).*psi_R(2:param,:) ...
    .* (vars_cc(3:param+1,:) - vars_cc(2:param,:)) ...
    + (1+kappa).*psi_L(1:param-1,:) .* (vars_cc(2:param,:) - vars_cc(1:param-1,:)));

% Outflow face: "zero gradient extrapolation to ghost"

% qR(imax,:) = (1/3)*(4*vars_cc(imax,:) - vars_cc(imax-1,:) );
qR(param,:) = vars_cc(end,:);

end

