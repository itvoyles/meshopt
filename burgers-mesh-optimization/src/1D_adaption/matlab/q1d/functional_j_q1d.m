function [ J ] = functional_j_q1d( x, TE )
%Functional "functional_j_q1d"
%   Takes current grid and TE, chooses type of functional
%   (cont, x-mtm, energy, or RSS), returns functional integrated over the
%   domain.

global functional TEmax
    TEmax1 = TEmax(1);
    TEmax2 = TEmax(2);
    TEmax3 = TEmax(3);
% Note: Removes ghost cell values [ TE(1,:) and TE(end,:) ]

if strcmp(functional, 'cont'),
    
%     J = sum( ( ( TE(2:end-1,1)/TEmax1).^2 ) .* ( x(2:end) - x(1:end-1)) );
    J = sqrt( sum( ( ( TE(2:end-1,1)/TEmax1).^2 ) .* (1) ) );
elseif strcmp(functional, 'xmtm'),
    
%     J = sum( ( ( TE(2:end-1,2)/TEmax2).^2 ) .* ( x(2:end) - x(1:end-1)) );
    J = sqrt( sum( ( ( TE(2:end-1,2)/TEmax2).^2 ) .* (1) ) );
elseif strcmp(functional, 'energy'),
    
%     J = sum( ( ( TE(2:end-1,3)/TEmax3).^2 ) .* ( x(2:end) - x(1:end-1)) );
J = sqrt( sum( ( ( TE(2:end-1,3)/TEmax3).^2 ) .* (1) ) );
    
elseif strcmp(functional, 'rss'),
    
    %     TEmax1 = max( abs(TE(:,1)) );
    %     TEmax2 = max( abs(TE(:,2)) );
    %     TEmax3 = max( abs(TE(:,3)) );
    
%     TEmax1 = TEmax(1);
%     TEmax2 = TEmax(2);
%     TEmax3 = TEmax(3);
%     
%     J = sqrt( ...
%         sum( ( (TE(2:end-1,1)/TEmax1).^2 ...
%         + (TE(2:end-1,2)/TEmax2).^2 ...
%         + (TE(2:end-1,3)/TEmax3).^2 ) ...
%         .* ( x(2:end) - x(1:end-1) )  ) );
    
    J = sqrt( ...
        sum( ( (TE(2:end-1,1)/TEmax1).^2 ...
        + (TE(2:end-1,2)/TEmax2).^2 ...
        + (TE(2:end-1,3)/TEmax3).^2 ) ...
        .* ( 1 )  ) );
else
    error( 'Enter a correct value for "functional"')
    
end
    


end

