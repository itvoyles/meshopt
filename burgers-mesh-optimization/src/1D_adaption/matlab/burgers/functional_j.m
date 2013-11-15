% !  
% !      N-1  
% !      --
% !  J = \  f(i)
% !      /
% !      --
% !      i=1
% !
% !          1 /                     \/               \
% !  f(i) = ---| TE(i+1)^2 + TE(i)^2 || x(i+1) - x(i) |
% !          2 \                     /\               /
% !
% !  dJ     -1  /                       \ 
% ! ----- = --- | TE(i+1)^2 - TE(i-1)^2 | 
% ! dx(i)    2  \                       /
% !
% !          /         dTE(i+1)        dTE(i)  \/             \
% !        + | TE(i+1) --------- + TE(i)------ || x(i+1)-x(i) |
% !          \           dx              dx    /\             /
% !
% !          /        dTE(i)           dTE(i-1) \/             \
% !        + | TE(i) -------- + TE(i-1)-------- || x(i)-x(i-1) |
% !          \          dx                dx    /\             /
% !
% !          /          dTE(i+1)  \/               \
% !        + | TE(i+1) ---------- || x(i+2)-x(i+1) |
% !          \            dx      /\               /
% !
% !          /          dTE(i-1)  \/               \
% !        + | TE(i-1) ---------- || x(i-1)-x(i-2) |
% !          \            dx      /\               /


function J = functional_j(x,TE)
global Lref

% phi = burgers1d_penaltyfunction(x);
% q = 2;
% C = 1/(length(x)-1)*0.01;
% phihat = (C*phi).^q;
% 
% 
% 
% J = sum((...
%         TE(2:end).^2+TE(1:end-1).^2 ...
%        +(phihat(2:end)+phihat(1:end-1))./( x(2:end)-x(1:end-1) ).^q ...
%        ).*( x(2:end)-x(1:end-1) )/2 );
% 
% a1 = TE(2:end).^2+TE(1:end-1).^2;
% Dx = (x(3:end)-x(1:end-2))/(2);
% a2 = (phihat(2:end)+phihat(1:end-1))./( x(2:end)-x(1:end-1) ).^q;
% xc = (x(2:end)+x(1:end-1))/2;
% semilogy(xc,a1,'-ok',xc,a2,'-or')

% 
% q = 10;
% phithresh = 0.99;
% 
% phi = burgers1d_penaltyfunction(x);
% for i = 1:length(phi)
%    if (phi(i)<phithresh)
%        phi(i) = 0;
%    end    
% end
% % phi = (burgers1d_penaltyfunction(x)).^q;
% C1 = 1;
% 
% dxuni = Lref/(length(x)-1);
% dx = x(2:end)-x(1:end-1);
% 
% 
% phibar = (phi(2:end)+phi(1:end-1))/2;
% phihat = C1.*phibar.*(1-min(dxuni,dx)/dxuni);
% 
% J = sum((...
%         TE(2:end).^2+TE(1:end-1).^2 ...
%        +phihat./(x(2:end)-x(1:end-1)) ...
%        ).*( x(2:end)-x(1:end-1) )/2 );
% 
% a1 = TE(2:end).^2+TE(1:end-1).^2;
% xc = (x(2:end)+x(1:end-1))/2;
% 
% plot(xc,phihat,'-ro',xc,a1,'-ko')
% pause(0.001)

J = sum((...
        TE(2:end).^2+TE(1:end-1).^2 ...
       ).*( x(2:end)-x(1:end-1) )/2 );

end