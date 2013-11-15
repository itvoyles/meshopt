% !****************************************************************************.^*
% !****************************************************************************.^*
% !                    function TE, dTE/dx, and test routines
% !****************************************************************************.^*
% !****************************************************************************.^*
% !verified accurate by comparing to second-order accurate finite difference
function dTEdx = dTE_dxi(u,ux,uxx,uxxx,uxxxx,uxxxxx,Xxi,Xxixi,...
                           dxi,nu,a,L,flag)
% ! dTE_dxi = pdTE_pdxi + dTE/dAi * dAi/dxi


%   !calculate metric derivatives based on flag
%   !solution derivatives are calculated outside of subroutine
  switch flag
      case {0} %!dTE(x_i)/dxi 
%     !solution and derivatives evaluated at xi
%     ! verified accurate with finite difference

%       !calculate metric derivatives
      dXxi_dx = 0./(2.*dxi);
      dXxixi_dx = -2./dxi.^2;

%       !calculate coefficients and derivatives
      A1 = (...
            1./6.*u.*uxxx     ...
           -1./12.*nu.*uxxxx  ...
            );

      dA1_dx = (1./6.*(ux.*uxxx+u.*uxxxx) ...
               -1./12.*nu.*uxxxxx);

      A2 = (...
           -1./3.*nu.*uxxx                    ...
           +1./2.*u.*uxx                      ...
           );

      dA2_dx = (1./2.*(ux.*uxx+u.*uxxx) ...
               -1./3.*nu.*uxxxx);

      A3 = nu.*uxx./4;

      dA3_dx = (1./4.*nu.*uxxx);

%       !calculate the partial of TE with respect to x
      pdTE_pdx = A1.*2.*Xxi.*dXxi_dx                                        ...
                +A2.*dXxixi_dx                                                  ...
                +A3.*2.*(Xxixi./Xxi).*(dXxixi_dx./Xxi-Xxixi./Xxi.^2.*dXxi_dx);

      dTE_dA1 = Xxi.^2;

      dTE_dA2 = Xxixi;

      dTE_dA3 = (Xxixi./Xxi).^2;

      dTEdx = dxi.^2.*(pdTE_pdx                  ...
                       +dTE_dA1.*dA1_dx           ...
                       +dTE_dA2.*dA2_dx            ...
                       +dTE_dA3.*dA3_dx);


      case {1} %!dTE(x_i+1)./dxi ****************************************************
%     ! solution and derivatives evaluated at xip1
%     !verified accurate with finite difference
    
%       !calculate metric derivatives
      dXxi_dx = -1./(2.*dxi);
      dXxixi_dx = 1./dxi.^2;

%       !calculate coefficients and derivatives
      A1 = (...
            1./6.*u.*uxxx     ...
           -1./12.*nu.*uxxxx  ...
            );

      dA1_dx = 0 ;

      A2 = (...
           -1./3.*nu.*uxxx                    ...
           +1./2.*u.*uxx                      ...
           );

      dA2_dx = 0 ;

      A3 = nu.*uxx./4;

      dA3_dx = 0;%!(1./4.*nu.*uxxx)

%       !calculate the partial of TE with respect to x
      pdTE_pdx = A1.*2.*Xxi.*dXxi_dx                                        ...
                +A2.*dXxixi_dx                                                  ...
                +A3.*2.*(Xxixi./Xxi).*(dXxixi_dx./Xxi-Xxixi./Xxi.^2.*dXxi_dx);

      dTE_dA1 = Xxi.^2;

      dTE_dA2 = Xxixi;

      dTE_dA3 = (Xxixi./Xxi).^2;

      dTEdx = dxi.^2.*(pdTE_pdx                  ...
                       +dTE_dA1.*dA1_dx           ...
                       +dTE_dA2.*dA2_dx            ...
                       +dTE_dA3.*dA3_dx);
    
    
      case {-1}% !dTE(x_i+1)./dxi ****************************************************
%     ! solution and derivatives evaluated at xip1
%     !verified accurate with finite difference
    
%       !calculate metric derivatives
      dXxi_dx = 1./(2.*dxi);
      dXxixi_dx = 1./dxi.^2;

%       !calculate coefficients and derivatives
      A1 = (...
            1./6.*u.*uxxx     ...
           -1./12.*nu.*uxxxx  ...
            );

      dA1_dx = 0 ;

      A2 = (...
           -1./3.*nu.*uxxx                    ...
           +1./2.*u.*uxx                      ...
           );

      dA2_dx = 0 ;

      A3 = nu.*uxx./4;

      dA3_dx = 0;%!(1./4.*nu.*uxxx)

%       !calculate the partial of TE with respect to x
      pdTE_pdx = A1.*2.*Xxi.*dXxi_dx                                        ...
                +A2.*dXxixi_dx                                                  ...
                +A3.*2.*(Xxixi./Xxi).*(dXxixi_dx./Xxi-Xxixi./Xxi.^2.*dXxi_dx);

      dTE_dA1 = Xxi.^2;

      dTE_dA2 = Xxixi;

      dTE_dA3 = (Xxixi./Xxi).^2;

      dTEdx = dxi.^2.*(pdTE_pdx                  ...
                       +dTE_dA1.*dA1_dx           ...
                       +dTE_dA2.*dA2_dx            ...
                       +dTE_dA3.*dA3_dx);
                                              
  end
  



end