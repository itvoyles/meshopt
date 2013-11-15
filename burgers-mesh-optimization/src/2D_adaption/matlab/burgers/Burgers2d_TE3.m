function TE = Burgers2d_TE3(x,y,argin)
%% test
% RE = 32;
% Lref = 8;
% 
% %create initial mesh
% imax = 41;
% jmax = imax;
% 
% x = repmat(linspace(-Lref./2,Lref./2,imax)',1,jmax);
% y = repmat(linspace(-Lref./2,Lref./2,jmax),imax,1);

%% 
%  Calculate metrics and Jacobians
%       Notes: *metrics are actually xsix./J, etax./J, etc.
%              *vol(i,j) is the cell volume or inverse Jacobian
%     Metrics./Jacobians at Interior Points

RE = argin(1);
Lref = argin(2);


[imax,jmax] = size(x);

i = 2: imax-1;
j = 2: jmax-1;
    xsix(i,j) =  0.5.*(y(i,j+1) - y(i,j-1));
    xsiy(i,j) = -0.5.*(x(i,j+1) - x(i,j-1));
    etax(i,j) = -0.5.*(y(i+1,j) - y(i-1,j));
    etay(i,j) =  0.5.*(x(i+1,j) - x(i-1,j));
    vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);


      
%     Metrics./Jacobians at Boundaries
i = 2:imax-1;
j = 1;
    xsix(i,j) =  0.5.*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2));
    xsiy(i,j) = -0.5.*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2));
    etax(i,j) = -0.5.*(y(i+1,j) - y(i-1,j));
    etay(i,j) =  0.5.*(x(i+1,j) - x(i-1,j));
    vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);
j = jmax;
    xsix(i,j) =  0.5.*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2));
    xsiy(i,j) = -0.5.*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2));
    etax(i,j) = -0.5.*(y(i+1,j) - y(i-1,j));
    etay(i,j) =  0.5.*(x(i+1,j) - x(i-1,j));
    vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);

    
j = 2:jmax-1;    
i = 1;
    xsix(i,j) =  0.5.*(y(i,j+1) - y(i,j-1));
    xsiy(i,j) = -0.5.*(x(i,j+1) - x(i,j-1));
    etax(i,j) = -0.5.*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j));
    etay(i,j) =  0.5.*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j));
    vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);
i = imax;
    xsix(i,j) =  0.5.*(y(i,j+1) - y(i,j-1));
    xsiy(i,j) = -0.5.*(x(i,j+1) - x(i,j-1));
    etax(i,j) = -0.5.*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j));
    etay(i,j) =  0.5.*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j));
    vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);

% Metrics./Jacobians at Corners
i = 1;
j = 1;
  xsix(i,j) =  0.5.*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2));
  xsiy(i,j) = -0.5.*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2));
  etax(i,j) = -0.5.*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j));
  etay(i,j) =  0.5.*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j));
  vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);
  
i = 1;
j = jmax;
  xsix(i,j) =  0.5.*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2));
  xsiy(i,j) = -0.5.*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2));
  etax(i,j) = -0.5.*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j));
  etay(i,j) =  0.5.*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j));
  vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);
  
i = imax;
j = 1;
  xsix(i,j) =  0.5.*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2));
  xsiy(i,j) = -0.5.*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2));
  etax(i,j) = -0.5.*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j));
  etay(i,j) =  0.5.*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j));
  vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);
  
i = imax;
j = jmax;
  xsix(i,j) =  0.5.*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2));
  xsiy(i,j) = -0.5.*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2));
  etax(i,j) = -0.5.*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j));
  etay(i,j) =  0.5.*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j));
  vol(i,j) = etay(i,j).*xsix(i,j) - xsiy(i,j).*etax(i,j);

  % Metric derivatives
d_xsix_dxsi = zeros(imax,jmax);
d_xsiy_dxsi = d_xsix_dxsi;
d_etax_dxsi = d_xsix_dxsi;
d_etay_dxsi = d_xsix_dxsi;
d_xsix_deta = d_xsix_dxsi;
d_xsiy_deta = d_xsix_dxsi;
d_etax_deta = d_xsix_dxsi;
d_etay_deta = d_xsix_dxsi;
  
i = 2: imax-1;
j = 2: jmax-1;
    d_xsix_dxsi(i,j) = 0.5.*( xsix(i+1,j)./vol(i+1,j) - xsix(i-1,j)./vol(i-1,j) );
    d_xsiy_dxsi(i,j) = 0.5.*( xsiy(i+1,j)./vol(i+1,j) - xsiy(i-1,j)./vol(i-1,j) );
    d_etax_dxsi(i,j) = 0.5.*( etax(i+1,j)./vol(i+1,j) - etax(i-1,j)./vol(i-1,j) );
    d_etay_dxsi(i,j) = 0.5.*( etay(i+1,j)./vol(i+1,j) - etay(i-1,j)./vol(i-1,j) );
    d_xsix_deta(i,j) = 0.5.*( xsix(i,j+1)./vol(i,j+1) - xsix(i,j-1)./vol(i,j-1) );
    d_xsiy_deta(i,j) = 0.5.*( xsiy(i,j+1)./vol(i,j+1) - xsiy(i,j-1)./vol(i,j-1) );
    d_etax_deta(i,j) = 0.5.*( etax(i,j+1)./vol(i,j+1) - etax(i,j-1)./vol(i,j-1) );
    d_etay_deta(i,j) = 0.5.*( etay(i,j+1)./vol(i,j+1) - etay(i,j-1)./vol(i,j-1) );




%solution derivatives .^.^.^.^.^.^.^.^.^.^.^.^.^.^.^.^
  u = uxn(x,y,RE,Lref,0);
  nu = 2.*Lref./RE;

i=2:imax-1;
j=2:jmax-1;

    dudxsi = 0.5.*( u(i+1,j) - u(i-1,j) );
    dudeta = 0.5.*( u(i,j+1) - u(i,j-1) );
    dudx = (xsix(i,j).*dudxsi + etax(i,j).*dudeta)./vol(i,j);
    dudy = (xsiy(i,j).*dudxsi + etay(i,j).*dudeta)./vol(i,j);

    d2udxsi2 =  u(i+1,j) - 2.0.*u(i,j) + u(i-1,j) ;
    d2udeta2 =  u(i,j+1) - 2.0.*u(i,j) + u(i,j-1) ;
    d2udxsideta = 0.5.*( 0.5.*( u(i+1,j+1) - u(i+1,j-1) )  ...
                           - 0.5.*( u(i-1,j+1) - u(i-1,j-1) ) );

% Terms for d2u/dx2

    term1 =   ( xsix(i,j).^2 ).*d2udxsi2  ...
            + 2.0.*( etax(i,j).*xsix(i,j) ).*d2udxsideta  ...
            + ( etax(i,j).^2 ).*d2udeta2;
    term1 = term1./( vol(i,j).^2 );
    term2 = dudxsi.*( xsix(i,j).*d_xsix_dxsi(i,j)  ...
                   + etax(i,j).*d_xsix_deta(i,j) ) ./ vol(i,j);
    term3 = dudeta.*( xsix(i,j).*d_etax_dxsi(i,j)  ...
                   + etax(i,j).*d_etax_deta(i,j) ) ./ vol(i,j);
    d2udx2 = term1 + term2 + term3;

    term1 =   ( xsiy(i,j).^2 ).*d2udxsi2  ...
            + 2.0.*( etay(i,j).*xsiy(i,j) ).*d2udxsideta  ...
            + ( etay(i,j).^2 ).*d2udeta2;
    term1 = term1./( vol(i,j).^2 );
    term2 = dudxsi.*( xsiy(i,j).*d_xsiy_dxsi(i,j)  ...
                   + etay(i,j).*d_xsiy_deta(i,j) ) ./ vol(i,j);
    term3 = dudeta.*( xsiy(i,j).*d_etay_dxsi(i,j)  ...
                   + etay(i,j).*d_etay_deta(i,j) ) ./ vol(i,j);
    d2udy2 = term1 + term2 + term3;

    TE = zeros(imax,jmax);
    TE(2:end-1,2:end-1) = (u(2:end-1,2:end-1).*dudx + u(2:end-1,2:end-1).*dudy - nu.*(d2udx2 + d2udy2));
    
    
  TE(1,:)=0;
  TE(end,:)=0;
  TE(:,1)=0;
  TE(:,end)=0;
  

end