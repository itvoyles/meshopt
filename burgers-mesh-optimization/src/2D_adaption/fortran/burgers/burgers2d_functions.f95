module burgers2d_functions

contains

elemental function uxn(x,y,a,v,l,flag)
  use select_precision, only : prec
  implicit none
  
  real(prec), intent(in) :: x,y,a,v,l
  integer,  intent(in)   :: flag
  real(prec) :: uxn
  real(prec) :: c1,c2,f,f1,f2,f3,f4,f5
  
 
      ! f = c1*tanh(c2x)
    c1 = -2._prec*v*a/L
    c2 = a/L

    f = c1*tanh(c2*(x+y))
    if (flag.eq.0) then
      uxn=f
      return
    endif

    f1 = -c2/c1*(f**2-c1**2)
    if (flag.eq.1) then
      uxn=f1
      return
    endif

    f2 = -c2/c1*(2._prec*f*f1)
    if (flag.eq.2) then
      uxn=f2
      return
    endif

    f3 = -2._prec*c2/c1*(f1*f1+f*f2)
    if (flag.eq.3) then
      uxn=f3
      return
    endif

    f4 = -2._prec*c2/c1*(3._prec*f1*f2+f*f3)
    if (flag.eq.4) then
      uxn=f4
      return
    endif
    
    f5 = -2._prec*c2/c1*(3._prec*(f2**2+f1*f3)+f1*f3+f*f4)
    if (flag.eq.5) then 
      uxn=f5
      return
    endif


    uxn=0._prec
    return
  
end function uxn

function Burgers_fd_cont_resid(x,y,a,v,l,imax,jmax)
  use select_precision, only : prec
  use derivatives
  use fd_derivatives
  use derivative_transform
  implicit none
  
  integer, intent(in) :: imax,jmax
  real(prec), intent(in) :: a,v,l
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  real(prec), dimension(imax,jmax):: Burgers_fd_cont_resid

  real(prec), dimension(imax,jmax)     :: u
  real(prec), dimension(imax,jmax) :: u_xi,u_xixi,u_eta,u_etaeta,u_xieta
  real(prec), dimension(imax,jmax) :: u_x, u_xx, u_y, u_yy
  real(prec), dimension(imax,jmax) :: x_xi,y_xi,x_eta,y_eta,Jinv
  real(prec), dimension(imax,jmax) :: xi_x,xi_y,eta_x,eta_y
  real(prec), dimension(imax,jmax) :: xi_xx, xi_yy, eta_xx, eta_yy
  real(prec), dimension(imax,jmax) :: xi_xxi, xi_xeta,eta_xxi,eta_xeta
  real(prec), dimension(imax,jmax) :: xi_yxi, xi_yeta,eta_yxi,eta_yeta
  real(prec), dimension(imax,jmax) :: term1, term2,temp
  real(prec) :: dxi, deta,Re
  integer :: order
  integer :: i,j


print*,'using wrong TE calc'
stop

  dxi = 1._prec !spacing in compuational space
  deta = 1._prec


  !metrics ******************
  do j = 1,jmax
  x_xi(:,j) = dudx(x(:,j),dxi,imax)
  y_xi(:,j) = dudx(y(:,j),dxi,imax)
  enddo
  
  do i = 1,imax
  x_eta(i,:) = dudx(x(i,:),deta,imax)
  y_eta(i,:) = dudx(y(i,:),deta,imax)
  enddo
  
  Jinv = x_xi*y_eta-x_eta*y_xi

  xi_x = y_eta/Jinv
  xi_y = -x_eta/Jinv
  eta_x = -y_xi/Jinv
  eta_y = x_xi/Jinv


  do j = 1,jmax-2
    xi_xxi(:,j) = dudx(xi_x(:,j),dxi,imax-2)
    eta_xxi(:,j) = dudx(eta_x(:,j),dxi,imax-2)
    xi_yxi(:,j) = dudx(xi_y(:,j),dxi,imax-2)
    eta_yxi(:,j) = dudx(eta_y(:,j),dxi,imax-2)
  enddo

  do j = 1,imax-2
    xi_xeta(j,:) = dudx(xi_x(j,:),deta,jmax-2)
    eta_xeta(j,:) = dudx(eta_x(j,:),deta,jmax-2)
    xi_yeta(j,:) = dudx(xi_y(j,:),deta,jmax-2)
    eta_yeta(j,:) = dudx(eta_y(j,:),deta,jmax-2)
  enddo

  xi_xx = xi_x*xi_xxi+eta_x*xi_xeta
  eta_xx = xi_x*eta_xxi+eta_x*eta_xeta
  xi_yy = xi_y*xi_yxi+eta_y*xi_yeta
  eta_yy = xi_y*eta_yxi+eta_y*eta_yeta


  !solution derivatives ****************
  u = uxn(x,y,a,v,L,0)

  u_xi(2:imax-1,2:imax-1) = fd_dudx(u(1:imax-2,2:jmax-1),u(3:imax,2:jmax-1),dxi)
  u_xixi(2:imax-1,2:imax-1) = fd_d2udx2(u(1:imax-2,2:jmax-1),u(2:imax-1,2:jmax-1),u(3:imax,2:jmax-1),dxi)
  u_eta(2:imax-1,2:imax-1) = fd_dudx(u(2:imax-1,1:jmax-2),u(2:imax-1,3:imax),deta)
  u_etaeta(2:imax-1,2:imax-1) = fd_d2udx2(u(2:imax-1,1:jmax-2),u(2:imax-1,2:jmax-1),u(2:imax-1,3:imax),deta)
  u_xieta(2:imax-1,2:imax-1) = fd_d2udxy(u(3:imax,3:jmax),u(3:imax,1:jmax-2),u(1:imax-2,3:jmax),u(1:imax-2,1:jmax-2),dxi,deta)


  u_x = dudx_from_dudxi(xi_x,eta_x,u_xi,u_eta)
  u_xx = d2udx2_from_d2udxi2(xi_x,eta_x,xi_xx,eta_xx,&
                             u_xi,u_eta,u_xixi,u_etaeta,u_xieta)


  u_y = dudy_from_dudeta(xi_y,eta_y,u_xi,u_eta)
  u_yy = d2udy2_from_d2udeta2(xi_y,eta_y,xi_yy,eta_yy,&
                              u_xi,u_eta,u_xixi,u_etaeta,u_xieta)


  Burgers_fd_cont_resid(2:imax-1,2:jmax-1) =  u(2:imax-1,2:jmax-1)&
                              *(u_x(2:imax-1,2:jmax-1)+u_y(2:imax-1,2:jmax-1))&
                            -v*(u_xx(2:imax-1,2:jmax-1)+u_yy(2:imax-1,2:jmax-1))


  return

  term1 = (xi_x+xi_y)*u*u_xi+(eta_y+eta_x)*u*u_eta

  term2 = (&
           ( xi_x**2+xi_y**2 )*u_xixi &
          -2._prec*( (xi_x*eta_x)+(xi_y*eta_y) )*u_xieta &
          +( eta_x**2+eta_y**2 )*u_etaeta &
          +(xi_xx-xi_yy)*u_xi &
          -(eta_xx-eta_yy)*u_eta &
          )
  

  !Burgers_fd_cont_resid = term1-v*term2




end function


function Truncation_error(x,y,a,v,l,imax,jmax)
  use select_precision, only : prec
  implicit none
  integer, intent(in) :: imax, jmax
  real(prec), intent(in) :: a,v,l
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  real(prec), dimension(imax,jmax) :: Truncation_error

  Real(kind=Prec),Dimension(imax,jmax) :: uloc   ! New Velocity


  Real(kind=Prec) ::  dudx            ! du/dx
  Real(kind=Prec) ::  dudy            ! du/dy
  Real(kind=Prec) ::  d2udx2          ! d2u/dx2
  Real(kind=Prec) ::  d2udy2          ! d2u/dy2
  Real(kind=Prec) ::  d2udx2CART          ! d2u/dx2
  Real(kind=Prec) ::  d2udy2CART          ! d2u/dy2

  Real(kind=Prec) ::  Laplacian         ! d2u/dx2 + d2u/dy2 (computed from transformed coordinates)
  Real(kind=Prec) ::  dudxsi            ! du/dxsi
  Real(kind=Prec) ::  dudeta            ! du/deta
  Real(kind=Prec) ::  d2udxsi2          ! d2u/dxsi2
  Real(kind=Prec) ::  d2udeta2          ! d2u/deta2
  Real(kind=Prec) ::  d2udxsideta       ! d2u/dxsi/deta (cross derivative term)
                                                    ! ***Note: Below, xsix is NOT xsix/J !!!***
  Real(kind=Prec) ::  term1             ! Temporary Variables
  Real(kind=Prec) ::  term2             ! Temporary Variables
  Real(kind=Prec) ::  term3             ! Temporary Variables


  Real(kind=Prec),Dimension(imax,jmax) ::  vol   ! volume of cell surrounding node
  Real(kind=Prec),Dimension(imax,jmax) ::  xsix  ! xsi_x/J grid metric
  Real(kind=Prec),Dimension(imax,jmax) ::  xsiy  ! xsi_y/J grid metric
  Real(kind=Prec),Dimension(imax,jmax) ::  etax  ! eta_x/J grid metric
  Real(kind=Prec),Dimension(imax,jmax) ::  etay  ! eta_y/J grid metric
  Real(kind=Prec),Dimension(imax,jmax) ::  d_xsix_dxsi  ! grid metric derivative
  Real(kind=Prec),Dimension(imax,jmax) ::  d_xsiy_dxsi  ! grid metric derivative
  Real(kind=Prec),Dimension(imax,jmax) ::  d_etax_dxsi  ! grid metric derivative
  Real(kind=Prec),Dimension(imax,jmax) ::  d_etay_dxsi  ! grid metric derivative
  Real(kind=Prec),Dimension(imax,jmax) ::  d_xsix_deta  ! grid metric derivative
  Real(kind=Prec),Dimension(imax,jmax) ::  d_xsiy_deta  ! grid metric derivative
  Real(kind=Prec),Dimension(imax,jmax) ::  d_etax_deta  ! grid metric derivative
  Real(kind=Prec),Dimension(imax,jmax) ::  d_etay_deta  ! grid metric derivative

  integer :: i,j

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

!  Calculate metrics and Jacobians
!       Notes: *metrics are actually xsix/J, etax/J, etc.
!              *vol(i,j) is the cell volume or inverse Jacobian
!     Metrics/Jacobians at Interior Points
  do i = 2, imax-1
  do j = 2, jmax-1
    xsix(i,j) =  0.5*(y(i,j+1) - y(i,j-1))
    xsiy(i,j) = -0.5*(x(i,j+1) - x(i,j-1))
    etax(i,j) = -0.5*(y(i+1,j) - y(i-1,j))
    etay(i,j) =  0.5*(x(i+1,j) - x(i-1,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
!    write(*,*) '**i,j: ',i,j
!    write(*,*) 'xsix,xsiy,etax,etay,vol: ',xsix(i,j),xsiy(i,j),
! &              etax(i,j),etay(i,j),vol(i,j)
  enddo
  enddo
      
!     Metrics/Jacobians at Boundaries
  do i = 2, imax-1
    j = 1
    xsix(i,j) =  0.5*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2))
    xsiy(i,j) = -0.5*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2))
    etax(i,j) = -0.5*(y(i+1,j) - y(i-1,j))
    etay(i,j) =  0.5*(x(i+1,j) - x(i-1,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
    j = jmax
    xsix(i,j) =  0.5*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2))
    xsiy(i,j) = -0.5*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2))
    etax(i,j) = -0.5*(y(i+1,j) - y(i-1,j))
    etay(i,j) =  0.5*(x(i+1,j) - x(i-1,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  enddo
  do j = 2, jmax-1    
    i = 1
    xsix(i,j) =  0.5*(y(i,j+1) - y(i,j-1))
    xsiy(i,j) = -0.5*(x(i,j+1) - x(i,j-1))
    etax(i,j) = -0.5*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j))
    etay(i,j) =  0.5*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
    i = imax
    xsix(i,j) =  0.5*(y(i,j+1) - y(i,j-1))
    xsiy(i,j) = -0.5*(x(i,j+1) - x(i,j-1))
    etax(i,j) = -0.5*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j))
    etay(i,j) =  0.5*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  enddo
! Metrics/Jacobians at Corners
  i = 1
  j = 1
  xsix(i,j) =  0.5*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2))
  xsiy(i,j) = -0.5*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2))
  etax(i,j) = -0.5*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j))
  etay(i,j) =  0.5*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  i = 1
  j = jmax
  xsix(i,j) =  0.5*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2))
  xsiy(i,j) = -0.5*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2))
  etax(i,j) = -0.5*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j))
  etay(i,j) =  0.5*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  i = imax
  j = 1
  xsix(i,j) =  0.5*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2))
  xsiy(i,j) = -0.5*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2))
  etax(i,j) = -0.5*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j))
  etay(i,j) =  0.5*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  i = imax
  j = jmax
  xsix(i,j) =  0.5*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2))
  xsiy(i,j) = -0.5*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2))
  etax(i,j) = -0.5*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j))
  etay(i,j) =  0.5*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)

  ! Metric derivatives
  do i = 2, imax-1
  do j = 2, jmax-1
    d_xsix_dxsi(i,j) = 0.5_Prec*( xsix(i+1,j)/vol(i+1,j) - xsix(i-1,j)/vol(i-1,j) )
    d_xsiy_dxsi(i,j) = 0.5_Prec*( xsiy(i+1,j)/vol(i+1,j) - xsiy(i-1,j)/vol(i-1,j) )
    d_etax_dxsi(i,j) = 0.5_Prec*( etax(i+1,j)/vol(i+1,j) - etax(i-1,j)/vol(i-1,j) )
    d_etay_dxsi(i,j) = 0.5_Prec*( etay(i+1,j)/vol(i+1,j) - etay(i-1,j)/vol(i-1,j) )
    d_xsix_deta(i,j) = 0.5_Prec*( xsix(i,j+1)/vol(i,j+1) - xsix(i,j-1)/vol(i,j-1) )
    d_xsiy_deta(i,j) = 0.5_Prec*( xsiy(i,j+1)/vol(i,j+1) - xsiy(i,j-1)/vol(i,j-1) )
    d_etax_deta(i,j) = 0.5_Prec*( etax(i,j+1)/vol(i,j+1) - etax(i,j-1)/vol(i,j-1) )
    d_etay_deta(i,j) = 0.5_Prec*( etay(i,j+1)/vol(i,j+1) - etay(i,j-1)/vol(i,j-1) )
  enddo
  enddo

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************


  Truncation_error = 0._prec
  
!  Save old u values at time level n
  uloc = uxn(x,y,a,v,l,0)

!  Apply the Euler explicit method for a Jacobi-like iteration

! Euler Explicit Method 
! Forward Sweep
  do j = 2, jmax - 1
  do i = 2, imax - 1
    dudxsi = 0.5_Prec*( uloc(i+1,j) - uloc(i-1,j) )
    dudeta = 0.5_Prec*( uloc(i,j+1) - uloc(i,j-1) )
    dudx = (xsix(i,j)*dudxsi + etax(i,j)*dudeta)/vol(i,j)
    dudy = (xsiy(i,j)*dudxsi + etay(i,j)*dudeta)/vol(i,j)

    d2udxsi2 =  uloc(i+1,j) - 2.0_Prec*uloc(i,j) + uloc(i-1,j) 
    d2udeta2 =  uloc(i,j+1) - 2.0_Prec*uloc(i,j) + uloc(i,j-1) 
    d2udxsideta = 0.5_Prec*( 0.5_Prec*( uloc(i+1,j+1) - uloc(i+1,j-1) )  &
                           - 0.5_Prec*( uloc(i-1,j+1) - uloc(i-1,j-1) ) )

    ! Terms for d2u/dx2
    term1 =   ( xsix(i,j)**2 )*d2udxsi2  &
            + 2.0_Prec*( etax(i,j)*xsix(i,j) )*d2udxsideta  &
            + ( etax(i,j)**2 )*d2udeta2
    term1 = term1/( vol(i,j)**2 )
    term2 = dudxsi*( xsix(i,j)*d_xsix_dxsi(i,j)  &
                   + etax(i,j)*d_xsix_deta(i,j) ) / vol(i,j)
    term3 = dudeta*( xsix(i,j)*d_etax_dxsi(i,j)  &
                   + etax(i,j)*d_etax_deta(i,j) ) / vol(i,j)
    d2udx2 = term1 + term2 + term3


    term1 =   ( xsiy(i,j)**2 )*d2udxsi2  &
            + 2.0_Prec*( etay(i,j)*xsiy(i,j) )*d2udxsideta  &
            + ( etay(i,j)**2 )*d2udeta2
    term1 = term1/( vol(i,j)**2 )
    term2 = dudxsi*( xsiy(i,j)*d_xsiy_dxsi(i,j)  &
                   + etay(i,j)*d_xsiy_deta(i,j) ) / vol(i,j)
    term3 = dudeta*( xsiy(i,j)*d_etay_dxsi(i,j)  &
                   + etay(i,j)*d_etay_deta(i,j) ) / vol(i,j)
    d2udy2 = term1 + term2 + term3

    Truncation_error(i,j) = uloc(i,j)*(dudx + dudy) - v*(d2udx2 + d2udy2)

  enddo
  enddo



End function



function burgers_fd_cont_resid_loc(x,y,a,v,l)
  use select_precision, only : prec
  use fd_derivatives
  use derivative_transform
  implicit none
  
  integer, parameter :: imax=5,jmax=5
  real(prec), intent(in) :: a,v,l
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  real(prec):: Burgers_fd_cont_resid_loc
  
  real(prec), dimension(imax,jmax) :: u
  real(prec) :: u_xi,u_xixi,u_eta,u_etaeta,u_xieta
  real(prec) :: u_x, u_xx, u_y, u_yy
  real(prec), dimension(imax,jmax) :: x_xi,y_xi,x_eta,y_eta,Jinv
  real(prec), dimension(imax,jmax) :: xi_x,xi_y,eta_x,eta_y
  real(prec), dimension(imax,jmax) :: xi_xx, xi_yy, eta_xx, eta_yy
  real(prec), dimension(imax,jmax) :: xi_xxi, xi_xeta,eta_xxi,eta_xeta
  real(prec), dimension(imax,jmax) :: xi_yxi, xi_yeta,eta_yxi,eta_yeta
  real(prec), dimension(imax,jmax) :: term1, term2,temp
  real(prec) :: dxi, deta,Re
  integer :: order
  integer :: i,j



  dxi = 1._prec !spacing in compuational space
  deta = 1._prec
  order = 2     !order of accuracy of derivative to calculate

x_xi = 0._prec
y_xi = 0._prec
x_eta = 0._prec
y_eta = 0._prec

!metrics ******************
x_xi(2,3) = fd_dudx(x(1,3),x(3,3),dxi)
x_xi(3,3) = fd_dudx(x(2,3),x(4,3),dxi)
x_xi(4,3) = fd_dudx(x(3,3),x(5,3),dxi)
x_xi(3,2) = fd_dudx(x(2,2),x(4,2),dxi)
x_xi(3,4) = fd_dudx(x(2,4),x(4,4),dxi)

y_xi(2,3) = fd_dudx(y(1,3),y(3,3),dxi)
y_xi(3,3) = fd_dudx(y(2,3),y(4,3),dxi)
y_xi(4,3) = fd_dudx(y(3,3),y(5,3),dxi)
y_xi(3,2) = fd_dudx(y(2,2),y(4,2),dxi)
y_xi(3,4) = fd_dudx(y(2,4),y(4,4),dxi)

x_eta(2,3) = fd_dudx(x(2,2),x(2,4),dxi)
x_eta(3,3) = fd_dudx(x(3,2),x(3,4),dxi)
x_eta(4,3) = fd_dudx(x(4,2),x(4,4),dxi)
x_eta(3,2) = fd_dudx(x(3,1),x(3,3),dxi)
x_eta(3,4) = fd_dudx(x(3,3),x(3,5),dxi)

y_eta(2,3) = fd_dudx(y(2,2),y(2,4),dxi)
y_eta(3,3) = fd_dudx(y(3,2),y(3,4),dxi)
y_eta(4,3) = fd_dudx(y(4,2),y(4,4),dxi)
y_eta(3,2) = fd_dudx(y(3,1),y(3,3),dxi)
y_eta(3,4) = fd_dudx(y(3,3),y(3,5),dxi)

Jinv = x_xi*y_eta-x_eta*y_xi

xi_x = y_eta/Jinv
xi_y = -x_eta/Jinv
eta_x = -y_xi/Jinv
eta_y = x_xi/Jinv





xi_xxi(3,3) = fd_dudx(xi_x(2,3),xi_x(4,3),dxi)
eta_xxi(3,3) = fd_dudx(eta_x(2,3),eta_x(4,3),dxi)
xi_yxi(3,3) = fd_dudx(xi_y(2,3),xi_y(4,3),dxi)
eta_yxi(3,3) = fd_dudx(eta_y(2,3),eta_y(4,3),dxi)

xi_xeta(3,3) = fd_dudx(xi_x(3,2),xi_x(3,4),deta)
eta_xeta(3,3) = fd_dudx(eta_x(3,2),eta_x(3,4),deta)
xi_yeta(3,3) = fd_dudx(xi_y(3,2),xi_y(3,4),deta)
eta_yeta(3,3) = fd_dudx(eta_y(3,2),eta_y(3,4),deta)

xi_xx = xi_x*xi_xxi+eta_x*xi_xeta
eta_xx = xi_x*eta_xxi+eta_x*eta_xeta
xi_yy = xi_y*xi_yxi+eta_y*xi_yeta
eta_yy = xi_y*eta_yxi+eta_y*eta_yeta




!solution derivatives ****************
  u = uxn(x,y,a,v,L,0)
  u_x = 0._prec
  u_xx = 0._prec
  u_y = 0._prec
  u_yy = 0._prec

  u_xi = fd_dudx(u(2,3),u(4,3),dxi)
  u_xixi = fd_d2udx2(u(2,3),u(3,3),u(4,3),dxi)
  u_eta = fd_dudx(u(3,2),u(3,4),deta)
  u_etaeta = fd_d2udx2(u(3,2),u(3,3),u(3,4),deta)
  u_xieta = fd_d2udxy(u(4,4),u(4,2),u(2,4),u(2,2),dxi,deta)


  u_x = dudx_from_dudxi(xi_x(3,3),eta_x(3,3),u_xi,u_eta)
  u_xx = d2udx2_from_d2udxi2(xi_x(3,3),eta_x(3,3),xi_xx(3,3),eta_xx(3,3),&
                             u_xi,u_eta,u_xixi,u_etaeta,u_xieta)


  u_y = dudy_from_dudeta(xi_y(3,3),eta_y(3,3),u_xi,u_eta)
  u_yy = d2udy2_from_d2udeta2(xi_y(3,3),eta_y(3,3),xi_yy(3,3),eta_yy(3,3),&
                              u_xi,u_eta,u_xixi,u_etaeta,u_xieta)


  Burgers_fd_cont_resid_loc =  u(3,3)*(u_x+u_y)-v*(u_xx+u_yy)






end function 


!test of function burgers_fd_cont_resid_loc
subroutine Burgers_fd_cont_resid_loc_test
  use select_precision, only : prec
  implicit none

  integer, parameter :: imax = 9,jmax = 9
  real(prec), dimension(imax,jmax) :: x,y,TE, xivec, etavec
  real(prec)  :: Tetest

  real(prec) :: a = 16._prec
  real(prec) :: uref = 2._prec
  real(prec) :: v
  real(prec) :: l = 8._prec
  real(prec) :: RE = 32._prec
  real(prec) :: pi

  integer :: i,j

  v = l*uref/Re
  pi = acos(-1._prec)

  do j = 1,jmax
  do i = 1,imax
  xivec(i,j) = real(i-1,prec)/real(imax-1,prec)
  enddo
  enddo

  do j = 1,jmax
  do i = 1,imax
  etavec(i,j) = real(j-1,prec)/real(jmax-1,prec)
  enddo
  enddo

  x = xivec**2-0.2_prec*sin(etavec*pi)
  y = etavec**2+0.5_prec*sin(xivec*pi)



  TE = Burgers_fd_cont_resid(x,y,a,v,l,imax,jmax)

  TEtest =  burgers_fd_cont_resid_loc(x(3:7,3:7),y(3:7,3:7),a,v,l)


  print*, TE(5,5), TE(5,5)-TEtest

end subroutine

!*******************************************************************************
!*******************************************************************************
!test of function Burgers_fd_cont_resid(x,y,a,v,l,imax,jmax)
!*******************************************************************************
!*******************************************************************************
subroutine Burgers_fd_cont_resid_test
  use select_precision, only : prec
  use derivative_transform
  implicit none

  integer, parameter :: imax=5, jmax=5
  real(prec), dimension(imax,jmax) :: u,u_x,u_y,u_xx,u_yy,u_xi,u_eta,u_xixi,u_etaeta, u_xieta
  real(prec),dimension(imax,jmax) :: x_xi,y_xi,x_eta,y_eta
  real(prec), dimension(imax,jmax) :: xi_x,xi_y,eta_x,eta_y,Jinv
  real(prec), dimension(imax,jmax) :: xi_xx, xi_yy, eta_xx, eta_yy
  real(prec), dimension(imax,jmax) :: xi_xxi, xi_xeta,eta_xxi,eta_xeta
  real(prec), dimension(imax,jmax) :: xi_yxi, xi_yeta,eta_yxi,eta_yeta,test
  real(prec) :: dxi,deta,dx
  integer    :: order = 2
  real(prec) :: a = 16._prec
  real(prec) :: uref = 2._prec
  real(prec) :: v
  real(prec) :: l = 8._prec
  real(prec) :: RE = 32._prec
  
  real(prec), dimension(2) :: dTE_dx_numerical, dTE_dx_analytic
  integer    :: n, i, j
  integer    :: nmax = 6
  real(prec) :: p
  real(prec) :: xi0,eta0,pi

  real(prec), dimension(imax,jmax) :: x,y,TE_exact,TE_test, xivec, etavec


  v = l*uref/Re
  pi = acos(-1._prec)

  write(*,*) '**************************************************'
  write(*,*) 'Test of the TE using the                          '
  write(*,*) 'continuous residual to calculate the exact TE.    '
  write(*,*) 'Expected: 2nd order accurate  '
  write(*,*) '**************************************************'
  write(*,*) ' I  |    dx      |    eps     |  p  '
  write(*,*) '--------------------------------------------------'

  100 format (I3,2e14.4,f7.2)
  xi0 = 0.5_prec
  eta0 = 0.5_prec
  dxi = 0.25_prec
  deta = 0.25_prec

  do i = 1,jmax
  xivec(:,i) = (/ 0._prec, 0.25_prec, 0.5_prec, 0.75_prec, 1._prec/)
  enddo
  do i = 1,imax
  etavec(i,:) = (/ 0._prec, 0.25_prec, 0.5_prec, 0.75_prec, 1._prec/)
  enddo



  do n = 1,nmax

  !shrink grid
  xivec = (xivec-xi0)/(2._prec)+xi0
  etavec = (etavec-eta0)/(2._prec)+eta0
  x = xivec**2-0.2_prec*sin(etavec*pi)
  y = etavec**2+0.5_prec*sin(xivec*pi)
  
  
  !write plot 3d file ******
  !write(n,*) 1
  !write(n,*) 5,5
  !write(*,'(4e23.15)')x
  !write(*,'(4e23.15)')y
  !*************************

  Te_test = Burgers_fd_cont_resid(x,y,a,v,l,imax,jmax)
  dTE_dx_numerical(1) = TE_test(3,3)

  !exact TE form continuous residual
  dTe_dx_analytic(1) = 0._prec
  
  !calculate observed order of accuracy
  if (n.gt.1) then

    p = log( (dTE_dx_analytic(2)-dTE_dx_numerical(2))/    &
             (dTE_dx_analytic(1)-dTE_dx_numerical(1)) )/  &
        log(2._prec)

    write(*,100) n, (x(4,3)-x(2,3))/2._prec, dTE_dx_numerical(1), p

  endif

  !transfer data
  dTE_dx_numerical(2) = dTE_dx_numerical(1)
  dTe_dx_analytic(2) = dTE_dx_analytic(1)

  enddo

  write(*,*) '--------------------------------------------------'
end subroutine


subroutine trap_sum2d_test
  use misc_func, only : trap_sum2d
  use select_precision, only : prec
  implicit none
  
  integer,parameter :: imax = 65
  integer,parameter :: jmax = 65
  integer,parameter :: imin = 5
  integer,parameter :: jmin = 5
  
  real(prec), dimension(2) :: dTE_dx_numerical, dTE_dx_analytic
  real(prec), dimension(imax,jmax) :: xi,eta,x,y,f
  real(prec) :: datatest, dataexact
  
  integer :: n,i,j,nmax,icur,jcur
  real(prec) :: pi, p
  
  
  write(*,*) '**************************************************'
  write(*,*) 'Test of trap_sum2d using an analytic integral on a'
  write(*,*) 'cartesian mesh.                                   '
  write(*,*) 'Expected: 2nd order accurate  '
  write(*,*) '**************************************************'
  write(*,*) ' I  |    dx      |    eps     |  p  '
  write(*,*) '--------------------------------------------------'
  
  
  
  
  100 format (I3,2e14.4,f7.2)
  
  !setup
  pi = acos(-1._prec)
  nmax = int(log(16._prec)/log(2._prec))
  
  icur = imin
  jcur = jmin
  
  !loop ************************
  do n = 1,nmax+1
  
  
  !create xi and eta coordinates
   xi = 0._prec
  eta = 0._prec
  do j = 1,jcur
    do i = 1,icur
       xi(i,j) = real(i-1,prec)/real(icur-1,prec)
      eta(i,j) = real(j-1,prec)/real(jcur-1,prec)
    enddo
  enddo
  
  !create physical grid
  !x = xi**2-0.1_prec*sin(eta*pi)
  !y = eta**2+0.1_prec*sin(xi*pi)
  !x = 0.6_prec*xi+0.4_prec*eta
  !y = 0.6_prec*eta+0.4_prec*xi
  x = xi
  y = eta
   
  !write plot 3d file ******
  !write(n,*) 1
  !write(n,*) icur,jcur
  !write(n,'(4e23.15)')x(1:icur,1:jcur)
  !write(n,'(4e23.15)')y(1:icur,1:jcur)
  !*************************

  f = x**3+y**4
  dataexact = y(icur,jcur)*x(icur,jcur)**4/4._prec&
             +x(icur,jcur)*y(icur,jcur)**5/5._prec

  !f = 1._prec
  !dataexact = x(icur,jcur)*y(icur,jcur)
  
  datatest=trap_sum2d(x(1:icur,1:jcur),y(1:icur,1:jcur),f(1:icur,1:jcur),icur,jcur)

  !print*, n,icur,jcur,dataexact,datatest

  dTE_dx_numerical(1) = datatest
  dTe_dx_analytic(1) = dataexact
  
  !calculate observed order of accuracy
  if (n.gt.1) then

    p = log( (dTE_dx_analytic(2)-dTE_dx_numerical(2))/    &
             (dTE_dx_analytic(1)-dTE_dx_numerical(1)) )/  &
        log(2._prec)

    write(*,100) n, (x(4,3)-x(2,3))/2._prec, &
                 dTE_dx_analytic(1)-dTE_dx_numerical(1), p

  endif




  !transfer data
  dTE_dx_numerical(2) = dTE_dx_numerical(1)
  dTe_dx_analytic(2) = dTE_dx_analytic(1)
  
  icur = (icur-1)*2+1
  jcur = (jcur-1)*2+1

  enddo
  write(*,*) '--------------------------------------------------'

end subroutine trap_sum2d_test


end module burgers2d_functions
