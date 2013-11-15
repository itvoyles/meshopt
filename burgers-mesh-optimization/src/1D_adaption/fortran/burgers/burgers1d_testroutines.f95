module burgers1d_testroutines
use burgers1d_functions
!contains subroutines to test
! uxn(x,a,v,l,derivative)::
!     analytic derivatives of burgers equation up to fifth derivative

! Burgers_TE(x,x_xi,x_xixi,dxi,v,a,L,f)
!     analytic truncation error for finite difference expression

! dTE_dxi(u,ux,uxx,uxxx,uxxxx,uxxxxx,xip1,xi,xim1,dxi,nu,a,L,flag)
!     analytic derivative of TE, for flag=0. dTE(xi)/dxi

! dTE_dxi(u,ux,uxx,uxxx,uxxxx,uxxxxx,xip1,xi,xim1,dxi,nu,a,L,flag)
!     analytic derivative of TE, for flag=1. dTE(xi)/dxim1


contains

subroutine run_burgers1d_tests
  implicit none
  

  call burgers_derivative_test
  call Burgers_TE_test
  call dTE_dxi_test
  call dTEip1_dxi_test
  call dTEim1_dxi_test  

end subroutine run_burgers1d_tests


!*******************************************************************************
!*******************************************************************************
!test of function uxn
!*******************************************************************************
!*******************************************************************************
subroutine burgers_derivative_test
  use select_precision
  implicit none
  
  real(prec) :: x,dx,xip1,xim1,f,fip1,fim1
  real(prec), dimension(2) :: fanalytic,fnumerical
  integer :: nmax = 5
  integer :: i,n,j
  real(prec) :: p
  
  integer :: curder = 4
  
  real(prec) :: L  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: a = 16._prec
  real(prec) :: Re    = 32._prec
  real(prec) :: v
  
  v = L*uref/Re
  
  write(*,*) '**************************************************'
  write(*,*) 'Test of the analytic derivatives of burgers       '
  write(*,*) 'equation using finite difference calculations.    '
  write(*,*) 'Expected: 2nd order accurate                      '
  write(*,*) '**************************************************'
  
    
    
  100 format (I3,e14.4,f7.2)
  

  
  
  do curder = 1,5
  
  write(*,*) ' Derivative test: ',curder
  write(*,*) '--------------------------------------------------'
  write(*,*) ' I  |    dx      |  p  '
  write(*,*) '--------------------------------------------------'
  
  x = 0.5_prec
  dx = 0.1_prec
  
  do n = 1,nmax
  
  xip1 = x+dx
  xim1 = x-dx
  
  fip1 = uxn(xip1,v,a,l,curder-1) 
  fim1 = uxn(xim1,v,a,l,curder-1)
  
  fnumerical(1) = (fip1-fim1)/(2._prec*dx)
  
  fanalytic(1) = uxn(x,v,a,l,curder)
  
   
  
  if (n.gt.1) then
  
    p = log( (fnumerical(2)-fanalytic(2))/(fnumerical(1)-fanalytic(1)) )&
       /log(2._prec)
      
    write(*,100) n,dx,p
  
  endif
  
  fnumerical(2) = fnumerical(1)
  fanalytic(2) = fanalytic(1)
  
  dx = dx/2._prec 
  
  
  enddo
  write(*,*) '**************************************************'
  enddo

end subroutine burgers_derivative_test


!*******************************************************************************
!*******************************************************************************
!test of function Burgers_TE(x,x_xi,x_xixi,dxi,v,a,L,f)
!*******************************************************************************
!*******************************************************************************
subroutine Burgers_TE_test
  use select_precision, only : prec
  use fd_derivative_calc
  implicit none

  real(prec) :: dx,xip1,xi,xim1,xip2,xim2
  real(prec) :: TEip1, TEi, TEim1
  real(prec), dimension(2) :: dTE_dx_numerical, dTE_dx_analytic
  integer    :: n
  integer    :: imax = 6
  real(prec) :: p
  real(prec) :: dxi,dxi0,x0,xi0

  real(prec), dimension(5) :: x,TE_exact,TE_test,x_xi,x_xixi,xivec
  real(prec) :: Lref  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: alpha = 16._prec
  real(prec) :: Re    = 32._prec
  real(prec) :: nu
  
  nu = Lref*uref/Re

  write(*,*) '**************************************************'
  write(*,*) 'Test of the analytic derivation of TE using the   '
  write(*,*) 'continuous residual to calculate the exact TE.    '
  write(*,*) 'Expected: 4th order accurate  '
  write(*,*) '**************************************************'
  write(*,*) ' I  |    dx      |  p  '
  write(*,*) '--------------------------------------------------'

  100 format (I3,e14.4,f7.2)
  x0 = 0.5_prec
  xi0 = 0.5_prec
  dxi0 = 0.25_prec
  xivec = (/ 0._prec, 0.25_prec, 0.5_prec, 0.75_prec, 1._prec/)
  dxi = dxi0

  do n = 1,imax

  !shrink grid
  xivec = (xivec-xi0)/(2._prec)+xi0
  x = xivec**4
  x0 = xi0**4
  dxi = dxi/2._prec

  !calculate grid metrics
  call dnfdxn(x_xi,x,dxi,1,2,5)
  call dnfdxn(x_xixi,x,dxi,2,2,5)

  !derived analytic expression for TE
  Te_test = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref)
  dTE_dx_numerical(1) = TE_test(3)

  !exact TE form continuous residual
  Te_exact = Burgers_fd_cont_resid(x,5,nu,alpha,Lref)
  dTE_dx_analytic(1) = Te_exact(3)

  !calculate observed order of accuracy
  if (n.gt.1) then

    p = log( (dTE_dx_analytic(2)-dTE_dx_numerical(2))/    &
             (dTE_dx_analytic(1)-dTE_dx_numerical(1)) )/  &
        log(2._prec)

    write(*,100) n, (x(4)-x(2))/2._prec, p

  endif

  !transfer data
  dTE_dx_numerical(2) = dTE_dx_numerical(1)
  dTe_dx_analytic(2) = dTE_dx_analytic(1)

  enddo

  write(*,*) '--------------------------------------------------'
end subroutine

!*******************************************************************************
!*******************************************************************************
!test for dTE_dxi(u,ux,uxx,uxxx,uxxxx,uxxxxx,xip1,xi,xim1,dxi,nu,a,L,flag)
! flag option 0: derivative of TE(xi) wrt xi
!*******************************************************************************
!*******************************************************************************
subroutine dTE_dxi_test
  use select_precision, only : prec
  implicit none

  real(prec) :: dx,xip1,xi,xim1,xip2,xim2,x,xxi
  real(prec) :: TEip1, TEi, TEim1
  real(prec), dimension(2) :: dTE_dxi_numerical, dTE_dxi_analytic,error
  integer    :: n
  integer    :: imax = 6
  real(prec) :: p
  real(prec) :: dxi,ddxi,xi0
  real(prec), dimension(5) :: xivec, xvec


  real(prec) :: Lref  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: alpha = 16._prec
  real(prec) :: Re    = 32._prec
  real(prec) :: nu
  
  real(prec) :: uex,ux,uxx,uxxx,uxxxx,uxxxxx
  
  nu = Lref*uref/Re

  write(*,*) '**************************************************'
  write(*,*) 'Test of the analytic derivative dTE(i)/dx(i)    '
  write(*,*) 'using finite difference calculations of the exact '
  write(*,*) 'TE.    '
  write(*,*) 'Expected: 2nd order accurate                      '
  write(*,*) '**************************************************'
  write(*,*) ' I  |    dx      |  p  '
  write(*,*) '--------------------------------------------------'

  100 format (I3,e14.4,f7.2)
  
  xivec = (/0._prec,0.25_prec,0.5_prec,0.75_prec,1._prec/)
  xi0 = 0.5_prec
  
  dxi = 0.25_prec
  ddxi = dxi/2._prec


  do n = 1,imax

  !set grid
  xvec = xivec**4
  xip1 = (xivec(3)+ddxi)**4
  xim1 = (xivec(3)-ddxi)**4
  xi = xvec(3)

  !calculate finite difference TE points and 2nd order accurate finite 
  !difference derivative
  xi = xip1
  TEip1 = Burgers_TE(xi,(xvec(4)-xvec(2))/(2._prec*dxi),&
                     (xvec(4)-2._prec*xi+xvec(2))/(dxi**2),&
                     dxi,nu,alpha,Lref)

  xi = xim1
  TEim1 = Burgers_TE(xi,(xvec(4)-xvec(2))/(2._prec*dxi),&
                     (xvec(4)-2._prec*xi+xvec(2))/(dxi**2),&
                     dxi,nu,alpha,Lref)

  dTE_dxi_numerical(1) = (TEip1-Teim1)/(xip1-xim1)



  !calculate analytic derivative
  xi = xvec(3)!+ddxi
  
  
  uex= uxn(xi,nu,alpha,Lref,0)
  ux = uxn(xi,nu,alpha,Lref,1)
  uxx = uxn(xi,nu,alpha,Lref,2)
  uxxx = uxn(xi,nu,alpha,Lref,3)
  uxxxx = uxn(xi,nu,alpha,Lref,4)
  uxxxxx = uxn(xi,nu,alpha,Lref,5)
  
  dTE_dxi_analytic(1) = dTE_dxi(uex,     &
                           ux,           &
                           uxx,          &
                           uxxx,         &
                           uxxxx,        &
                           uxxxxx,       &
                           (xvec(4)-xvec(2))/(2._prec*dxi),&
                           (xvec(4)-2._prec*xvec(3)+xvec(2))/(dxi**2),&
                           dxi,                           &
                           nu,                            &
                           alpha,                         &
                           Lref,                          &
                           0)

  !calculate observed order of accuracy
  error(1) = (dTE_dxi_numerical(2)-dTE_dxi_numerical(1))
  !print*, n,ddxi,dTE_dxi_numerical(1),dTE_dxi_analytic(1),&
  !error(2)/error(1)
  !(dTE_dxi_numerical(2)/dTE_dxi_analytic(2))!/    &
  !          (dTE_dxi_numerical(1)-dTE_dxi_analytic(1))

  if (n.gt.1) then

    p = log( (dTE_dxi_numerical(2)-dTE_dxi_analytic(2))/    &
            (dTE_dxi_numerical(1)-dTE_dxi_analytic(1)) )/  &
        log(2._prec)

    !p = log(error(2)/error(1))/log(2._prec)

    write(*,100) n, dxi, p!, TEip1, TEim1

  endif

  !transfer data
  dTE_dxi_numerical(2) = dTE_dxi_numerical(1)
  dTE_dxi_analytic(2) = dTE_dxi_analytic(1)
  error(2) = error(1)
  
  !xivec = (xivec-xi0)/(2._prec)+xi0
  ddxi = ddxi/2._prec

  enddo

  write(*,*) '--------------------------------------------------'
end subroutine



!*******************************************************************************
!*******************************************************************************
!test for dTE_dxi(u,ux,uxx,uxxx,uxxxx,uxxxxx,xip1,xi,xim1,dxi,nu,a,L,flag)
! flag option 0: derivative of TE(xi) wrt xim1 (or TE(xip1) wrt xi)
!*******************************************************************************
!*******************************************************************************
subroutine dTEip1_dxi_test
  use select_precision, only : prec
  implicit none

  real(prec) :: dx,xip1,xi,xim1,xip2,xim2,x,xxi
  real(prec) :: TEip1, TEi, TEim1
  real(prec), dimension(2) :: dTE_dxi_numerical, dTE_dxi_analytic,error
  integer    :: n
  integer    :: imax = 6
  real(prec) :: p
  real(prec) :: dxi,ddxi,xi0
  real(prec), dimension(5) :: xivec, xvec


  real(prec) :: Lref  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: alpha = 16._prec
  real(prec) :: Re    = 32._prec
  real(prec) :: nu
  
  real(prec) :: uex,ux,uxx,uxxx,uxxxx,uxxxxx
  
  nu = Lref*uref/Re

  write(*,*) '**************************************************'
  write(*,*) 'Test of the analytic derivative dTE(i+1)/dx(i)    '
  write(*,*) 'using finite difference calculations of the exact '
  write(*,*) 'TE.    '
  write(*,*) 'Expected: 2nd order accurate                      '
  write(*,*) '**************************************************'
  write(*,*) ' I  |    dx      |  p  '
  write(*,*) '--------------------------------------------------'

  100 format (I3,e14.4,f7.2)
  
  xivec = (/0._prec,0.25_prec,0.5_prec,0.75_prec,1._prec/)
  xi0 = 0.5_prec
  
  dxi = 0.25_prec
  ddxi = dxi/2._prec


  do n = 1,imax

  !set grid
  xvec = xivec**4
  xip1 = (xivec(3)+ddxi)**4
  xim1 = (xivec(3)-ddxi)**4
  xi = xvec(3)

  !calculate finite difference TE points and 2nd order accurate finite 
  !difference derivative
  xi = xip1
  TEip1 = Burgers_TE(xvec(4),(xvec(5)-xi)/(2._prec*dxi),&
                     (xvec(5)-2._prec*xvec(4)+xi)/(dxi**2),&
                     dxi,nu,alpha,Lref)

  xi = xim1
  TEim1 = Burgers_TE(xvec(4),(xvec(5)-xi)/(2._prec*dxi),&
                     (xvec(5)-2._prec*xvec(4)+xi)/(dxi**2),&
                     dxi,nu,alpha,Lref)

  xi = xvec(3)
  TEi = Burgers_TE(xvec(4),(xvec(5)-xi)/(2._prec*dxi),&
                     (xvec(5)-2._prec*xvec(4)+xi)/(dxi**2),&
                     dxi,nu,alpha,Lref)

  dTE_dxi_numerical(1) = (TEip1-Teim1)/(xip1-xim1)



  !calculate analytic derivative
  xi = xvec(4)
  
  
  uex= uxn(xi,nu,alpha,Lref,0)
  ux = uxn(xi,nu,alpha,Lref,1)
  uxx = uxn(xi,nu,alpha,Lref,2)
  uxxx = uxn(xi,nu,alpha,Lref,3)
  uxxxx = uxn(xi,nu,alpha,Lref,4)
  uxxxxx = uxn(xi,nu,alpha,Lref,5)
  
  dTE_dxi_analytic(1) = dTE_dxi(uex,     &
                           ux,           &
                           uxx,          &
                           uxxx,         &
                           uxxxx,        &
                           uxxxxx,       &
                           (xvec(5)-xvec(3))/(2._prec*dxi),&
                           (xvec(5)-2._prec*xvec(4)+xvec(3))/(dxi**2),&
                           dxi,                           &
                           nu,                            &
                           alpha,                         &
                           Lref,                          &
                           1)

  !calculate observed order of accuracy
  error(1) = (dTE_dxi_numerical(2)-dTE_dxi_numerical(1))
  !print*, n,ddxi,dTE_dxi_numerical(1),dTE_dxi_analytic(1),&
  !error(2)/error(1)
  !(dTE_dxi_numerical(2)-dTE_dxi_analytic(2))/    &
  !          (dTE_dxi_numerical(1)-dTE_dxi_analytic(1))

  if (n.gt.1) then

    p = log( (dTE_dxi_numerical(2)-dTE_dxi_analytic(2))/    &
            (dTE_dxi_numerical(1)-dTE_dxi_analytic(1)) )/  &
        log(2._prec)

    !p = log(error(2)/error(1))/log(2._prec)

    write(*,100) n, dxi, p!, TEip1, TEim1

  endif

  !transfer data
  dTE_dxi_numerical(2) = dTE_dxi_numerical(1)
  dTE_dxi_analytic(2) = dTE_dxi_analytic(1)
  error(2) = error(1)
  
  !xivec = (xivec-xi0)/(2._prec)+xi0
  ddxi = ddxi/2._prec

  enddo

  write(*,*) '--------------------------------------------------'
end subroutine


!*******************************************************************************
!*******************************************************************************
!test for dTE_dxi(u,ux,uxx,uxxx,uxxxx,uxxxxx,xip1,xi,xim1,dxi,nu,a,L,flag)
! flag option 0: derivative of TE(xi) wrt xim1 (or TE(xip1) wrt xi)
!*******************************************************************************
!*******************************************************************************
subroutine dTEim1_dxi_test
  use select_precision, only : prec
  implicit none

  real(prec) :: dx,xip1,xi,xim1,xip2,xim2,x,xxi
  real(prec) :: TEip1, TEi, TEim1
  real(prec), dimension(2) :: dTE_dxi_numerical, dTE_dxi_analytic,error
  integer    :: n
  integer    :: imax = 6
  real(prec) :: p
  real(prec) :: dxi,ddxi,xi0
  real(prec), dimension(5) :: xivec, xvec


  real(prec) :: Lref  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: alpha = 16._prec
  real(prec) :: Re    = 32._prec
  real(prec) :: nu
  
  real(prec) :: uex,ux,uxx,uxxx,uxxxx,uxxxxx
  
  nu = Lref*uref/Re

  write(*,*) '**************************************************'
  write(*,*) 'Test of the analytic derivative dTE(i-1)/dx(i)    '
  write(*,*) 'using finite difference calculations of the exact '
  write(*,*) 'TE.    '
  write(*,*) 'Expected: 2nd order accurate                      '
  write(*,*) '**************************************************'
  write(*,*) ' I  |    dx      |  p  '
  write(*,*) '--------------------------------------------------'

  100 format (I3,e14.4,f7.2)
  
  xivec = (/0._prec,0.25_prec,0.5_prec,0.75_prec,1._prec/)
  xi0 = 0.5_prec
  
  dxi = 0.25_prec
  ddxi = dxi/2._prec


  do n = 1,imax

  !set grid
  xvec = xivec**1
  xip1 = (xivec(3)+ddxi)**1
  xim1 = (xivec(3)-ddxi)**1
  xi = xvec(3)

  !calculate finite difference TE points and 2nd order accurate finite 
  !difference derivative
  xi = xip1
  TEip1 = Burgers_TE(xvec(2),(xi-xvec(1))/(2._prec*dxi),&
                     (xi-2._prec*xvec(2)+xvec(1))/(dxi**2),&
                     dxi,nu,alpha,Lref)

  xi = xim1
  TEim1 = Burgers_TE(xvec(2),(xi-xvec(1))/(2._prec*dxi),&
                     (xi-2._prec*xvec(2)+xvec(1))/(dxi**2),&
                     dxi,nu,alpha,Lref)

  xi = xvec(3)
  TEi = Burgers_TE(xvec(2),(xi-xvec(1))/(2._prec*dxi),&
                     (xi-2._prec*xvec(2)+xvec(1))/(dxi**2),&
                     dxi,nu,alpha,Lref)

  dTE_dxi_numerical(1) = (TEip1-Teim1)/(xip1-xim1)



  !calculate analytic derivative
  xi = xvec(2)


  uex= uxn(xi,nu,alpha,Lref,0)
  ux = uxn(xi,nu,alpha,Lref,1)
  uxx = uxn(xi,nu,alpha,Lref,2)
  uxxx = uxn(xi,nu,alpha,Lref,3)
  uxxxx = uxn(xi,nu,alpha,Lref,4)
  uxxxxx = uxn(xi,nu,alpha,Lref,5)
  
  dTE_dxi_analytic(1) = dTE_dxi(uex,     &
                           ux,           &
                           uxx,          &
                           uxxx,         &
                           uxxxx,        &
                           uxxxxx,       &
                           (xvec(5)-xvec(3))/(2._prec*dxi),&
                           (xvec(5)-2._prec*xvec(4)+xvec(3))/(dxi**2),&
                           dxi,                           &
                           nu,                            &
                           alpha,                         &
                           Lref,                          &
                           -1)

  !calculate observed order of accuracy
  error(1) = (dTE_dxi_numerical(2)-dTE_dxi_numerical(1))
  !print*, n,ddxi,dTE_dxi_numerical(1),dTE_dxi_analytic(1),&
  !error(2)/error(1)
  !(dTE_dxi_numerical(2)-dTE_dxi_analytic(2))/    &
  !          (dTE_dxi_numerical(1)-dTE_dxi_analytic(1))

  if (n.gt.1) then

    p = log( (dTE_dxi_numerical(2)-dTE_dxi_analytic(2))/    &
            (dTE_dxi_numerical(1)-dTE_dxi_analytic(1)) )/  &
        log(2._prec)

    !p = log(error(2)/error(1))/log(2._prec)

    write(*,100) n, dxi, p!, TEip1, TEim1

  endif

  !transfer data
  dTE_dxi_numerical(2) = dTE_dxi_numerical(1)
  dTE_dxi_analytic(2) = dTE_dxi_analytic(1)
  error(2) = error(1)
  
  !xivec = (xivec-xi0)/(2._prec)+xi0
  ddxi = ddxi/2._prec

  enddo

  write(*,*) '--------------------------------------------------'
end subroutine




end module burgers1d_testroutines


