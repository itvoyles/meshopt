module functional


contains

function calc_J(x,imax,nu,alpha,Lref)
  use select_precision, only : prec
  use burgers1d_functions
  use derivatives
  use functional_integral_te_squared
  implicit none
  
  integer, intent(in) :: imax
  real(prec), intent(in) :: nu,alpha,Lref
  real(prec), dimension(imax), intent(in) :: x
  real(prec) :: calc_J
  
  real(prec), dimension(imax) :: TE
  
  
  TE = Truncation_error(x,nu,alpha,Lref,imax)
  calc_J = functional_J(x,TE,imax)



end function calc_J




!*******************************************************************************
!*******************************************************************************
! Subroutine to calculate the functional sensitivity
!*******************************************************************************
!*******************************************************************************
subroutine calc_dJdx(djdx,x,imax,nu,alpha,Lref)
  use select_precision, only : prec
  use burgers1d_functions, only : dTE_dxi, burgers_te, uxn
  use fd_derivative_calc
  !use functional_integral_TE_squared
  use functional_integral_TE_4th
  !use functional_integral_te_1st
  implicit none
  
  integer, intent(in) :: imax
  real(prec), dimension(imax),intent(in) :: x
  real(prec), intent(in)                 :: nu, alpha, Lref
  real(prec), dimension(imax),intent(out)::djdx

  
  real(prec), dimension(imax) :: TE, dTEi_dxvec, dTEip1_dxvec, dTEim1_dxvec
  real(prec), dimension(imax) :: x_xi, x_xixi
  real(prec)                  :: dxi
  real(prec) :: uex, ux,uxx,uxxx,uxxxx,uxxxxx
  integer    :: i
  
  !compuational space grid spacing
  dxi = 1._prec

  !calculate grid metrics ------------------------------------------------------
  call dnfdxn(x_xi,x,dxi,1,2,imax)
  call dnfdxn(x_xixi,x,dxi,2,2,imax)
  
  !Derived truncation error
  TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref)
  !TE = Burgers_fd_cont_resid(x,imax,nu,alpha,lref)
  TE(1) = 0._prec
  TE(imax) = 0._prec

  !Derivative of derived truncation error
  dTEi_dxvec = 0._prec
  dTEip1_dxvec = 0._prec
  dTEim1_dxvec = 0._prec
  do i = 2,imax-1
    uex=     uxn(x(i),nu,alpha,Lref,0)
    ux =     uxn(x(i),nu,alpha,Lref,1)
    uxx =    uxn(x(i),nu,alpha,Lref,2)
    uxxx =   uxn(x(i),nu,alpha,Lref,3)
    uxxxx =  uxn(x(i),nu,alpha,Lref,4)
    uxxxxx = uxn(x(i),nu,alpha,Lref,5)

    dTEi_dxvec(i) = dTE_dxi(uex,             &
                           ux,               &
                           uxx,              &
                           uxxx,             &
                           uxxxx,            &
                           uxxxxx,           &
                           x_xi(i),          &
                           x_xixi(i),        &
                           dxi,              &
                           nu,               &
                           alpha,            &
                           Lref,             &
                           0)

    uex=     uxn(x(i+1),nu,alpha,Lref,0)
    ux =     uxn(x(i+1),nu,alpha,Lref,1)
    uxx =    uxn(x(i+1),nu,alpha,Lref,2)
    uxxx =   uxn(x(i+1),nu,alpha,Lref,3)
    uxxxx =  uxn(x(i+1),nu,alpha,Lref,4)
    uxxxxx = uxn(x(i+1),nu,alpha,Lref,5)
    
    dTEip1_dxvec(i) = dTE_dxi(uex,           &
                           ux,               &
                           uxx,              &
                           uxxx,             &
                           uxxxx,            &
                           uxxxxx,           &
                           x_xi(i+1),        &
                           x_xixi(i+1),      &
                           dxi,              &
                           nu,               &
                           alpha,            &
                           Lref,             &
                           1)


    uex=     uxn(x(i-1),nu,alpha,Lref,0)
    ux =     uxn(x(i-1),nu,alpha,Lref,1)
    uxx =    uxn(x(i-1),nu,alpha,Lref,2)
    uxxx =   uxn(x(i-1),nu,alpha,Lref,3)
    uxxxx =  uxn(x(i-1),nu,alpha,Lref,4)
    uxxxxx = uxn(x(i-1),nu,alpha,Lref,5)
    
    dTEim1_dxvec(i) = dTE_dxi(uex,           &
                           ux,               &
                           uxx,              &
                           uxxx,             &
                           uxxxx,            &
                           uxxxxx,           &
                           x_xi(i-1),        &
                           x_xixi(i-1),      &
                           dxi,              &
                           nu,               &
                           alpha,            &
                           Lref,             &
                           -1)

                                           
  enddo                     



  call analytic_dJdx(djdx,x,&
                     TE,dTEi_dxvec,dTEip1_dxvec,dTEim1_dxvec,&
                     imax)

end subroutine calc_dJdx


!*******************************************************************************
!*******************************************************************************
! test of functional sensitivity calculation
!*******************************************************************************
!*******************************************************************************
subroutine test_J_eq_integral_TE_squared
  use select_precision, only : prec
  use burgers1d_functions, only : dTE_dxi, burgers_te, uxn
  use fd_derivative_calc
  !use functional_integral_TE_squared  
  use functional_integral_TE_4th 
  !use functional_integral_te_1st
  implicit none

  real(prec) :: dx,xip1,xim1
  integer,parameter    :: imax = 17
  real(prec), dimension(imax) :: x,xtemp,xi,TE, TEip1,TEim1
  real(prec), dimension(imax) :: dTEi_dxvec,dTEip1_dxvec,dTEim1_dxvec
  real(prec), dimension(imax) :: dJ_dx,x_xi,x_xixi
  real(prec) :: xi0,x0
  real(prec) :: Jip1, Ji, Jim1

  real(prec), dimension(2) :: dJ_dx_numerical, dJ_dx_analytic
  integer    :: n,i,j
  integer    :: nmax = 4
  real(prec) :: p
  real(prec) :: dxi,ddxi
  real(prec), dimension(imax) :: temp
  integer :: order = 2
  real(prec) :: Lref  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: alpha = 16._prec
  real(prec) :: Re    = 32._prec
  real(prec) :: nu
  integer    :: pwr=1
  
  real(prec) :: uex,ux,uxx,uxxx,uxxxx,uxxxxx

  nu = Lref*uref/Re

  write(*,*) '**************************************************'
  write(*,*) 'Test of the analytic derivative dJ/dx(i) using a  '
  write(*,*) 'finite difference calculation of the exact TE.    '
  write(*,*) '**************************************************'





  100 format (I3,e14.4,f7.2,2e12.4)
  
  
  j = (imax-1)/2+1 !index to test J
  !do j = 2,2
  do j = 2,imax-1
  
  print*,
  write(*,*) '--------------------------------------------------'
  write(*,'(X,A,I4)') 'node:',j
  write(*,*) '--------------------------------------------------'
  write(*,*) ' I  |    dx      |  p  | numerical | analytic  '
  write(*,*) '--------------------------------------------------'

  
  
  dxi = 1._prec/real(imax-1,prec)
  !create original grid
  do i = 1,imax
    xi(i) = (0.5_prec-real(i-1,prec)*dxi)*8._prec
  enddo
  ddxi = dxi/1._prec
  x0 = xi0
  x = xi**pwr
  
  
  !start do loop for order of accuracy checks **********************************
  do n = 1,nmax
  x = xi**pwr
  ddxi = ddxi/2._prec


  !calculate analytic sensitivies
  call calc_dJdx(temp,x,imax,nu,alpha,Lref)
  dJ_dx_analytic(1) = temp(j)

  !write(*,*) imax,nu,alpha,Lref
  !write(*,'(1f8.4)')x
  

!functional with perturbed point by +ddxi
!calculate grid metrics ------------------------------------------------------
  xip1 = (xi(j)+ddxi)**pwr
  x(j) = xip1
  call dnfdxn(x_xi,x,dxi,1,order,imax)
  call dnfdxn(x_xixi,x,dxi,2,order,imax)
  !Derived truncation error
  TEip1 = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref)
  TEip1(1) = 0._prec
  TEip1(imax) = 0._prec
  Jip1 = functional_J(x,TEip1,imax)

!functional with perturbed point by -ddxi
!calculate grid metrics ------------------------------------------------------
  xim1 = (xi(j)-ddxi)**pwr
  x(j) = xim1
  call dnfdxn(x_xi,x,dxi,1,order,imax)
  call dnfdxn(x_xixi,x,dxi,2,order,imax)
  !Derived truncation error
  TEim1 = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref)
  TEim1(1) = 0._prec
  TEim1(imax) = 0._prec
  Jim1 = functional_J(x,TEim1,imax)

  dJ_dx_numerical(1) = (Jip1-Jim1)/(xip1-xim1) !2th order term

  !calculate observed order of accuracy
  if (n.gt.1) then 
    p = log( (dJ_dx_analytic(2)-dJ_dx_numerical(2))/   &
             (dJ_dx_analytic(1)-dJ_dx_numerical(1)))/  &
       log(2._prec)

    write(*,100) n, ddxi, p,dJ_dx_numerical(1),dJ_dx_analytic(1)
  endif

  !transfer data
  dJ_dx_numerical(2) = dJ_dx_numerical(1)
  dJ_dx_analytic(2) = dJ_dx_analytic(1)



  enddo
  
  write(*,*) '--------------------------------------------------'
  enddo

end subroutine test_J_eq_integral_TE_squared


end module functional
