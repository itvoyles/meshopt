!Module of functions related to steady, 1D Burgers' Equation
module burgers1d_functions


contains
!*******************************************************************************
!*******************************************************************************
!           EXACT SOLUTION AND DERIVATIVES TO 1D BURGERS' EQUATION
!*******************************************************************************
!*******************************************************************************
  !exact solution to 1D burgers' equation
  ! u = -2va/L*tanh(xa/L)
  !
  !derivatives
  ! ux    =  2.0v(a/L)^2(tanh(ax/L)^2-1.0)
  ! uxx   = -4.0v(a/L)^3*tanh(ax/L)*(tanh(ax/L)^3-1.0)
  ! uxxx  =  4.0v(a/L)^4*(tanh(ax/L)^2-1.0)
  !         +8.0v(a/L)^4*tanh(ax/L)^2*(tanh(ax/L)^2-1.0)
  ! uxxxx = -32.0v(a/L)^5*tanh(ax/L)*(tanh(ax/L)^2-1.0)^2
  !         -16.0v(a/L)^5*tanh(ax/L)^3*(tanh(ax/L)^2-1.0)
  ! INPUTS
  !   u => exact solution
  !   v => viscosity
  !   a => scaling factor
  !   L => domain length
  elemental function uexact(x,v,a,L)
    use select_precision, only : prec
    implicit none
    real(prec), intent(in):: x,v,a,L
    real(prec)            :: uexact

    uexact = -2._prec*v*a/L*tanh(x*a/L)
  end function uexact

  elemental function uxn(x,v,a,L,flag)
  !up to 5th derivative of burgers' exact solution
  !tested using finite difference, derivative 1-5 are verified accurate.
    use select_precision, only : prec
    implicit none

    real(prec), intent(in) :: x,v,a,L
    integer, intent(in) :: flag
    real(prec)             :: uxn

    real(prec) :: c1,c2,f,f1,f2,f3,f4,f5

    ! f = c1*tanh(c2x)
    c1 = -2._prec*v*a/L
    c2 = a/L

    f = c1*tanh(c2*x)
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

    select case(flag)
      case(1)
        uxn = f1
      case(2)
        uxn = f2
      case(3)
        uxn = f3
      case(4)
        uxn = f4
      case(5)
        uxn = f5
    end select

  end function
  

  

!*******************************************************************************
!*******************************************************************************
! Burgers' Truncation Error with discrete metrics
!*******************************************************************************
!*******************************************************************************
function Truncation_error(x,v,a,L,imax)
  use select_precision, only : prec
  use derivatives
  implicit none
  
  integer, intent(in) :: imax
  real(prec), intent(in) :: v,a,L
  real(prec), dimension(imax) :: x
  real(prec), dimension(imax) :: truncation_error
  
  real(prec), dimension(imax) :: x_xi, x_xixi
  
  
  x_xi =  dudx(x,1._prec,imax)
  x_xixi =  d2udx2(x,1._prec,imax)
  truncation_error = Burgers_TE(x,x_xi,x_xixi,1._prec,v,a,L)
  
  
end function


!verified 4th order accurate on generic grid with transformation x = xi^4
elemental function Burgers_TE(x,x_xi,x_xixi,dxi,v,a,L,f)
  use select_precision, only : prec
  implicit none

  real(prec), intent(in) :: x,x_xi,x_xixi,dxi !grid related inputs
  real(prec), intent(in) :: v,a,L             !solution specifi inputs
  integer, intent(in), optional :: f          !flag to return specific terms
  real(prec) :: Burgers_Te

  real(prec) :: TE_stretch, Te_stand, Te_mixed
  real(prec) :: u,uxx,uxxx,uxxxx
  integer :: flag

  if (present(f)) then
    flag=f
  else
    flag = 0
  endif

  u = uxn(x,v,a,L,0)
  uxx = uxn(x,v,a,L,2)
  uxxx = uxn(x,v,a,L,3)
  uxxxx = uxn(x,v,a,L,4)
  
  !standard term
  TE_stand = dxi**2*x_xi**2*(&
                              1._prec/6._prec*u*uxxx    &
                             -1._prec/12._prec*v*uxxxx  &
                              )

  TE_stretch = dxi**2*x_xixi*(                                             &
                               -1._prec/3._prec*v*uxxx                     &
                               +0.5_prec*u*uxx                             &
                               )
 
  !mixed term
  TE_mixed = dxi**2*v*(x_xixi/x_xi)**2*uxx/4._prec


  If (flag .eq. 1) then
    Burgers_TE = TE_stand
  elseif (flag .eq. 2) then
    Burgers_TE = TE_stretch
  elseif (flag .eq. 3) then
    Burgers_TE = TE_mixed
  else
    Burgers_TE = TE_stand+Te_stretch+TE_mixed
  endif


end function Burgers_TE




!*******************************************************************************
!*******************************************************************************
!                    function TE, dTE/dx, and test routines
!*******************************************************************************
!*******************************************************************************
!verified accurate by comparing to second-order accurate finite difference
elemental function dTE_dxi(u,ux,uxx,uxxx,uxxxx,uxxxxx,Xxi,Xxixi,&
                           dxi,nu,a,L,flag)
! dTE_dxi = pdTE_pdxi + dTE/dAi * dAi/dxi

  use select_precision, only : prec

  implicit none
  integer, intent(in)    :: flag !0 for dTE(x_i)/dxi, 1 for dTE(x_i+1)/dxi
  real(prec), intent(in) :: nu, a, L !burgers solution specifics
  real(prec), intent(in) :: Xxi,Xxixi, dxi !grid points required for derivative
  real(prec), intent(in) :: u,ux,uxx,uxxx,uxxxx,uxxxxx !solution and derivatives
  
  real(prec) :: dTE_dxi !derivative of truncation error with respect to x

  real(prec) :: pdTE_pdx, A1, A2, A3
  real(prec) :: dA1_dx, dA2_dx, dA3_dx
  real(prec) :: dTE_dA1, dTE_dA2, dTE_dA3

  real(prec) :: dXxi_dx,dXxixi_dx ! grid metrics

  !calculate metric derivatives based on flag
  !solution derivatives are calculated outside of subroutine
  select case(flag)
    case(0) !dTE(x_i)/dxi ******************************************************
    !solution and derivatives evaluated at xi
    ! verified accurate with finite difference

      !calculate metric derivatives
      dXxi_dx = 0._prec/(2._prec*dxi)
      dXxixi_dx = -2._prec/dxi**2

      !calculate coefficients and derivatives
      A1 = (&
            1._prec/6._prec*u*uxxx     &
           -1._prec/12._prec*nu*uxxxx  &
            )

      dA1_dx = (1._prec/6._prec*(ux*uxxx+u*uxxxx) &
               -1._prec/12._prec*nu*uxxxxx)

      A2 = (&
           -1._prec/3._prec*nu*uxxx                    &
           +1._prec/2._prec*u*uxx                      &
           )

      dA2_dx = (1._prec/2._prec*(ux*uxx+u*uxxx) &
               -1._prec/3._prec*nu*uxxxx)

      A3 = nu*uxx/4._prec

      dA3_dx = (1._prec/4._prec*nu*uxxx)

      !calculate the partial of TE with respect to x
      pdTE_pdx = A1*2._prec*Xxi*dXxi_dx                                        &
                +A2*dXxixi_dx                                                  &
                +A3*2._prec*(Xxixi/Xxi)*(dXxixi_dx/Xxi-Xxixi/Xxi**2*dXxi_dx)

      dTE_dA1 = Xxi**2

      dTE_dA2 = Xxixi

      dTE_dA3 = (Xxixi/Xxi)**2

      dTE_dxi = dxi**2*(pdTE_pdx                  &
                       +dTE_dA1*dA1_dx           &
                       +dTE_dA2*dA2_dx            &
                       +dTE_dA3*dA3_dx)


    case(1) !dTE(x_i+1)/dxi ****************************************************
    ! solution and derivatives evaluated at xip1
    !verified accurate with finite difference
    
      !calculate metric derivatives
      dXxi_dx = -1._prec/(2._prec*dxi)
      dXxixi_dx = 1._prec/dxi**2

      !calculate coefficients and derivatives
      A1 = (&
            1._prec/6._prec*u*uxxx     &
           -1._prec/12._prec*nu*uxxxx  &
            )

      dA1_dx = 0._prec 

      A2 = (&
           -1._prec/3._prec*nu*uxxx                    &
           +1._prec/2._prec*u*uxx                      &
           )

      dA2_dx = 0._prec 

      A3 = nu*uxx/4._prec

      dA3_dx = 0._prec!(1._prec/4._prec*nu*uxxx)

      !calculate the partial of TE with respect to x
      pdTE_pdx = A1*2._prec*Xxi*dXxi_dx                                        &
                +A2*dXxixi_dx                                                  &
                +A3*2._prec*(Xxixi/Xxi)*(dXxixi_dx/Xxi-Xxixi/Xxi**2*dXxi_dx)

      dTE_dA1 = Xxi**2

      dTE_dA2 = Xxixi

      dTE_dA3 = (Xxixi/Xxi)**2

      dTE_dxi = dxi**2*(pdTE_pdx                  &
                       +dTE_dA1*dA1_dx           &
                       +dTE_dA2*dA2_dx            &
                       +dTE_dA3*dA3_dx)
    
    
    case(-1) !dTE(x_i+1)/dxi ****************************************************
    ! solution and derivatives evaluated at xip1
    !verified accurate with finite difference
    
      !calculate metric derivatives
      dXxi_dx = 1._prec/(2._prec*dxi)
      dXxixi_dx = 1._prec/dxi**2

      !calculate coefficients and derivatives
      A1 = (&
            1._prec/6._prec*u*uxxx     &
           -1._prec/12._prec*nu*uxxxx  &
            )

      dA1_dx = 0._prec 

      A2 = (&
           -1._prec/3._prec*nu*uxxx                    &
           +1._prec/2._prec*u*uxx                      &
           )

      dA2_dx = 0._prec 

      A3 = nu*uxx/4._prec

      dA3_dx = 0._prec!(1._prec/4._prec*nu*uxxx)

      !calculate the partial of TE with respect to x
      pdTE_pdx = A1*2._prec*Xxi*dXxi_dx                                        &
                +A2*dXxixi_dx                                                  &
                +A3*2._prec*(Xxixi/Xxi)*(dXxixi_dx/Xxi-Xxixi/Xxi**2*dXxi_dx)

      dTE_dA1 = Xxi**2

      dTE_dA2 = Xxixi

      dTE_dA3 = (Xxixi/Xxi)**2

      dTE_dxi = dxi**2*(pdTE_pdx                  &
                       +dTE_dA1*dA1_dx           &
                       +dTE_dA2*dA2_dx            &
                       +dTE_dA3*dA3_dx)
                                              
  end select
  



end function dTE_dxi


!******************************************************************************
!******************************************************************************
! Calculates the exact truncation error for finite difference Burgers' equation
!******************************************************************************
!******************************************************************************
function Burgers_fd_cont_resid(x,imax,v,a,l)
  use select_precision
  use fd_derivative_calc
  use derivative_transform
  implicit none
  
  integer, intent(in) :: imax
  real(prec), intent(in) :: v,a,l
  real(prec), dimension(imax), intent(in) :: x
  real(prec), dimension(imax):: Burgers_fd_cont_resid
  
  real(prec), dimension(imax) :: u,u_xi,u_xixi,u_x,u_xx,x_xi,x_xixi
  real(prec) :: dxi
  integer :: order
  
  dxi = 1._prec !spacing in compuational space
  order = 2     !order of accuracy of derivative to calculate
  
  !calculate grid metrics ------------------------------------------------------
  call dnfdxn(x_xi,x,dxi,1,order,imax)
  call dnfdxn(x_xixi,x,dxi,2,order,imax)


  !calculate exact solution and derivatives ------------------------------------
  u      = uxn(x,v,a,l,0)
  
  ! derivatives in computational space
  call dnfdxn(u_xi,u,dxi,1,order,imax)
  call dnfdxn(u_xixi,u,dxi,2,order,imax)


  ! calculate the continuous residual for finite difference method
  Burgers_fd_cont_resid = ((u/x_xi) + v*(x_xixi/x_xi**3))*u_xi &
                        - (v/x_xi**2)*u_xixi


end function 







end module burgers1d_functions


