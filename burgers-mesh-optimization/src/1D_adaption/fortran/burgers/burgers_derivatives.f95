	!Module of functions related to steady, 1D Burgers' Equation
module burgers1d_functions
!Burgers derivatives, but replaced by a single better generic function, uxn()

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

  elemental function ux(x,v,a,L)
    use select_precision, only : prec
    implicit none
    real(prec), intent(in):: x,v,a,L
    real(prec)            :: ux

    ux = 2.0_prec*v*(a/L)**2*(tanh(a*x/L)**2 - 1.0_prec)
  end function

  elemental function uxx(x,v,a,L)
    use select_precision, only : prec
    implicit none
    real(prec), intent(in):: x,v,a,L
    real(prec)            :: uxx

    uxx = -4.0_prec*v*(a/L)**3*tanh(a*x/L)&
          *(tanh(a*x/L)**2 - 1.0_prec)
  end function

  function uxxx(x,v,a,L)
    use select_precision, only : prec
    implicit none
    real(prec), intent(in):: x,v,a,L
    real(prec)            :: uxxx

    print*, 'Error! function uxxx in Burgers1d_functions has an error. Use uxn() isntead'
    uxxx=0._prec
    return
    uxxx =  -4.0_prec*v*(a/L)**4*(tanh(a*x/L)**2-1.0_prec) &
           + 8.0_prec*v*(a/L)**4*tanh(a*x/L)**2*(tanh(a*x/L)**2-1.0_prec)
  end function

  elemental function uxxxx(x,v,a,L)
    use select_precision, only : prec
    implicit none
    real(prec), intent(in):: x,v,a,L
    real(prec)            :: uxxxx

    uxxxx =  - 32.0_prec*v*(a/L)**5*tanh(a*x/L)*(tanh(a*x/L)**2-1.0_prec)**2 &
             - 16.0_prec*v*(a/L)**5*tanh(a*x/L)**3*(tanh(a*x/L)**2-1.0_prec)
  end function


  elemental function uxxxxx(x,v,a,L)
  !5th derivative of burgers' exact solution
    use select_precision, only : prec
    implicit none

    real(prec), intent(in) :: x,v,a,L
    real(prec)             :: uxxxxx

    real(prec) :: c1,c2,f,f1,f2,f3,f4,f5

    ! f = c1*tanh(c2x)
    c1 = -2._prec*v*a/L
    c2 = a/L

    f = c1*tanh(c2*x)
    f1 = -c2/c1*(f**2-c1**2)
    f2 = -c2/c1*(2._prec*f*f1)
    f3 = -2._prec*c2/c1*(f1*f1+f*f2)
    f4 = -2._prec*c2/c1*(3._prec*f1*f2+f*f3)
    f5 = -2._prec*c2/c1*(3._prec*(f2**2+f1*f3)+f1*f4+f*f4)

    uxxxxx = f5

  end function


  elemental function uxn(x,v,a,L,flag)
  !5th derivative of burgers' exact solution
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
    f1 = -c2/c1*(f**2-c1**2)
    f2 = -c2/c1*(2._prec*f*f1)
    f3 = -2._prec*c2/c1*(f1*f1+f*f2)
    f4 = -2._prec*c2/c1*(3._prec*f1*f2+f*f3)
    f5 = -2._prec*c2/c1*(3._prec*(f2**2+f1*f3)+f1*f4+f*f4)

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


end module burgers1d_functions





 
