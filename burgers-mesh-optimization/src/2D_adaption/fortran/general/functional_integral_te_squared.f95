!*******************************************************************************
!*******************************************************************************
!                         Functional J and dJ/dx
!*******************************************************************************
!*******************************************************************************
module functional_integral_te_squared
!  
!      N-1  
!      --
!  J = \  f(i)
!      /
!      --
!      i=1
!
!          1 /                     \/               \
!  f(i) = ---| TE(i+1)^2 + TE(i)^2 || x(i+1) - x(i) |
!          2 \                     /\               /
!

contains

function functional_J(x,y,TE,imax,jmax)
  use select_precision, only : prec
  use misc_func, only : trap_sum2d
  implicit none
  
  integer, intent(in) :: imax,jmax
  real(prec), dimension(imax,jmax),intent(in) :: x,y, TE
  real(prec) :: functional_J

  functional_J = trap_sum2d(x,y,TE**2,imax,jmax)
  
end function functional_J


end module functional_integral_TE_squared
