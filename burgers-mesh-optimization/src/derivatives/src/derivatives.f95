module derivatives
  use select_precision, only : prec
  use fd_derivatives

contains
function dudx(x,dx,imax)
  implicit none

  integer, intent(in) :: imax
  real(prec),intent(in) :: dx
  real(prec), dimension(imax),intent(in) :: x
  real(prec), dimension(imax) :: dudx

  dudx(2:imax-1) = (x(3:imax)-x(1:imax-2))/(2._prec*dx)
  dudx(1) = (-3._prec*x(1)+4._prec*x(2)-x(3))/(2._prec*dx)
  dudx(imax) = (3._prec*x(imax)-4._prec*x(imax-1)+x(imax-2))/2._prec*dx


end function dudx

function d2udx2(x,dx,imax)
  implicit none

  integer, intent(in) :: imax
  real(prec),intent(in) :: dx
  real(prec), dimension(imax),intent(in) :: x
  real(prec), dimension(imax) :: d2udx2

  d2udx2(2:imax-1) = (x(3:imax)-2._prec*x(2:imax-1)+x(1:imax-2))/dx**2
  d2udx2(1) = (2._prec*x(1)-5._prec*x(2)+4._prec*x(3)-x(4))/dx**2
  d2udx2(imax) = (2._prec*x(imax)-5._prec*x(imax-1)+4._prec*x(imax-2)-x(imax-3))/dx**2

end function



end module
