!*******************************************************************************
!*******************************************************************************
!                         Functional J and dJ/dx
!*******************************************************************************
!*******************************************************************************
module functional_integral_te_4th
use select_precision, only : prec
!  
!      N-1  
!      --
!  J = \  f(i)
!      /
!      --
!      i=1
!
!          1 /                     \/               \
!  f(i) = ---| TE(i+1)^4 + TE(i)^4 || x(i+1) - x(i) |
!          2 \                     /\               /
!
!  dJ     -1  /                       \ 
! ----- = --- | TE(i+1)^4 - TE(i-1)^4 | 
! dx(i)    2  \                       /
!
!          /         dTE(i+1)        dTE(i)  \/             \
!        + | TE(i+1) --------- + TE(i)------ || x(i+1)-x(i) |
!          \           dx              dx    /\             /
!
!          /        dTE(i)           dTE(i-1) \/             \
!        + | TE(i) -------- + TE(i-1)-------- || x(i)-x(i-1) |
!          \          dx                dx    /\             /
!
!          /          dTE(i+1)  \/               \
!        + | TE(i+1) ---------- || x(i+2)-x(i+1) |
!          \            dx      /\               /
!
!          /          dTE(i-1)  \/               \
!        + | TE(i-1) ---------- || x(i-1)-x(i-2) |
!          \            dx      /\               /

integer, private, parameter :: pwr = 2
real(prec), private :: norm_coef = 1._prec!20._prec**16
contains



function functional_J(x,TE,imax)
  use select_precision, only : prec
  use misc_func, only : trap_sum
  implicit none
  
  integer, intent(in) :: imax
  real(prec), dimension(imax), intent(in) :: x, TE
  real(prec) :: functional_J
  
  functional_J = trap_sum(x,TE**pwr,imax)*norm_coef

  
end function functional_J


subroutine analytic_dJdx(djdx,x,TE,dTEi_dx,dTEip1_dx,dTEim1_dx,imax)
  use select_precision, only : prec
  implicit none
  
  integer, intent(in) :: imax
  real(prec),dimension(imax),intent(in)   :: TE
  real(prec),dimension(imax),intent(in)   :: dTEi_dx,dTEip1_dx,dTEim1_dx
  real(prec),dimension(imax),intent(in)   :: x
  !real(prec),dimension(imax)   :: dJdx_eq_integral_TE_squared
  real(prec), dimension(imax), intent(out):: djdx
  real(prec),dimension(imax) :: pdf_pdx
  integer :: i
  real(prec) :: two = 2._prec
  real(prec) :: p1,p2,p3
  

  do i = 2,imax-1

    pdf_pdx(i) =-(TE(i+1)**pwr-TE(i-1)**pwr)/2._prec &
                +(pwr*TE(i+1)**(pwr-1)*dTEip1_dx(i)&
                +pwr*TE(i)**(pwr-1)*dTEi_dx(i))*(x(i+1)-x(i))/2._prec&
                +(pwr*TE(i)**(pwr-1)*dTEi_dx(i)&
                +pwr*TE(i-1)**(pwr-1)*dTEim1_dx(i))*(x(i)-x(i-1))/2._prec

    if (i.gt.2) then
      pdf_pdx(i) = pdf_pdx(i)+pwr*TE(i-1)**(pwr-1)*dTEim1_dx(i)*(x(i-1)-x(i-2))/2._prec
    endif
    
    if (i.lt.imax-1) then
      pdf_pdx(i) = pdf_pdx(i)+pwr*TE(i+1)**(pwr-1)*dTEip1_dx(i)*(x(i+2)-x(i+1))/2._prec
    endif

  enddo

  djdx = 0._prec
  djdx(2:imax-1) = pdf_pdx(2:imax-1)*norm_coef


end subroutine analytic_dJdx

end module functional_integral_TE_4th
