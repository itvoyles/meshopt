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
!  dJ     -1  /                       \ 
! ----- = --- | TE(i+1)^2 - TE(i-1)^2 | 
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

contains

function functional_J(x,TE,imax)
  use select_precision, only : prec
  use misc_func, only : trap_sum
  implicit none
  
  integer, intent(in) :: imax
  real(prec), dimension(imax), intent(in) :: x, TE
  real(prec) :: functional_J
  
  functional_J = trap_sum(x,TE**2,imax)
  
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

    pdf_pdx(i) =-(TE(i+1)**2-TE(i-1)**2)/2._prec &
                +(TE(i+1)*dTEip1_dx(i)+TE(i)*dTEi_dx(i))*(x(i+1)-x(i))&
                +(TE(i)*dTEi_dx(i)+TE(i-1)*dTEim1_dx(i))*(x(i)-x(i-1))  

    if (i.gt.2) then
      pdf_pdx(i) = pdf_pdx(i)+TE(i-1)*dTEim1_dx(i)*(x(i-1)-x(i-2))
    endif
    
    if (i.lt.imax-1) then
      pdf_pdx(i) = pdf_pdx(i)+TE(i+1)*dTEip1_dx(i)*(x(i+2)-x(i+1))
    endif

  enddo

  djdx = 0._prec
  djdx(2:imax-1) = pdf_pdx(2:imax-1)


end subroutine analytic_dJdx

end module functional_integral_TE_squared
