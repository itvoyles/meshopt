module spring_functional
  use functional_integral_te_squared
  use spring_simple
  
  
contains

function calc_J(k,v,a,l,imax,jmax)
  use select_precision, only : prec
  use burgers2d_functions,  only : truncation_error	
  implicit none
  
  integer, intent(in) :: imax,jmax
  real(prec), dimension((imax-1)*j+imax*(jmax-1)), intent(in) :: k
  real(prec), dimension(imax,jmax) :: x,y
  real(prec), intent(in) :: v,a,l
  real(prec) :: calc_J
  
  real(prec), dimension(imax,jmax) :: TE
  
  call springsystem(x,y,k,imax,jmax,Lref)
  
  TE = Truncation_error(x,y,a,v,l,imax,jmax)
  calc_J = functional_J(x,y,TE,imax,jmax)



end function calc_J


function calc_djdx(k,v,a,L,imax,jmax)
  use select_precision, only : prec
  use burgers2d_functions
  implicit none

  integer, intent(in) :: imax,jmax
  real(prec), dimension((imax-1)*j+imax*(jmax-1)), intent(in) :: k
  real(prec), intent(in) :: v,a,L
  real(prec), dimension((imax-1)*j+imax*(jmax-1)) :: calc_djdx
  
  real(prec) :: err_lim !error tolerance (normalized)
  
  integer :: i
  real(prec) :: dx
  real(prec) :: J0

    
  err_lim = 0.0001_prec
  J0 = calc_J(k,imax,jmax,a,v,l)
  calc_djdx = 0._prec
  do i = 1,(imax-1)*j+imax*(jmax-1)
  
    dx = (maxval(k))/100000._prec
    calc_djdx(i) = djdx_loc(k,i,dx,err_lim,v,a,L,imax,jmax)
    
  enddo
  


  !remove NaN
  where (calc_djdx.ne.calc_djdx)
  calc_djdx = 0._prec
  end where

end function calc_djdx


function djdx_loc(k,i,J0,dx,err_lim,v,a,L,imax,jmax)
!given a mesh, this function calculates the mesh sensitivity of the functional
!J by perturbing the mesh by dx. Richardson extrapolation is then used to 
!estimate the error. If the error is below the defined tolerance (err_lim) then
!djdx is returned, if not the required refinement factor is calculated to reach 
!this limit and the derivative is calculated again.
  use select_precision, only : prec
  use burgers2d_functions
  implicit none
  
  integer, intent(in) :: imax,jmax
  integer, intent(in) :: ii
  real(prec) :: J0
  real(prec), dimension((imax-1)*j+imax*(jmax-1)), intent(in) :: k
  real(prec), intent(in) ::dxin,err_lim,v,a,L
  real(prec) :: djdx_loc
  
  real(prec), dimension((imax-1)*j+imax*(jmax-1)) :: kip1
  real(prec) :: Jip1!,Jim1


  dx = dxin
  xip1 = x
  
  !calc 1
  kip1(ii) = k(ii)+dx
  Jip1 = calc_J(kip1,imax,jmax,a,v,l)


  dJdx_loc = (Jip1-J0)/(dx)




end function djdx_loc

end module spring_functional
