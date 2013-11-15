module calc_grid

contains

subroutine calc_initial_grid(xout,imax,a,v,l)
  use select_precision, only : prec
  use fd_derivative_calc
  use burgers1d_functions
  implicit none
  
  integer,  intent(in) :: imax
  real(prec), intent(in) :: a,v,l
  real(prec), dimension(imax), intent(out) :: xout

  real(prec), dimension(imax) :: x,xnew,xold, TEn
  real(prec), dimension(imax-1) :: TE
  real(prec), dimension(imax-1) :: dx, dxold, dxnew
  real(prec) :: Leng
  integer :: i,iter
  real(prec) :: tol
  real(prec) :: dxstd


  dxstd = l/real(imax-1,prec)

  dxold = l/real(imax-1,prec)
  do i = 1,imax
  xold(i) = -l/2._prec+(i-1)*l/real(imax-1,prec)
  enddo



  !calculates the exact truncation error from the finite-difference residual
  TEn = abs(Burgers_fd_cont_resid(xold,imax,v,a,l))
  TE = (TEn(1:imax-1)+TEn(2:imax))/2._prec




  !where (TE.gt.100._prec/dxstd)
  !  TE = 100._prec/dxstd
  !end where
  
  !where (TE.lt.0.01_prec*maxval(TE))
  !  TE = 0.01_prec*maxval(TE)
  !end where



  dx = (1._prec/TE)
  dx(1) = dx(2)
  dx(imax-1) = dx(imax-2)
  where (dx.lt.0.1_prec*maxval(dx))
    dx = 0.1_prec*maxval(dx)
  end where
  
  
  leng = sum(dx)
  dx = dx/leng
  




  xnew(1) = -l/2._prec
  do i = 2,imax
    xnew(i) = xnew(i-1)+dx(i-1)*l
  enddo

  

  
  !print*,
  !print*, 'x'  
 ! write(*,'(1e23.14)') xnew
  
  !print*,
  !print*, 'dx'
  !write(*,'(1e23.14)') dx
  

  
  xout = xnew




end subroutine calc_initial_grid


end module calc_grid
