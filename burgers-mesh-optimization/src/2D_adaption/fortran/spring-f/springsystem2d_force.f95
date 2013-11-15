module springsystem2d_force
  use select_precision, only : prec
  
  real(prec), allocatable, dimension(:,:) :: Kx, Ky, Kxinv, Kyinv
  
contains
  subroutine springsystem_setup(imax,jmax)
  use misc_func, only : inv
  
  implicit none
  integer, intent(in) :: imax, jmax
  integer :: cnt,i,j
  
  real(prec), dimension(imax*jmax,imax*jmax) :: K,testx, testy

  
  allocate(Kx(imax*jmax,imax*jmax),    &
           Kxinv(imax*jmax,imax*jmax), &
           Ky(imax*jmax,imax*jmax),    &
           Kyinv(imax*jmax,imax*jmax))
          
          
          
  ! Setup coefficient matrix
  K = 0._prec
  cnt = 1 
  do j = 1,jmax
    do i = 1,imax
      if (i>1)  K(cnt,cnt-1) = -1._prec
      if (j>1)  K(cnt,cnt-imax) = -1._prec
      if (i<imax) K(cnt,cnt+1) = -1._prec
      if (j<jmax) K(cnt,cnt+imax) = -1._prec
      
      K(cnt,cnt) = -sum(K(cnt,:))
      cnt = cnt + 1
      
    enddo
  enddo


          
  ! Setup boundary conditions for each direction
  Kx = K
  Ky = K

  ! x-bc
  do i = 1,imax*jmax,imax
    Kx(i,:) = 0._prec
    Kx(i,i) = 1._prec
  enddo

  do i = imax,imax*jmax,imax
    Kx(i,:) = 0._prec
    Kx(i,i) = 1._prec
  enddo

  ! y-bc
  do i = 1,imax
    Ky(i,:) = 0._prec
    Ky(i,i) = 1._prec
  enddo

  do i = imax*(jmax-1)+1,imax*jmax
    Ky(i,:) = 0._prec
    Ky(i,i) = 1._prec
  enddo
  

  Kxinv = inv(Kx,imax*jmax)
  Kyinv = inv(Ky,imax*jmax)

   
  
  end subroutine springsystem_setup  
  
  !*****************************************************************************
  function forwardspringsystem(F,imax,jmax)
  ! inputs a vector of forces with size(2*imax*jmax,1)
  ! returns an array of nodal locations of size(imax,jmax,2)
  use select_precision, only : prec
  implicit none
  
  integer, intent(in) :: imax,jmax
  real(prec),dimension(imax*jmax), intent(in) :: F
  real(prec), dimension(imax,jmax,2) :: forwardspringsystem
  real(prec), dimension(imax*jmax) :: xvec, yvec
  integer :: i,j
  
  xvec = matmul(Kxinv,F(1:imax*jmax))
  yvec = matmul(Kyinv,F(imax*jmax+1:2*imax*jmax))
  
  forwardspringsystem(:,:,1) = reshape(xvec,(/imax,jmax/))
  forwardspringsystem(:,:,2) = reshape(yvec,(/imax,jmax/))  
  
  end function forwardspringsystem
  
  !*****************************************************************************
  function reversespringsystem(x,y,imax,jmax)
    !input an array of x and y nodal locations and sizes imax and jmax
    !returns a vector of forces with size(2*imax*jmax,1)
    use select_precision, only : prec
    implicit none
    
    integer, intent(in) :: imax,jmax
    real(prec), dimension(imax,jmax),intent(in) :: x,y
    real(prec), dimension(imax*jmax*2) :: reversespringsystem
    
    reversespringsystem(1:imax*jmax) = matmul(Kx,reshape(x,(/imax*jmax/)))
    reversespringsystem(imax*jmax+1:2*imax*jmax) = matmul(Ky,reshape(y,(/imax*jmax/)))
  
  
  end function reversespringsystem
  
  
  
end module springsystem2d_force
