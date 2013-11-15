module select_precision 

! use selected_real_kind to select desired kinds of real variables in 
! a processor-independent manner 

implicit none

save

! Declare Parameters:
integer, parameter :: sngl  = selected_real_kind(p=6,  r=37)
integer, parameter :: dbl   = selected_real_kind(p=13, r=200)
integer, parameter :: extnd = selected_real_kind(p=17, r=2000)
integer, parameter :: quad  = selected_real_kind(p=26, r=200)
integer, parameter :: prec = dbl 

end module select_precision

module simple_spring
use select_precision, only : prec




contains 

subroutine springsystem(x,y,kin,imax,jmax,Lref)
  implicit none
  integer,intent(in) :: imax,jmax
  real(prec), intent(in) :: Lref
  real(prec), intent(in), dimension((imax-1)*jmax+imax*(jmax-1)) :: kin
  real(prec), dimension(imax,jmax), intent(out) :: x,y
  
  real(prec), dimension( (imax-1)*jmax ) :: k
  real(prec), dimension( (jmax-1)*imax ) :: j
  
  real(prec), dimension(imax,jmax) :: ytemp
  real(prec), dimension(imax*jmax) :: ytempvec
  real(prec), dimension(imax,jmax-1) :: jtemp
  real(prec), dimension(jmax) :: xl,xr
  real(prec), dimension(imax) :: yb,yt
  
  integer :: i,jj,cnt
  
  k = kin( 1:(imax-1)*jmax )
  xl = -Lref/2._prec
  xr =  Lref/2._prec
  x = springsystem_loc(k,imax,jmax,(imax-1)*jmax,xl,xr)
  
  j = kin( (imax-1)*jmax+1 : (imax-1)*jmax+imax*(jmax-1) )
!converts to a matrix and flips left/right
  cnt = 1
  do jj = 1,jmax-1
  do i = 1,imax
    jtemp(imax-i+1,jj) = j(cnt)
    cnt = cnt + 1
  enddo
  enddo
  

  
  cnt = 1
  do i = 1,imax
  do jj = 1,jmax-1
    j(cnt) = jtemp(i,jj)
    cnt = cnt + 1
    !print*, j(cnt)
  enddo
  enddo
  !stop
!calculates y in the rotated system
  yb = -Lref/2._prec
  yt =  Lref/2._prec
  y = springsystem_loc(j,jmax,imax,(jmax-1)*imax,yb,yt)

  y = transpose(y)
  cnt = 1
  do jj = 1,jmax
  do i = 1,imax
    ytemp(imax+1-i,jj) = y(i,jj)
  enddo
  enddo
  
 ! y = transpose(ytemp)

end subroutine springsystem




function springsystem_loc(k,imax,jmax,kn,xl,xr)
  implicit none
                
  integer, intent(in) :: imax, jmax,kn
 
  real(prec), intent(in), dimension(kn) :: k
  real(prec), dimension(imax), intent(in) :: xl,xr
  real(prec), dimension(imax,jmax) :: springsystem_loc
  real(prec), dimension(imax*jmax) :: x
                       
  real(prec), dimension(imax*jmax) :: LL,D,UU,B
  integer :: ii,jj,nk,nj,ecnt,wcnt,cnt

  
  nk = (imax-1)*jmax 
  ecnt = 1
  wcnt = 1
  cnt  = 1
  !setup primary matrix
  do jj = 1,jmax
    do ii = 1,imax
        
        D(cnt) = 0._prec
        if (ii<imax) then
            UU(cnt) = -k(ecnt)
            D(cnt) = k(ecnt)
            ecnt = ecnt + 1;
        else
            UU(cnt) = 0._prec
        endif
        
        if (ii>1) then
            LL(cnt) = -k(wcnt)
            D(cnt) = D(cnt) + k(wcnt)
            wcnt = wcnt + 1;
        else
            LL(cnt) = 0._prec
        endif
        
        
        
        cnt = cnt + 1;
        
    enddo
  enddo
  
  !setup boundary conditions
  B = 0._prec
  cnt = 1
  !left boundary
  do ii = 1,imax*jmax,imax
    UU(ii) = 0._prec
    LL(ii) = 0._prec
    D(ii) = 1._prec
    B(ii) = xl(cnt)
    cnt = cnt + 1
  enddo
  
  !right boundary
  cnt = 1
  do ii = imax,imax*jmax,imax
    UU(ii) = 0._prec
    LL(ii) = 0._prec
    D(ii) = 1._prec
    B(ii) = xr(cnt)
    cnt = cnt + 1
  enddo
  

  
  x = trisolve(LL,D,UU,B,imax*jmax)

  cnt = 1
  do jj = 1,jmax
  do ii = 1,imax
    springsystem_loc(ii,jj) = x(cnt)
    cnt = cnt + 1
  enddo
  enddo
  
end function springsystem_loc


!Tri-diagonal linear algebra solver --------------------------------------------
function trisolve(LLin,Din,UUin,Bin,m)

implicit none
  integer, intent(in) :: m
  real(prec), dimension(m), intent(in) :: LLin,Din,UUin,Bin
  real(prec), dimension(m) :: LL,D,UU,B
  real(prec), dimension(m) :: trisolve

  real(prec) :: k
  integer :: i,j
  
  LL = LLin
  D = Din
  UU = UUin
  B  = Bin
  
  do i= 1,m-1
     k = LL(i+1)/D(i)
     LL(i+1) = LL(i+1)-k*D(i)
     D(i+1) = D(i+1)-k*UU(i)
     B(i+1) = B(i+1)-k*B(i)   
  enddo

  trisolve = 0._prec
  trisolve(m) = B(m)/D(m)
  do j = 1,m-1
    i = m-j
    trisolve(i) = ( B(i) - UU(i)*trisolve(i+1) )/D(i) 
  enddo

end function trisolve


subroutine springsystemtest
  implicit none
  
  integer, parameter :: imax = 9
  integer, parameter :: jmax = 9
  integer, parameter :: nk = (imax-1)*jmax
  integer, parameter :: nj = imax*(jmax-1)
  real(prec), dimension(nk+nj) :: k
  
  real(prec), dimension(imax,jmax) :: x,y
  real(prec) :: Lref = 1._prec
  
  integer :: i,j,cnt
  real(prec) :: pi,piloc
  
  pi = acos(-1._prec)
  
  cnt = 1
  do j = 1,jmax
  do i = 1,imax-1
    piloc = pi*(real(i-1,prec)/real(imax-2,prec))/2._prec
    k(cnt) = 10._prec*sin(piloc)+1._prec
    cnt = cnt + 1
  enddo
  enddo
  
  cnt = nk+1
  do j = 1,jmax-1
  do i = 1,imax
    piloc = pi*(real(j-1,prec)/real(jmax-2,prec))/2._prec
    k(cnt) = 10._prec*sin(piloc)+1._prec
    cnt = cnt + 1
  enddo
  enddo

  !do i = 1,nk+nj
  !  write(*,'(I4,1g12.4)')i,k(i)
  !enddo
  
  
  call springsystem(x,y,k,imax,jmax,Lref)
  
  write(*,'(9e12.4)')x
  print*,
  write(*,'(9e12.4)')y
  

end subroutine



end module simple_spring

