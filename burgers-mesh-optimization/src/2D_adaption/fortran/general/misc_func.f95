module misc_func
!module contains misc. functions for various simple calculations
!
!FUNCTIONS:
! factoriali(var) => elemental function, integer factorial of var
! gaussian_elimination(x,A,B,m) => subroutine, solves system Ax=B of size m
!   function is broken up into two parts for optimization
!   gaussian_elimination_p1 => subroutine, reduces A to diag matrix
!   gaussian_elimination_p2 => subroutine, reduces b matrix from p1 history

contains
!*******************************************************************************
! integer factorial of var
!*******************************************************************************
  elemental function factoriali(var)
    implicit none
    integer, intent(in) :: var
    integer :: factoriali

    integer :: i
  
    factoriali=1
    do i = 2,var
      factoriali=factoriali*i
    enddo
    return
  end function
  
!integration using trapezoidal rule ********************************************
function trap_sum(x,datain,n)
  use select_precision, only : prec
  implicit none
  
  integer, intent(in)                 :: n
  real(prec), dimension(n),intent(in) :: x,datain
  real(prec)              :: trap_sum
  
  real(prec):: J
  integer   :: i
  
  J = 0._prec
  do i = 1,n-1
    J = (x(i+1)-x(i))*(datain(i)+datain(i+1))/2._prec + J
  enddo
  trap_sum = J

end function trap_sum

function trap_sum2dopt(x,y,f,n,m)
  use select_precision, only : prec
  implicit none

  integer, intent(in) :: n,m
  real(prec), dimension(n,m), intent(in) :: x,y,f
  real(prec)                             :: trap_sum2dopt


  real(prec),dimension(n-1,m-1) :: Area,fave
  real(prec),dimension(n-1,m-1) :: fold
  real(prec), dimension(n-1,m-1,2) :: v1,v2,v3,v4
  real(prec), dimension(n-1,m-1)   :: x1,x2,x3,x4,y1,y2,y3,y4
  integer, dimension(n-1,m-1) :: ii,jj
  integer :: i,j
  continue





  fave(:,:) = (f(1:n-1,1:m-1)+f(2:n,2:m)+f(1:n-1,2:m)+f(2:n,1:m-1))/4._prec




 ! if (sum(Area)==0._prec) then
  where (fave.ne.fold)

  x1(:,:) = x(1:n-1,1:m-1)
  x2(:,:) = x(2:n,1:m-1)
  x3(:,:) = x(1:n-1,2:m)
  x4(:,:) = x(2:n,2:m)

  y1(:,:) = y(1:n-1,1:m-1)
  y2(:,:) = y(2:n,1:m-1)
  y3(:,:) = y(1:n-1,2:m)
  y4(:,:) = y(2:n,2:m)

  !calculates area by calculating the area of two half-quads 
  !A = (v1 x v2)/2+(v3 x v4)/2  [v1 x v2 => crossproduct]

  v1(:,:,1) = x2-x1
  v1(:,:,2) = y2-y1

  v2(:,:,1) = x3-x1
  v2(:,:,2) = y3-y1

  v3(:,:,1) = x3-x4
  v3(:,:,2) = y3-y4

  v4(:,:,1) = x2-x4
  v4(:,:,2) = y2-y4

  Area = 0.5_prec*abs(v1(:,:,1)*v2(:,:,2)-v1(:,:,2)*v2(:,:,1)) &
       + 0.5_prec*abs(v3(:,:,1)*v4(:,:,2)-v3(:,:,2)*v4(:,:,1))

 ! endif

  end where
  fold = fave
  trap_sum2dopt = sum(Area*fave)


end function


function trap_sum2d(x,y,f,n,m)
  use select_precision, only : prec
  implicit none

  integer, intent(in) :: n,m
  real(prec), dimension(n,m), intent(in) :: x,y,f
  real(prec)                             :: trap_sum2d


  real(prec),dimension(n-1,m-1) :: Area,fave
  real(prec), dimension(n-1,m-1,2) :: v1,v2,v3,v4



  v1(:,:,1) = x(2:n,1:m-1)-x(1:n-1,1:m-1)
  v1(:,:,2) = y(2:n,1:m-1)-y(1:n-1,1:m-1)

  v2(:,:,1) = x(1:n-1,2:m)-x(1:n-1,1:m-1)
  v2(:,:,2) = y(1:n-1,2:m)-y(1:n-1,1:m-1)

  v3(:,:,1) = x(1:n-1,2:m)-x(2:n,2:m)
  v3(:,:,2) = y(1:n-1,2:m)-y(2:n,2:m)

  v4(:,:,1) = x(2:n,1:m-1) - x(2:n,2:m)
  v4(:,:,2) = y(2:n,1:m-1) - y(2:n,2:m)

  !calculates area by calculating the area of two half-quads 
  !A = (v1 x v2)/2+(v3 x v4)/2  [v1 x v2 => crossproduct]

  Area = 0.5_prec*abs(v1(:,:,1)*v2(:,:,2)-v1(:,:,2)*v2(:,:,1)) &
       + 0.5_prec*abs(v3(:,:,1)*v4(:,:,2)-v3(:,:,2)*v4(:,:,1))

  fave = (f(1:n-1,1:m-1)+f(2:n,2:m)+f(1:n-1,2:m)+f(2:n,1:m-1))/4._prec
  trap_sum2d= sum(Area*fave)


end function


!*******************************************************************************
! Matrix inverse calculation
!*******************************************************************************
function inv(A,n)
  use select_precision, only : prec
  implicit none
  
  integer, intent(in) :: n
  real(prec), intent(in), dimension(n,n) :: A
  real(prec), dimension(n,n) :: inv
  real(prec), dimension(n,n) :: test
  real(prec), dimension(n,2*n) :: Aug
  integer :: row,col,i
  real(prec) :: k
  
  Aug = 0._prec
  Aug(:,1:n) = A
  do col = 1,n
    Aug(col,col+n) = 1._prec
  enddo
  



  ! use Gaussian elimination on Augmented matrix to reduce left matrix to diagonal
  ! and resulting in a fully populated right matrix
  do col = 1,n
    do row = 1,n

      if (row/=col) then
        k = Aug(row,col)/Aug(col,col)
        Aug(row,:) = Aug(row,:)-k*Aug(col,:)
      endif

    enddo
  enddo

  
  ! divide by diagonal to reduce left matrix to identity matrix and finalize
  ! inverse calulations
  do row = 1,n
    Aug(row,:)=Aug(row,:)/Aug(row,row)
  enddo
  
 
  inv = Aug(:,n+1:2*n)


end function inv



!*******************************************************************************
! gaussian elimination to solve linear system Ax=B of size m
!*******************************************************************************
  subroutine gaussian_elimination(x,A,B,m)
  ! Solves linear system Ax=B
    use select_precision  
    implicit none
    
    Integer, intent(in) :: m !m x m size matrix to solve
    real(prec), dimension(m,m), intent(in) :: A
    real(prec), dimension(m), intent(in)   :: B
    real(prec), dimension(m), intent(out)  :: x
    
    integer :: col,row,j
    real(prec), dimension(m,m) :: Aold, Anew
    real(prec), dimension(m,m) :: kk
    real(prec), dimension(m) :: normalize
    real(prec) :: k
    
    call gaussian_elimination_p1(kk,normalize,A,m)
    call gaussian_elimination_p2(x,kk,normalize,B,m)

	
  end subroutine gaussian_elimination

!*******************************************************************************

!*******************************************************************************
  subroutine gaussian_elimination_p1(history,diag,A,m)
    use select_precision  
    implicit none
    
    Integer, intent(in) :: m !m x m size matrix to solve
    real(prec), dimension(m,m), intent(in) :: A
    real(prec), dimension(m,m), intent(out):: history
    real(prec), dimension(m), intent(out)  :: diag
    
    integer :: col,row,j
    real(prec), dimension(m,m) :: Aold, Anew
    real(prec) :: k
    

    !diagonalize A matrix and record reduction history
    history = 0._prec
    Aold = A
    Do col = 1,m

			Do row = 1,m
				k = Aold(row,col)/Aold(col,col)
				
				if (row.eq.col) then
				  k = 0._prec
				endif
			
				Anew(row,1:m) = Aold(row,1:m)-k*Aold(col,1:m)
				history(row,col) = k

			EndDo	
				
			
			Aold = Anew
			diag(col) = Aold(col,col)			
	  EndDo


	   
  end subroutine gaussian_elimination_p1
  
  
!*******************************************************************************

!*******************************************************************************
  subroutine gaussian_elimination_p2(x,history,diag,B,m)
    use select_precision, only : prec
    implicit none
    
    Integer, intent(in) :: m !m x m size matrix to solve
    real(prec), dimension(m,m), intent(in):: history
    real(prec), dimension(m), intent(in)  :: diag,B
    real(prec), dimension(m), intent(out) :: x
    
    integer :: col
    
  	!simplify B matrix using A matrix reduction history and diagonal
	  x = B
	  Do col = 1,m
			x = x-history(1:m,col)*x(col)
		enddo
		x = x/diag
	 
  
  end subroutine gaussian_elimination_p2
!*******************************************************************************
!*******************************************************************************
!                          END GUASSIAN ELIMINATION
!*******************************************************************************
!*******************************************************************************

end module misc_func
