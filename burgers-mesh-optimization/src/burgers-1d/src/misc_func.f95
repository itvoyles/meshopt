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
