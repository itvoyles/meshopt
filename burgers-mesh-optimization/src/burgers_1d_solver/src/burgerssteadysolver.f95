!******************************************************************************
!******************************************************************************
!**************** BURGERS EQN IMPLICIT SOLVER : SIMPLIFIED ********************
!******************************************************************************
!******************************************************************************
! Original Developer: Chris Roy
! Modified 8/6/2009 : Ed Alyanak 

!*******************************************************************************
!*************************** MODULES *******************************************
!*******************************************************************************


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
MODULE select_precision

IMPLICIT NONE

SAVE

INTEGER, PARAMETER :: DOUBLE = Selected_Real_Kind(p=13,r=200)
INTEGER, PARAMETER :: Prec = DOUBLE  

END MODULE select_precision


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
MODULE set_inputs

USE select_precision

IMPLICIT NONE

SAVE

INTEGER, PARAMETER 		::      imax = 17	! Nr of points along X-direction
INTEGER, PARAMETER 		::    	jmax = 1			! Number of points in the y-direction
INTEGER   						::      i                   ! Looping index for x-direction
INTEGER   						::      j = 1               ! Looping index for y-direction, but fixed in the 1D case
INTEGER   						::      iter                ! Looping index for iterative loop
INTEGER   						::  	itermax = 1000000    ! Maximum number of iterations
REAL(kind=Prec)					::      alpha = 32.0 !8.0_Prec !32.0_Prec    ! Scaling Factor
REAL(kind=Prec)					::    	dt 					! Time step
REAL(kind=Prec) 				::      tol = 1e-10_Prec 	! Tolerance for convergence
REAL(kind=Prec) 				::      RE = 32
REAL(kind=Prec) 				::      Lref = 8
REAL(kind=Prec) 				::      nu
REAL(kind=Prec) 				::      cfl = 1._Prec		! CFL number     
REAL(kind=Prec) 				::		t=0.0_Prec			! time
REAL(kind=Prec) 				::      theta = 1.0_Prec   	! Theta value determines the solver type
REAL(kind=Prec) 				::      alfa     			! 1st Coeffficient in the discretized equation
REAL(kind=Prec) 				::      beta 			   	! 2nd Coefficient in the discretized equation
REAL(kind=Prec) 				::      gamma			   	! 3Rd Coefficient in the discretized equation
REAL(KIND=Prec),DIMENSION(imax) :: 		uexact				! Exact Solution
REAL(kind=Prec) 				::      temp1               ! temp value to read
INTEGER							::		istatus  			! status of file being read
REAL(KIND=Prec),DIMENSION(imax) :: 		Jacob  				! Jacobian, a geometric property of the mesh
REAL(KIND=Prec),DIMENSION(imax) :: 		dJdxsi				! Derivative of the Jacobian
REAL(KIND=Prec),DIMENSION(imax) :: 		dt_vec				! time vector


END MODULE set_inputs

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! GEOMETRY
MODULE geometry

USE select_precision
USE set_inputs

IMPLICIT NONE

SAVE

REAL(kind=Prec)					::      L = 4.0_Prec		! Half Length of the domain
REAL(kind=Prec)					::      dx 				 	! Cell Size
REAL(kind=Prec)					::      dxsi=1.0_Prec	 	! Cell Size in the Computational domain
REAL(kind=Prec) 				::      xbegin = -4.0_Prec 	! First point in the domain
REAL(kind=Prec) 				::      xend = 4.0_Prec  	! Final point in the domain
REAL(KIND=Prec),DIMENSION(imax,jmax) ::	 	x   		 	! x - coordinates
REAL(KIND=Prec),DIMENSION(imax,jmax) ::	 	y   		 	! y - coordinates
REAL(KIND=Prec),DIMENSION(imax) 	 ::	 	xdash		 	! for calculating global time step
!REAL(KIND=Prec),DIMENSION(imax,jmax) ::	 	xold			!x - coordinates to check adaption criterion

END MODULE geometry
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! RESIDUAL CALC INPUT
MODULE residual_calc_input

USE select_precision

IMPLICIT NONE

SAVE

REAL(kind=Prec) ::      Res									! Variable used to calculate residue
REAL(kind=Prec) ::      Residue								! Residual
REAL(kind=Prec) ::      L2norm								! L2norm
REAL(kind=Prec) ::      R1									! Residue at the first iteration for normalization
REAL(KIND=Prec),DIMENSION(257) :: check2					! First Diagonal in the Tridiagonal system

END MODULE residual_calc_input
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! DISC ERROR INPUT
MODULE disc_error_input

USE select_precision
USE geometry

IMPLICIT NONE

SAVE

REAL(kind=Prec),DIMENSION(imax) :: Disc_error				! discretization error
REAL(kind=Prec) ::      disc_error_L1						! L1norm
REAL(kind=Prec) ::      disc_error_L2						! L2norm 

END MODULE disc_error_input
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! ITER CALC INPUT
MODULE it_calc_input

USE select_precision
USE set_inputs

IMPLICIT NONE

SAVE

REAL(kind=Prec),DIMENSION(IMAX) :: finalsoln				! Converged Velocity
REAL(kind=Prec),DIMENSION(IMAX) :: workingline				! Working value of velocity at a particular iteration
REAL(kind=Prec),DIMENSION(IMAX) :: iterror					! Iteration Error
REAL(kind=Prec)					:: temp_it					! Temporary Variable to store data		
REAL(kind=Prec) 				:: res_iterror_1			! Summation term for the L1 norm
REAL(kind=Prec) 				:: res_iterror_2			! Summation term for the L2 norm
REAL(kind=Prec) 				:: iterror_L2				! L2 norm of iterative errors
REAL(kind=Prec) 				:: iterror_L1				! L1 norm of the iterative errors
REAL(kind=Prec) 				:: iterror_Linf				! L1 norm of the iterative errors
INTEGER                         :: nvalue					! Integer to determine the location of the variable wrt to total 
INTEGER                         :: arraypos					! Position of the solution array from the solution history file
INTEGER							:: k						!	

END MODULE it_calc_input
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! SOLUTION
MODULE solution

USE select_precision
USE geometry
USE set_inputs

IMPLICIT NONE

SAVE

REAL(KIND=Prec),DIMENSION(imax) :: uold						! Original Velocity
REAL(KIND=Prec),DIMENSION(imax) :: unew						! New Velocity
REAL(KIND=Prec),DIMENSION(imax) :: y1						! Numerical Solution
REAL(KIND=Prec),DIMENSION(imax) :: ub						! Intermediate Velocity while solving

END MODULE solution
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
! TRIDIAG
MODULE tridiag

USE select_precision
USE set_inputs

IMPLICIT NONE

SAVE

REAL(KIND=Prec),DIMENSION(imax) :: AA						! First Diagonal in the Tridiagonal system
REAL(KIND=Prec),DIMENSION(imax) :: BB						! Second Diagonal in the Tridiagonal system  
REAL(KIND=Prec),DIMENSION(imax) :: CC						! Third Diagonal in the Tridiagonal system  
REAL(KIND=Prec),DIMENSION(imax) :: G						! Constant matrix in the RHS of MX=G
REAL(KIND=Prec),DIMENSION(3) :: checking					! First Diagonal in the Tridiagonal system

END MODULE tridiag


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

!*********************************************************************** 
!********************** SUBROUTINES  ***********************************
!***********************************************************************

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE initial_conditions

USE select_precision 
USE set_inputs
USE geometry
USE solution

IMPLICIT NONE

INTEGER							:: k				 ! Dummy index
INTEGER							:: iskip			 ! Number of entries to skip in the grid.dat file
INTEGER							:: imax_grid		 ! Number of points in the grid.dat file

! READ GRID
! Grid Format
! x1
! x2 
! ...
! xN
!write(*,*) 'Reading Grid'
OPEN(unit=12,file='grid.grd',STATUS='OLD',ACTION='READ',IOSTAT=istatus)
read(12,*) !skip line
read(12,*) !skip line
read(12,*) (x(i,j),i=1,imax)

CLOSE(12)

if (x(imax,j) /= 4.0_Prec) then
	write(*,*) 'GRID ERROR X(imax,j) NOT EQUAL TO 4.0'
endif

xbegin	= x(1,j)
xend	= x(imax,j)

! Initial Estimate and Exact Solution
DO i = 1,imax
	uold(i)		= (-2.0_Prec*nu*alpha/L)*dsinh(x(i,j)*alpha/L)/(dcosh(x(i,j)*alpha/L)+dexp(-t*nu*alpha**2/L**2))  
    uexact(i)	= -2.0_Prec*nu*alpha/L*dtanh(x(i,j)*alpha/L)   ! Exact Solution      
	y(i,j)		= 0.0_Prec                       
END DO 

! Set up metric output file
!open(50,file='gridinfo.dat',status='unknown')
!write(50,*)'TITLE = "BE Grid Transform. Data"'
!write(50,*)'variables="x(m)""y(m)""Jacobian""djdxsi""dxdxsi"'

CALL metrics

END SUBROUTINE initial_conditions



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE boundary_conditions

USE select_precision
USE set_inputs
USE geometry
USE solution

IMPLICIT NONE

! Setting up of the Boundary Conditions
uold(1) = -2.0_Prec*alpha*Nu/L*dtanh(xbegin*alpha/L)
uold(imax)=-2.0_Prec*alpha*Nu/L*dtanh(xend*alpha/L)

END SUBROUTINE boundary_conditions



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE metrics
 
 USE select_precision 
 USE set_inputs
 USE geometry
  
 IMPLICIT NONE
 
 !**************************************
 !********** Metric Calculations *******
 !**************************************
 
 ! Metrics being calculated are: 
 !       - J : Jacobian
 !       - dJdxsi: derivative of the Jacobian wrt xsi
 !       - dxdxsi: derivative of x wrt xsi
 
 ! Calculation of the Jacobian Array
 DO i=2,imax-1
	Jacob(i)=2._Prec/(x(i+1,j)-x(i-1,j))
 END DO    
 Jacob(1)=2._Prec/(-3._Prec*x(1,j)+4._Prec*x(2,j)-x(3,j))                  !2nd Order
 Jacob(imax)=2._Prec/(3._Prec*x(imax,j)-4._Prec*x(imax-1,j)+x(imax-2,j))   !
 
 ! Calculation of the derivative of the jacobian
 DO i=2,imax-1
	dJdxsi(i)=-1*((Jacob(i))**2)*(x(i+1,j)-2._Prec*x(i,j)+x(i-1,j))    ! Chain rule
 END DO   
 dJdxsi(1)=0._Prec
 dJdxsi(imax)=0._Prec
 
 ! Calculating xdash for the calculation of global timestep
 DO i=1,imax
 	xdash(i)=1._Prec/Jacob(i)
 END DO    

 !write(50,*) 'zone T="n = 0 " '
 !write(50,*) 'I=',imax,' J=',jmax
 !write(50,*) 'DATAPACKING=POINT'
 !write(50,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE )'
 !do i = 1, imax
 !  write(50,*)x(i,j),y(i,j),Jacob(i),dJdxsi(i),xdash(i)
 !enddo
 
END SUBROUTINE metrics



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tridiag_setup

USE select_precision
USE set_inputs
USE geometry
USE tridiag
USE solution

IMPLICIT NONE

DO i=3,imax-2
! Original Code
!          AA(i-1)=-1.0_Prec*theta*(alfa*ub(i)+beta)
!          BB(i-1)=(1.0_Prec+2.0_Prec*theta*beta)
!          CC(i-1)=theta*(alfa*ub(i)-beta)

!	Gen.Coord Code with theta
          AA(i-1) = theta*( -alfa*Jacob(i)*ub(i) - beta*(Jacob(i)**2) + gamma*dJdxsi(i)*Jacob(i))
          BB(i-1) = ( 1.0_Prec + 2.0_Prec*theta*beta*(Jacob(i)**2) )
          CC(i-1) = theta*alfa*Jacob(i)*ub(i) - theta*beta*(Jacob(i)**2) - theta*gamma*dJdxsi(i)*Jacob(i)
          G(i-1)  = uold(i)-((1.0_Prec-theta)*alfa*ub(i)*(uold(i+1)-uold(i-1))) + &
         & ((1.0_Prec-theta)*beta*(uold(i+1) - 2.0_Prec*uold(i) + uold(i-1)) )

! Gen. Coord Code w/o theta
!    	   	AA(i-1)=(-alfa*Jacob(i)*ub(i)-beta*(Jacob(i)**2)+gamma*dJdxsi(i))
!			BB(i-1)=(1.0_Prec+2.0_Prec*beta*(Jacob(i)**2))
!	        CC(i-1)= alfa*Jacob(i)*ub(i)-beta*(Jacob(i)**2)-gamma*dJdxsi(i)
!           G(i-1)=uold(i)
END DO 

! Setting up the end points of the diagonals 

! Original Burger's Eqn Tridiagonal Endpoints
!       AA(1)=0._prec;
!       AA(imax-2)=-1.0_Prec*(theta*alfa*ub(imax-1)+theta*beta)
!        BB(1)=(1.0_Prec+2.0_Prec*beta*J(2)**2)
!       BB(imax-2)=(1.0_Prec+2.0_Prec*theta*beta)
!       CC(1)= (theta*alfa*ub(2)-theta*beta)
!       CC(imax-2)=0._prec;

! Gen Coord Code Tridiagonal Endpoints
       AA(1)=0._prec;
       AA(imax-2)=( -theta*alfa*Jacob(imax-1)*ub(imax-1) - theta*beta*Jacob(imax-1)**2 & 
       				+ gamma*dJdxsi(imax-1)*Jacob(imax-1))
       BB(1)=(1.0_Prec + 2.0_Prec*theta*beta*Jacob(2)**2)    
       BB(imax-2)=( 1.0_Prec + 2.0_Prec*theta*beta*Jacob(imax-1)**2)
       CC(1)= (theta*alfa*Jacob(2)*ub(2) - theta*beta*Jacob(2)**2 - gamma*dJdxsi(2)*Jacob(2))
       CC(imax-2)=0._prec;

       G(1)=uold(2) - ((1.0-theta)*alfa*ub(2)*(uold(3) - uold(1)))+ &
       & ((1.0 - theta)*beta*(uold(3) - 2.0_Prec*uold(2) + uold(1)))+ &
       & (theta*alfa*ub(2)*Jacob(2) + theta*beta*(Jacob(2))**2 - theta*gamma*dJdxsi(2)*Jacob(2))*uold(1)                  
       G(imax-2)=uold(imax-1) - ((1.0-theta)*alfa*ub(imax-1)*(uold(imax) - uold(imax-2)))&
       & +((1.0 - theta)*beta*(uold(imax) - 2.0_Prec*uold(imax-1) + uold(imax-2)))&
       & - (theta*alfa*Jacob(imax-1)*ub(imax-1) - theta*beta*(Jacob(imax-1))**2 & 
       - theta*gamma*dJdxsi(imax-1)*Jacob(imax-1))*uold(imax)

       checking(1)=((1.0-theta)*alfa*ub(2)*(uold(3)-uold(1)))
       checking(2)=((1.0-theta)*beta*(uold(3)-2.0_Prec*uold(2)+uold(1)))
       checking(3)=-(theta*alfa*Jacob(imax-1)*ub(imax-1)-theta*beta*(Jacob(imax-1))**2)*uold(imax)

END SUBROUTINE tridiag_setup



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE tridiag_solver (n,AA,B,CC,X,G)

USE set_inputs

IMPLICIT NONE

INTEGER,INTENT(IN) 						:: n			! Number of equations to be solved
INTEGER 								:: k    		! Variable to determine the position of X
REAL(KIND=Prec)							:: Tx			! Row Multiplier in Gaussian Elmnsn
REAL(KIND=Prec),DIMENSION(N),INTENT(IN) :: AA,B,CC,G   	! Diagonals
REAL(KIND=Prec),DIMENSION(N),INTENT(OUT):: X 			! Solution 
REAL(KIND=Prec),DIMENSION(N)			:: BB,GG		! Working Arrays

!.....THIS SUBROUTINE SOLVES TRIDIAGONAL SYSTEMS OF EQUATIONS
!.....BY GAUSS ELIMINATION
!.....THE PROBLEM SOLVED IS MX=G WHERE M=TRI(A,B,C)
!.....THIS ROUTINE DOES NOT DESTROY THE ORIGINAL MATRIX
!.....AND MAY BE CALLED A NUMBER OF TIMES WITHOUT REDEFINING
!.....THE MATRIX
!.....N = NUMBER OF EQUATIONS SOLVED (UP TO 1000)
!.....FORWARD ELIMINATION
!.....BB IS A SCRATCH ARRAY NEEDED TO AVOID DESTROYING B ARRAY

!Transfering values to the working arrays
BB = B
GG=G

DO I=2,n
 Tx = AA(I)/BB(I-1)
 BB(I) = BB(I) - CC(I-1)*Tx
 GG(I) = GG(I) - GG(I-1)*Tx
END DO

!.....BACK SUBSTITUTION
X(n) = GG(n)/BB(n)
DO  I=1,n-1
 k = n-I
 X(k) = (GG(k)-CC(k)*X(k+1))/BB(k)
END DO

RETURN
END SUBROUTINE tridiag_solver



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE residual_calc

USE select_precision
USE set_inputs
USE geometry
USE solution
USE residual_calc_input

IMPLICIT NONE

Residue = 0._Prec
DO i=2,imax-1
	 Res = unew(i)*Jacob(i)/(2.0_Prec)*(unew(i+1)-unew(i-1))-Nu*(Jacob(i)**2)* (unew(i+1)-2.0_Prec*unew(i)+unew(i-1))&
     & -nu*dJdxsi(i)*Jacob(i)*(unew(i+1)-unew(i-1))/2.0_Prec
	 Residue = Residue + Res**2
END DO

L2Norm =sqrt(Residue/(imax-2))
IF (iter==1) r1 = L2Norm
L2Norm=L2Norm/r1

if(mod(iter,10)==0) then
  ! Printing Output on screen
  WRITE(*,11) iter,dt*REAL(iter),dt,L2Norm
  11 FORMAT ( 1X, I6, e24.16,e24.16,e24.16)
  
  ! Printing Output onto a file
  WRITE (57,134) iter, L2norm

endif
  134 FORMAT (I9,e24.16)
END SUBROUTINE residual_calc


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
SUBROUTINE discretization_error

USE set_inputs
USE geometry
USE disc_error_input
USE solution

IMPLICIT NONE

	 DO i = 1,imax
       uexact(i)=-2.0_Prec*Nu*alpha/L*dtanh(x(i,j)*alpha/L)   ! Exact Solution                              
     END DO 

     Disc_error=unew-uexact
     
     disc_error_L1= 0.0_Prec
     disc_error_L2= 0.0_Prec
     
	 DO i=1,imax
       disc_error_L1 = disc_error_L1 + abs(Disc_error(i))
       disc_error_L2 = disc_error_L2 + (abs(Disc_error(i)))**2     
     END DO

     disc_error_L1 = disc_error_L1/REAL(imax)
     disc_error_L2 = dsqrt(disc_error_L2/REAL(imax))

!     WRITE(*,81) disc_error_L1
!81    FORMAT(1X,'L1Norm of discretisation error=',E24.16)      
!     WRITE(*,8) disc_error_L2
!8     FORMAT(1X,'L2Norm of discretisation error=',e24.16)

!     WRITE(30,83) disc_error_L1,disc_error_L2
!83    FORMAT (1X,'variables="L1 norm ="',e24.16,/,'"L2norm =',e24.16)
     
     !DO I=1,IMAX
     !  WRITE(30,334)x(i), unew(I),uexact(i),unew(i)-uexact(i)  
     !  WRITE(30,334) (unew(I),I=1,imax)
     !334 FORMAT (1X, 4E14.5)
     !END DO

END SUBROUTINE discretization_error


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************

!*********************************************************************** 
!********************** MAIN PROGRAM ***********************************
!***********************************************************************
 
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
PROGRAM steady_implicit_original_burgers_eqn
 
 USE select_precision
 USE set_inputs
 USE geometry
 USE solution
 USE tridiag
 USE residual_calc_input
 USE disc_error_input
 USE it_calc_input

IMPLICIT NONE 


alpha = RE/2._prec
nu = 2._prec*Lref/RE


CALL initial_conditions  
CALL boundary_conditions ! set based on exact solution
 
! Set Input values
!dx = L*2._Prec/REAL(imax-1) ! Cell Distance

!xxxxxxxxxxxxxxxxxxxxxxxxxxx Main loop start xxxxxxxxxxxxxxxxxxxxxx 
iteration:DO iter=1,itermax
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
			
           ub = uold   ! interchanging velocity values
           
		   ! determine local time step for given CFL value
           DO i=1,imax
             dt_vec(i) = cfl*(xdash(i)**2)/(abs(uold(i))*xdash(i)+2._Prec*nu)
           END DO
		   ! select minimum local time step for calculation  
           dt=MINVAL(dt_vec)
        
           alfa = dt/2._Prec !*dxsi!*Jacob(i)          ! 1st Coefficient in the discretized equation
           beta = nu*dt  !/dxsi**2!*Jacob(i)**2        ! 2nd Coefficient in the discretized equation
           gamma=nu*dt/2.0_Prec !*dJdxsi(i)
            
           ! Setting up the tridiagonal Matrix
           CALL tridiag_setup
        
           !Calling Tridiagonal Solver
           CALL tridiag_solver (imax-2,AA,BB,CC,y1,G)
 
           ! Total Velocity Vector
           DO i=1,imax-2
             unew(i+1)=y1(i)
           END DO
           unew(1)=uold(1)
           unew(imax)=uold(imax)
 
           CALL residual_calc 

           IF (L2Norm.LT.Tol) then
		   write(*,*) 'TOL MET, EXIT AT ITERATION ', iter
		   EXIT  
		   END IF
           
		   uold = unew

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx    
END DO iteration
!xxxxxxxxxxxxxxxxxxxxxxxxxxx Main loop ends xxxxxxxxxxxxxxxxxxxxxxx

if (iter < itermax) then 
	WRITE(*,*) 'converged'
else
	WRITE(*,*) 'ITERMAX REACHED SOLUTION MAY NOT BE CONVERGED', L2Norm
	if (Maxval(abs(unew)) > 10000.0_Prec) then
		WRITE(*,*) 'SOLUTION DIVERGED, SETING SOLUTION TO -1'
		do i=1,imax
			unew(i) = -1.0_Prec
		enddo
	endif
endif
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
CALL discretization_error
 
!xxxxxxxxxxxxxxxxxxxxx Writing out Solution xxxxxxxxxxxxxxxxxxxxxx
        
OPEN(40,file='output.dat',status='unknown')
WRITE(40,*) 'VARIABLES="x(i,j)""unew(i)""uexact(i)""unew(i)-uexact(i)"'
DO i=1,imax
WRITE(40,43) x(i,j),unew(i),uexact(i),unew(i)-uexact(i)
43  FORMAT ( 1X, 4e24.16 )
END DO
CLOSE(40)
OPEN(50,file='finished.dat',status='unknown')
Write(50,*) '1'
close(50)
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
END PROGRAM steady_implicit_original_burgers_eqn

