!******************************************************************************
!***************************    MODULES   *************************************
!******************************************************************************

Module Select_Precision 

! Use Selected_Real_Kind to select desired kinds of real variables in a processor-independent manner 

Implicit None

Save

! Declare Parameters:
Integer, Parameter :: Single = Selected_Real_Kind(p=6,r=37)
Integer, Parameter :: Double = Selected_Real_Kind(p=13,r=200)
Integer, Parameter :: Quad   = Selected_Real_Kind(p=26,r=200)
Integer, Parameter :: Prec = Double     ! Precision for the computations is set on this line
                                        ! Prec = Single gives single precision
                                        ! Prec = Double gives double precision
End Module

!******************************************************************************
Module Set_Inputs

Use Select_Precision

! Set the code inputs for the 2D heat conduction code

Implicit None

Save

Integer, Parameter ::   imax = 9 !65			 ! Number of points in the x-direction
Integer, Parameter ::   jmax = 9 !81			 ! Number of points in the y-direction
Integer ::              nmax = 100000      	 ! max # of steps
Integer ::              iterout = 50000      ! number of time steps between output
Integer ::              imms = 1             ! Exact sol. flag ( = 1 for exact sol., = 0 otherwise) 
Integer ::              idebug = 0           ! debug flag ( = 1 to write debug info, = 0 otherwise) 
Real(kind=Prec) ::      cfl  = 0.1_Prec !225_Prec    ! maximum CFL for each i,j location
Real(kind=Prec) ::      toler = 1.0e-10_Prec  ! tolerance for convergence
Real(kind=Prec) ::      RE = 8
Real(kind=Prec) ::      Uref = 2
Real(kind=Prec) ::      alpha = 16.00_Prec    ! scaling parameter
Real(kind=Prec) ::      rnu = 0.5_Prec   	 ! viscosity (m^2/s)
Real(kind=Prec) ::      fsmall = 1.e-12_Prec ! small parameter
Real(kind=Prec) ::      dt_max = 1.e99_Prec  ! maximum delta t
Real(kind=Prec) ::      rtime = 0.0_Prec     ! set initial time to zero
Real(kind=Prec) ::      Pi 

End Module

!******************************************************************************
Module Geometry

Use Select_Precision
Use Set_Inputs

! Set the geometry for the 2D heat conduction code

Implicit None

Save

! Declare Parameters:
Integer ::              i               ! Looping index for xsi-direction
Integer ::              j               ! Looping index for eta-direction
Real(kind=Prec) ::      xmin = -4.0_Prec     ! Domain dimensions (Cartesian only)
Real(kind=Prec) ::      xmax = 4.0_Prec
Real(kind=Prec) ::      ymin = -4.0_Prec
Real(kind=Prec) ::      ymax = 4.0_Prec
Real(kind=Prec) ::      rLength=8._prec 	 ! Characteristic length
Real(kind=Prec) ::      dx		 ! Node spacing in x-direction
Real(kind=Prec) ::      dy 		 ! Node spacing in y-direction
Real(kind=Prec),Dimension(imax,jmax) ::  x = -99.9 ! x-coordinates for grid
Real(kind=Prec),Dimension(imax,jmax) ::  y = -99.9 ! y-coordinates for grid
Real(kind=Prec),Dimension(imax,jmax) ::  vol = -99.9 ! volume of cell surrounding node
Real(kind=Prec),Dimension(imax,jmax) ::  xsix = -99.9 ! xsi_x/J grid metric
Real(kind=Prec),Dimension(imax,jmax) ::  xsiy = -99.9 ! xsi_y/J grid metric
Real(kind=Prec),Dimension(imax,jmax) ::  etax = -99.9 ! eta_x/J grid metric
Real(kind=Prec),Dimension(imax,jmax) ::  etay = -99.9 ! eta_y/J grid metric

Real(kind=Prec),Dimension(imax,jmax) ::  d_xsix_dxsi = -99.9 ! grid metric derivative
Real(kind=Prec),Dimension(imax,jmax) ::  d_xsiy_dxsi = -99.9 ! grid metric derivative
Real(kind=Prec),Dimension(imax,jmax) ::  d_etax_dxsi = -99.9 ! grid metric derivative
Real(kind=Prec),Dimension(imax,jmax) ::  d_etay_dxsi = -99.9 ! grid metric derivative
Real(kind=Prec),Dimension(imax,jmax) ::  d_xsix_deta = -99.9 ! grid metric derivative
Real(kind=Prec),Dimension(imax,jmax) ::  d_xsiy_deta = -99.9 ! grid metric derivative
Real(kind=Prec),Dimension(imax,jmax) ::  d_etax_deta = -99.9 ! grid metric derivative
Real(kind=Prec),Dimension(imax,jmax) ::  d_etay_deta = -99.9 ! grid metric derivative

End Module

!******************************************************************************
Module MMS_Constants

Use Select_Precision

Implicit None

Save
! Note: these are currently not used
Real(kind=Prec) :: T0 = 400.0_Prec
Real(kind=Prec) :: Tx = 45.0_Prec
Real(kind=Prec) :: Ty = 35.0_Prec
Real(kind=Prec) :: Txy = 27.5_Prec
Real(kind=Prec) :: ax = 1.0_Prec/3.0_Prec
Real(kind=Prec) :: ay = 1.0_Prec/4.0_Prec
Real(kind=Prec) :: axy = 1.0_Prec/2.0_Prec

End Module

!******************************************************************************
Module Solution

Use Select_Precision
Use Set_Inputs

! Dimension solution variable arrays

Implicit None

Save

! Declare Arrays:

Real(kind=Prec),Dimension(imax,jmax) :: u = -99.9_Prec  ! Velocity
Real(kind=Prec),Dimension(imax,jmax) :: u_Old = -99.9_Prec  ! Old velocity
Real(kind=Prec),Dimension(imax,jmax) :: dt = -99.9_Prec  ! Time step
Real(kind=Prec),Dimension(imax,jmax) :: u_MMS = -99.9_Prec  ! MMS Exact Solution (Velocity)
Real(kind=Prec),Dimension(imax,jmax) :: Source_MMS = -99.9_Prec  ! MMS Source Terms
!Real(kind=Prec),Dimension(imax,jmax) :: F1_Flux = -99.9_Prec  ! Flux for xsi-direction at i+1/2,j
!Real(kind=Prec),Dimension(imax,jmax) :: G1_Flux = -99.9_Prec  ! Flux for eta-direction at i,j+1/2

End Module

!******************************************************************************
!INFORMATION ON ALL MODULES
! 1 --> Select_Precision
! 2 --> Set_Inputs   (contains Select_Precision)
! 3 --> Geometry     (contains Select_Precision, Set_Inputs)
! 4 --> MMS_Constants(contains Select_Precision)
! 5 --> Solution     (contains Select_Precision, Set_Inputs)  
!******************************************************************************

!******************************************************************************
!**************************    SUBROUTINES    *********************************
!******************************************************************************

Subroutine Write_Precision_Info

Use Select_Precision

  Implicit None

  ! Declare variables of each type
  Real(kind=Single) :: var1 = 0.0
  Real(kind=Double) :: var2 = 0.0_Double

  ! Write characteristics of selected variables
  Write (*,*) '********************************************************'
  Write (*,*) '********************************************************'
  Write (*,*)
  Write (*,100) 'Single precision',Kind(var1),Precision(var1),Range(var1)
  Write (*,100) 'Double precision',Kind(var2),Precision(var2),Range(var2)
  ! Write (*,100) 'var3',Kind(var3),Precision(var3),Range(var3)
  100 Format(1X,A,': kind = ',I2,', Precision = ',I2,', Range = ',I3)

  ! Let user know if quad precision is available
  If (Quad == -1) then
    Write (*,*) '** Quad precision not available **'
  Else
    Write (*,*) '** Quad precision is available on this system **'
  End If

  ! Let user know precision for current computations
  Write (*,*) 
  Write (*,*) '********************************'
  If (Prec == Single) then
    Write (*,*) 'Currently using Single precision'
  Else If (Prec == Double) then
    Write (*,*) 'Currently using Double precision'
  Else
    Write (*,*) 'Error: Precision must be specified as single or double!'
  Endif
  Write (*,*) '********************************'
  Write (*,*)

End Subroutine

!******************************************************************************

Subroutine Initialize_Constants

!Use Select_Precision
Use Set_Inputs

Implicit None

Save

Pi = ACOS(-1.0_Prec) ! Set Pi

End Subroutine

!******************************************************************************

Subroutine Set_Geometry

!Use Select_Precision
Use Geometry

  Implicit None

  Integer         ::  i_Cartesian = 0  ! Cartesian switch: =0 to read in grid, =1 for Cartesian
  Integer         ::  nzones = 1       ! Number of zones in the grid (currently must be = 1)
  Integer         ::  imax_grid        ! Maximum i index for grid file grid.dat
  Integer         ::  jmax_grid        ! Maximum j index for grid file grid.dat
  Integer         ::  kmax_grid        ! Maximum k index for grid file grid.dat (not used for 2D)
  Integer         ::  k                ! k index (not used for 2D)
  Real(kind=Prec) :: zztemp = 0.0_Prec ! Default z-coordinate
  Real(kind=Prec) :: theta = 0.0_Prec  ! Cartesian grid rotation angle
  Real(kind=Prec) :: theta_local = 0.0_Prec  ! Grid curvature angle
  Real(kind=Prec) :: skew = 0.0_Prec   ! Cartesian grid skewness factor
  Real(kind=Prec) :: xold = 0.0_Prec   ! Temporary storage for old x coordinate
  Real(kind=Prec) :: yold = 0.0_Prec   ! Temporary storage for old y coordinate



  dx = (xmax - xmin)/dfloat(imax - 1)
  dy = (ymax - ymin)/dfloat(jmax - 1)

  !Write out info
  write(*,*) 'imax,jmax: ',imax,jmax
  write(*,*) 'xmin,xmax: ',xmin,xmax
  write(*,*) 'ymin,ymax: ',ymin,ymax
  write(*,*) 'dx,dy: ',dx,dy
  write(*,*) 'rLength: ',rLength

  !Cartesian
  !Specify x and y coordinates

  if(i_Cartesian == 1) then
    ! Create Cartesian grid from zero to unity
    do i = 1, imax
    do j = 1, jmax
      x(i,j) = xmin + dfloat(i - 1)/dfloat(imax - 1)*(xmax - xmin)
      y(i,j) = ymin + dfloat(j - 1)/dfloat(jmax - 1)*(ymax - ymin)
    enddo
    enddo

    ! Grid skewing by factor "skew"
    skew = 0.3  ! <<<<<<<<<**************  GRID SKEW !!!!!!!!
    do i = 1, imax
    do j = 1, jmax
      xold = x(i,j)
      yold = y(i,j)
      x(i,j) = xold + skew*yold
      y(i,j) = yold + skew*xold
    enddo
    enddo

    ! Coordinate rotation by angle theta
    theta = Pi/7.0_Prec ! Radians  <<<<<<<<<**************  GRID ROTATION !!!!!!!!
    do i = 1, imax
    do j = 1, jmax
      xold = x(i,j)
      yold = y(i,j)
      x(i,j) = xold*cos(theta) - yold*sin(theta)
      y(i,j) = xold*sin(theta) + yold*cos(theta)
    enddo
    enddo
    
    ! Grid curvature by angle theta
    theta = -Pi/9.0_Prec ! Radians  <<<<<<<<<**************  GRID ROTATION !!!!!!!!
    do i = 1, imax
    do j = 1, jmax
      theta_local = theta*float(i-1)/float(imax-1)
      xold = x(i,j)
      yold = y(i,j)
      x(i,j) = xold*cos(theta_local) - yold*sin(theta_local)
      y(i,j) = xold*sin(theta_local) + yold*cos(theta_local)
    enddo
    enddo
    
  else if (i_Cartesian == 0) then
    ! Or, read in grid coordinates from file grid.dat
    open(unit=12,file='grid.grd')
    read(12,*) nzones
    if(nzones.ne.1) then
      write(*,*) 'ERROR: Only one zone allowed!!!'
      stop
    endif
    read(12,*) imax_grid,jmax_grid,kmax_grid
!    if( (imax /= imax_grid).or.(jmax /= jmax_grid).or.(kmax_grid /= 1) ) then	!Original
    if( (imax /= imax_grid).or.(jmax /= jmax_grid) ) then  !Removed kmax = 1 requirement
      write(*,*) 'Error: size of grid in grid.dat differs from input value !!!'
      write(*,*) 'grid.dat: imax,jmax = ',imax_grid,jmax_grid
      write(*,*) 'Input: imax,jmax = ',imax,jmax
      write(*,*) 'kmax_grid = ',kmax_grid
      stop
    else
      read(12,*) (((x(i,j),i=1,imax),j=1,jmax),k=1,kmax_grid),  &
                 (((y(i,j),i=1,imax),j=1,jmax),k=1,kmax_grid),  &
                 (((zztemp,i=1,imax),j=1,jmax),k=1,kmax_grid)
    endif
  else
    write(*,*) 'Error: i_Cartesian must be zero or one !!!'
    stop
  endif
  
!  Calculate metrics and Jacobians
!       Notes: *metrics are actually xsix/J, etax/J, etc.
!              *vol(i,j) is the cell volume or inverse Jacobian
!     Metrics/Jacobians at Interior Points
  do i = 2, imax-1
  do j = 2, jmax-1
    xsix(i,j) =  0.5*(y(i,j+1) - y(i,j-1))
    xsiy(i,j) = -0.5*(x(i,j+1) - x(i,j-1))
    etax(i,j) = -0.5*(y(i+1,j) - y(i-1,j))
    etay(i,j) =  0.5*(x(i+1,j) - x(i-1,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
!    write(*,*) '**i,j: ',i,j
!    write(*,*) 'xsix,xsiy,etax,etay,vol: ',xsix(i,j),xsiy(i,j),
! &              etax(i,j),etay(i,j),vol(i,j)
  enddo
  enddo
      
!     Metrics/Jacobians at Boundaries
  do i = 2, imax-1
    j = 1
    xsix(i,j) =  0.5*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2))
    xsiy(i,j) = -0.5*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2))
    etax(i,j) = -0.5*(y(i+1,j) - y(i-1,j))
    etay(i,j) =  0.5*(x(i+1,j) - x(i-1,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
    j = jmax
    xsix(i,j) =  0.5*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2))
    xsiy(i,j) = -0.5*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2))
    etax(i,j) = -0.5*(y(i+1,j) - y(i-1,j))
    etay(i,j) =  0.5*(x(i+1,j) - x(i-1,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  enddo
  do j = 2, jmax-1    
    i = 1
    xsix(i,j) =  0.5*(y(i,j+1) - y(i,j-1))
    xsiy(i,j) = -0.5*(x(i,j+1) - x(i,j-1))
    etax(i,j) = -0.5*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j))
    etay(i,j) =  0.5*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
    i = imax
    xsix(i,j) =  0.5*(y(i,j+1) - y(i,j-1))
    xsiy(i,j) = -0.5*(x(i,j+1) - x(i,j-1))
    etax(i,j) = -0.5*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j))
    etay(i,j) =  0.5*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j))
    vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  enddo
! Metrics/Jacobians at Corners
  i = 1
  j = 1
  xsix(i,j) =  0.5*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2))
  xsiy(i,j) = -0.5*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2))
  etax(i,j) = -0.5*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j))
  etay(i,j) =  0.5*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  i = 1
  j = jmax
  xsix(i,j) =  0.5*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2))
  xsiy(i,j) = -0.5*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2))
  etax(i,j) = -0.5*(-3.*y(i,j)+4.*y(i+1,j)-y(i+2,j))
  etay(i,j) =  0.5*(-3.*x(i,j)+4.*x(i+1,j)-x(i+2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  i = imax
  j = 1
  xsix(i,j) =  0.5*(-3.*y(i,j)+4.*y(i,j+1)-y(i,j+2))
  xsiy(i,j) = -0.5*(-3.*x(i,j)+4.*x(i,j+1)-x(i,j+2))
  etax(i,j) = -0.5*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j))
  etay(i,j) =  0.5*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)
  i = imax
  j = jmax
  xsix(i,j) =  0.5*(3.*y(i,j)-4.*y(i,j-1)+y(i,j-2))
  xsiy(i,j) = -0.5*(3.*x(i,j)-4.*x(i,j-1)+x(i,j-2))
  etax(i,j) = -0.5*(3.*y(i,j)-4.*y(i-1,j)+y(i-2,j))
  etay(i,j) =  0.5*(3.*x(i,j)-4.*x(i-1,j)+x(i-2,j))
  vol(i,j) = etay(i,j)*xsix(i,j) - xsiy(i,j)*etax(i,j)

  ! Metric derivatives
  do i = 2, imax-1
  do j = 2, jmax-1
    d_xsix_dxsi(i,j) = 0.5_Prec*( xsix(i+1,j)/vol(i+1,j) - xsix(i-1,j)/vol(i-1,j) )
    d_xsiy_dxsi(i,j) = 0.5_Prec*( xsiy(i+1,j)/vol(i+1,j) - xsiy(i-1,j)/vol(i-1,j) )
    d_etax_dxsi(i,j) = 0.5_Prec*( etax(i+1,j)/vol(i+1,j) - etax(i-1,j)/vol(i-1,j) )
    d_etay_dxsi(i,j) = 0.5_Prec*( etay(i+1,j)/vol(i+1,j) - etay(i-1,j)/vol(i-1,j) )
    d_xsix_deta(i,j) = 0.5_Prec*( xsix(i,j+1)/vol(i,j+1) - xsix(i,j-1)/vol(i,j-1) )
    d_xsiy_deta(i,j) = 0.5_Prec*( xsiy(i,j+1)/vol(i,j+1) - xsiy(i,j-1)/vol(i,j-1) )
    d_etax_deta(i,j) = 0.5_Prec*( etax(i,j+1)/vol(i,j+1) - etax(i,j-1)/vol(i,j-1) )
    d_etay_deta(i,j) = 0.5_Prec*( etay(i,j+1)/vol(i,j+1) - etay(i,j-1)/vol(i,j-1) )
  enddo
  enddo

! Write out initial metrics
  open(50,file='gridinfo.dat',status='unknown')
  write(50,*)'TITLE = "2D Heat Conduction Grid Transform. Data"'
  write(50,*)'variables="x(m)""y(m)""Vol""xsix""xsiy""etax""etay"'
  write(50,*) 'zone T="n = 0 " '
  write(50,*) 'I=',imax,' J=',jmax
  write(50,*) 'DATAPACKING=POINT'
  write(50,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE)'
  do j = 1, jmax
  do i = 1, imax
    write(50,*)x(i,j),y(i,j),vol(i,j),xsix(i,j),xsiy(i,j),etax(i,j),etay(i,j)
  enddo
  enddo
  close(50)

! Write out additional metric info
  open(51,file='highermetricinfo.dat',status='unknown')
  write(51,*)'TITLE = "2D Burgers Equation Higher Metric Info"'
  write(51,*)'variables="x(m)""y(m)""d_xsix_dxsi""d_xsiy_dxsi" ', &
             '"d_etax_dxsi""d_etay_dxsi""d_xsix_deta""d_xsiy_deta" '
  write(51,*)'"d_etax_deta""d_etay_deta"'

  write(51,*) 'zone T="n = 0 " '
  write(51,*) 'I=',imax-2,' J=',jmax-2
  write(51,*) 'DATAPACKING=POINT'
  write(51,*) 'DT=(DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE DOUBLE ', & 
                   'DOUBLE DOUBLE DOUBLE )'

  do j = 2, jmax-1
  do i = 2, imax-1
    write(51,*)x(i,j),y(i,j),d_xsix_dxsi(i,j),d_xsiy_dxsi(i,j) &
              ,d_etax_dxsi(i,j),d_etay_dxsi(i,j),d_xsix_deta(i,j) &
              ,d_xsiy_deta(i,j),d_etax_deta(i,j),d_etay_deta(i,j) 
  enddo
  enddo
  close(51)

End Subroutine

!******************************************************************************

Subroutine Set_Initial_Conditions

!Use Select_Precision
!Use Set_Inputs
Use Geometry
Use Solution

  Implicit None 

  ! Set initial conditions (u = 0.0)
  do i = 1, imax
  do j = 1, jmax
    u(i,j) = 2.0_Prec - 4.0_Prec*float(i-1)/float(imax-1) ! Linear IC in xsi direction
   ! u(i,j) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(i,j) + y(i,j))*alpha/rLength) ! Exact IC
    u_Old(i,j) = u(i,j)
  enddo
  enddo

End Subroutine

!******************************************************************************

Subroutine Set_Boundary_Conditions

!Use Select_Precision
!Use Set_Inputs
Use Geometry
Use Solution

  Implicit None

  ! Set boundary conditions (top and bottom boundaries)
  do i = 2, imax-1
    u(i,1) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(i,1) + y(i,1))*alpha/rLength)
    u(i,jmax) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(i,jmax) + y(i,jmax))*alpha/rLength)
  enddo
  ! Set boundary conditions (left and right boundaries, including corners)
  do j = 1, jmax
    u(1,j) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(1,j) + y(1,j))*alpha/rLength)
    u(imax,j) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(imax,j) + y(imax,j))*alpha/rLength)
  enddo

End Subroutine

!******************************************************************************

Subroutine Output_Solution(n)

!Use Select_Precision
!Use Set_Inputs
Use Geometry
Use Solution

  Implicit None
  
  Integer ::  n    ! Iteration number
  Real(kind=Prec),dimension(imax,jmax) :: Trunc_Error = 0.0_Prec     ! Truncation Error
  
  call Truncation_error(Trunc_error)

! Write out solution at iteration n
  write(40,*) 'zone T="',n,'" '
  write(40,*) 'I=',imax,' J=',jmax
  write(40,*) 'DATAPACKING=POINT'
  do j = 1, jmax
  do i = 1, imax
    write(40,*)x(i,j),y(i,j),u(i,j),u_MMS(i,j),(u(i,j)-u_MMS(i,j)),Trunc_error(i,j)
  enddo
  enddo

End Subroutine

!******************************************************************************

Subroutine Set_Local_Time_Step(dt_min)

!Use Select_Precision
!Use Set_Inputs
Use Geometry
Use Solution

  Implicit None

  Real(kind=Prec) ::  dt_min                    ! Minimum time step in domain
  Real(kind=Prec) ::  dt_diffusive = 1.e99_Prec ! Diffusive time step
  Real(kind=Prec) ::  dt_convective = 1.e99_Prec ! Diffusive time step
  Real(kind=Prec) ::  dist_xsi                  ! Temporary variable
  Real(kind=Prec) ::  dist_eta                  ! Temporary variable
  Real(kind=Prec) ::  dist	                  	! Temporary variable
  
!  Calculate time step        
  dt_min = 1.e99_Prec      
  do j = 2, jmax - 1
  do i = 2, imax - 1
    dist_eta = sqrt( xsix(i,j)**2 + xsiy(i,j)**2 )
    dist_xsi = sqrt( etax(i,j)**2 + etay(i,j)**2 )
    dist = amin1(dist_eta, dist_xsi)
    dt_diffusive = cfl*0.5_Prec*dist*dist/rnu  
    dt_convective = cfl*dist/amax1(abs(u(i,j)),0.000001_Prec)
    dt(i,j) = amin1(dt_max,amin1(dt_diffusive, dt_convective))
    dt_min = amin1(dt_min,dt(i,j))
  enddo
  enddo

End Subroutine

!******************************************************************************

Subroutine Euler_Explicit_Iteration

!Use Select_Precision
!Use Set_Inputs
Use Geometry
Use Solution

  Implicit None

  Real(kind=Prec),Dimension(imax,jmax) :: u_New = -99.9_Prec  ! New Velocity
  Real(kind=Prec) ::  dudt     	  = -99.0_Prec     	! du/dt derivative

  Real(kind=Prec) ::  dudx      = -99.0_Prec      ! du/dx
  Real(kind=Prec) ::  dudy      = -99.0_Prec      ! du/dy
  Real(kind=Prec) ::  d2udx2    = -99.0_Prec      ! d2u/dx2
  Real(kind=Prec) ::  d2udy2    = -99.0_Prec      ! d2u/dy2
  Real(kind=Prec) ::  d2udx2CART    = -99.0_Prec      ! d2u/dx2
  Real(kind=Prec) ::  d2udy2CART    = -99.0_Prec      ! d2u/dy2

  Real(kind=Prec) ::  Laplacian   = -99.0_Prec      ! d2u/dx2 + d2T/dy2 (computed from transformed coordinates)
  Real(kind=Prec) ::  dudxsi      = -99.0_Prec      ! du/dxsi
  Real(kind=Prec) ::  dudeta      = -99.0_Prec      ! du/deta
  Real(kind=Prec) ::  d2udxsi2    = -99.0_Prec      ! d2u/dxsi2
  Real(kind=Prec) ::  d2udeta2    = -99.0_Prec      ! d2u/deta2
  Real(kind=Prec) ::  d2udxsideta = -99.0_Prec      ! d2u/dxsi/deta (cross derivative term)
                                                    ! ***Note: Below, xsix is NOT xsix/J !!!***
  Real(kind=Prec) ::  term1       = -99.0_Prec      ! Temporary Variables
  Real(kind=Prec) ::  term2       = -99.0_Prec      ! Temporary Variables
  Real(kind=Prec) ::  term3       = -99.0_Prec      ! Temporary Variables

!  Save old u values at time level n
  do j = 1, jmax
  do i = 1, imax
    u_Old(i,j) = u(i,j)
  enddo
  enddo

!  Apply the Euler explicit method for a Jacobi-like iteration

! Euler Explicit Method 
! Forward Sweep
  do j = 2, jmax - 1
  do i = 2, imax - 1
    dudxsi = 0.5_Prec*( u_Old(i+1,j) - u_Old(i-1,j) )
    dudeta = 0.5_Prec*( u_Old(i,j+1) - u_Old(i,j-1) )
    dudx = (xsix(i,j)*dudxsi + etax(i,j)*dudeta)/vol(i,j)
    dudy = (xsiy(i,j)*dudxsi + etay(i,j)*dudeta)/vol(i,j)

    d2udxsi2 =  u_Old(i+1,j) - 2.0_Prec*u_Old(i,j) + u_Old(i-1,j) 
    d2udeta2 =  u_Old(i,j+1) - 2.0_Prec*u_Old(i,j) + u_Old(i,j-1) 
    d2udxsideta = 0.5_Prec*( 0.5_Prec*( u_Old(i+1,j+1) - u_Old(i+1,j-1) )  &
                           - 0.5_Prec*( u_Old(i-1,j+1) - u_Old(i-1,j-1) ) )

    ! Terms for d2u/dx2
    term1 =   ( xsix(i,j)**2 )*d2udxsi2  &
            + 2.0_Prec*( etax(i,j)*xsix(i,j) )*d2udxsideta  &
            + ( etax(i,j)**2 )*d2udeta2
    term1 = term1/( vol(i,j)**2 )
    term2 = dudxsi*( xsix(i,j)*d_xsix_dxsi(i,j)  &
                   + etax(i,j)*d_xsix_deta(i,j) ) / vol(i,j)
    term3 = dudeta*( xsix(i,j)*d_etax_dxsi(i,j)  &
                   + etax(i,j)*d_etax_deta(i,j) ) / vol(i,j)
    d2udx2 = term1 + term2 + term3

    term1 =   ( xsiy(i,j)**2 )*d2udxsi2  &
            + 2.0_Prec*( etay(i,j)*xsiy(i,j) )*d2udxsideta  &
            + ( etay(i,j)**2 )*d2udeta2
    term1 = term1/( vol(i,j)**2 )
    term2 = dudxsi*( xsiy(i,j)*d_xsiy_dxsi(i,j)  &
                   + etay(i,j)*d_xsiy_deta(i,j) ) / vol(i,j)
    term3 = dudeta*( xsiy(i,j)*d_etay_dxsi(i,j)  &
                   + etay(i,j)*d_etay_deta(i,j) ) / vol(i,j)
    d2udy2 = term1 + term2 + term3

    dudt = -(u_Old(i,j)*dudx + u_Old(i,j)*dudy - rnu*(d2udx2 + d2udy2))
    u(i,j) = u_Old(i,j) + dudt*dt(i,j)
  enddo
  enddo

!$$$$$$   Call Set_MMS_Boundary_Conditions
!$$$$$$ !  Save old u values at time level n
!$$$$$$   do j = 1, jmax
!$$$$$$   do i = 1, imax
!$$$$$$     u_Old(i,j) = u(i,j)
!$$$$$$   enddo
!$$$$$$   enddo
!$$$$$$ 
!$$$$$$ ! Backward Sweep
!$$$$$$   do i = imax - 1, 2, -1  
!$$$$$$   do j = jmax - 1, 2, -1 
!$$$$$$     dudxsi = 0.5_Prec*( u_Old(i+1,j) - u_Old(i-1,j) )
!$$$$$$     dudeta = 0.5_Prec*( u_Old(i,j+1) - u_Old(i,j-1) )
!$$$$$$     dudx = (xsix(i,j)*dudxsi + etax(i,j)*dudeta)/vol(i,j)
!$$$$$$     dudy = (xsiy(i,j)*dudxsi + etay(i,j)*dudeta)/vol(i,j)
!$$$$$$ 
!$$$$$$     d2udxsi2 =  u_Old(i+1,j) - 2.0_Prec*u_Old(i,j) + u_Old(i-1,j) 
!$$$$$$     d2udeta2 =  u_Old(i,j+1) - 2.0_Prec*u_Old(i,j) + u_Old(i,j-1) 
!$$$$$$     d2udxsideta = 0.5_Prec*( 0.5_Prec*( u_Old(i+1,j+1) - u_Old(i+1,j-1) )  &
!$$$$$$                            - 0.5_Prec*( u_Old(i-1,j+1) - u_Old(i-1,j-1) ) )
!$$$$$$     ! Terms for d2u/dx2
!$$$$$$     term1 =   ( xsix(i,j)**2 )*d2udxsi2  &
!$$$$$$             + 2.0_Prec*( etax(i,j)*xsix(i,j) )*d2udxsideta  &
!$$$$$$             + ( etax(i,j)**2 )*d2udeta2
!$$$$$$     term1 = term1/( vol(i,j)**2 )
!$$$$$$     term2 = dudxsi*( xsix(i,j)*d_xsix_dxsi(i,j)  &
!$$$$$$                    + etax(i,j)*d_xsix_deta(i,j) ) / vol(i,j)
!$$$$$$     term3 = dudeta*( xsix(i,j)*d_etax_dxsi(i,j)  &
!$$$$$$                    + etax(i,j)*d_etax_deta(i,j) ) / vol(i,j)
!$$$$$$     d2udx2 = term1 + term2 + term3
!$$$$$$ 
!$$$$$$     term1 =   ( xsiy(i,j)**2 )*d2udxsi2  &
!$$$$$$             + 2.0_Prec*( etay(i,j)*xsiy(i,j) )*d2udxsideta  &
!$$$$$$             + ( etay(i,j)**2 )*d2udeta2
!$$$$$$     term1 = term1/( vol(i,j)**2 )
!$$$$$$     term2 = dudxsi*( xsiy(i,j)*d_xsiy_dxsi(i,j)  &
!$$$$$$                    + etay(i,j)*d_xsiy_deta(i,j) ) / vol(i,j)
!$$$$$$     term3 = dudeta*( xsiy(i,j)*d_etay_dxsi(i,j)  &
!$$$$$$                    + etay(i,j)*d_etay_deta(i,j) ) / vol(i,j)
!$$$$$$     d2udy2 = term1 + term2 + term3
!$$$$$$ 
!$$$$$$     dudt = -(u_Old(i,j)*dudx + u_Old(i,j)*dudy - rnu*(d2udx2 + d2udy2))
!$$$$$$     u(i,j) = u_Old(i,j) + dudt*dt(i,j)
!$$$$$$   enddo
!$$$$$$   enddo


End Subroutine


Subroutine Truncation_Error(TE)

!Use Select_Precision
!Use Set_Inputs
Use Geometry
Use Solution

  Implicit None

  Real(kind=Prec),Dimension(imax,jmax) :: uloc = -99.9_Prec  ! New Velocity
  Real(kind=Prec),Dimension(imax,jmax),intent(out) :: TE !truncation error

  Real(kind=Prec) ::  dudx      = -99.0_Prec      ! du/dx
  Real(kind=Prec) ::  dudy      = -99.0_Prec      ! du/dy
  Real(kind=Prec) ::  d2udx2    = -99.0_Prec      ! d2u/dx2
  Real(kind=Prec) ::  d2udy2    = -99.0_Prec      ! d2u/dy2
  Real(kind=Prec) ::  d2udx2CART    = -99.0_Prec      ! d2u/dx2
  Real(kind=Prec) ::  d2udy2CART    = -99.0_Prec      ! d2u/dy2

  Real(kind=Prec) ::  Laplacian   = -99.0_Prec      ! d2u/dx2 + d2u/dy2 (computed from transformed coordinates)
  Real(kind=Prec) ::  dudxsi      = -99.0_Prec      ! du/dxsi
  Real(kind=Prec) ::  dudeta      = -99.0_Prec      ! du/deta
  Real(kind=Prec) ::  d2udxsi2    = -99.0_Prec      ! d2u/dxsi2
  Real(kind=Prec) ::  d2udeta2    = -99.0_Prec      ! d2u/deta2
  Real(kind=Prec) ::  d2udxsideta = -99.0_Prec      ! d2u/dxsi/deta (cross derivative term)
                                                    ! ***Note: Below, xsix is NOT xsix/J !!!***
  Real(kind=Prec) ::  term1       = -99.0_Prec      ! Temporary Variables
  Real(kind=Prec) ::  term2       = -99.0_Prec      ! Temporary Variables
  Real(kind=Prec) ::  term3       = -99.0_Prec      ! Temporary Variables

!  Save old u values at time level n
  do j = 1, jmax
  do i = 1, imax
    uloc(i,j) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(i,j) + y(i,j))*alpha/rLength) ! Exact IC
  enddo
  enddo

!  Apply the Euler explicit method for a Jacobi-like iteration

! Euler Explicit Method 
! Forward Sweep
  do j = 2, jmax - 1
  do i = 2, imax - 1
    dudxsi = 0.5_Prec*( uloc(i+1,j) - uloc(i-1,j) )
    dudeta = 0.5_Prec*( uloc(i,j+1) - uloc(i,j-1) )
    dudx = (xsix(i,j)*dudxsi + etax(i,j)*dudeta)/vol(i,j)
    dudy = (xsiy(i,j)*dudxsi + etay(i,j)*dudeta)/vol(i,j)

    d2udxsi2 =  uloc(i+1,j) - 2.0_Prec*uloc(i,j) + uloc(i-1,j) 
    d2udeta2 =  uloc(i,j+1) - 2.0_Prec*uloc(i,j) + uloc(i,j-1) 
    d2udxsideta = 0.5_Prec*( 0.5_Prec*( uloc(i+1,j+1) - uloc(i+1,j-1) )  &
                           - 0.5_Prec*( uloc(i-1,j+1) - uloc(i-1,j-1) ) )

    ! Terms for d2u/dx2
    term1 =   ( xsix(i,j)**2 )*d2udxsi2  &
            + 2.0_Prec*( etax(i,j)*xsix(i,j) )*d2udxsideta  &
            + ( etax(i,j)**2 )*d2udeta2
    term1 = term1/( vol(i,j)**2 )
    term2 = dudxsi*( xsix(i,j)*d_xsix_dxsi(i,j)  &
                   + etax(i,j)*d_xsix_deta(i,j) ) / vol(i,j)
    term3 = dudeta*( xsix(i,j)*d_etax_dxsi(i,j)  &
                   + etax(i,j)*d_etax_deta(i,j) ) / vol(i,j)
    d2udx2 = term1 + term2 + term3


    term1 =   ( xsiy(i,j)**2 )*d2udxsi2  &
            + 2.0_Prec*( etay(i,j)*xsiy(i,j) )*d2udxsideta  &
            + ( etay(i,j)**2 )*d2udeta2
    term1 = term1/( vol(i,j)**2 )
    term2 = dudxsi*( xsiy(i,j)*d_xsiy_dxsi(i,j)  &
                   + etay(i,j)*d_xsiy_deta(i,j) ) / vol(i,j)
    term3 = dudeta*( xsiy(i,j)*d_etay_dxsi(i,j)  &
                   + etay(i,j)*d_etay_deta(i,j) ) / vol(i,j)
    d2udy2 = term1 + term2 + term3

    TE(i,j) = (uloc(i,j)*dudx + uloc(i,j)*dudy - rnu*(d2udx2 + d2udy2))

  enddo
  enddo



End Subroutine

!******************************************************************************

Function MMS_Exact_Solution(rLength,x,y) ! This function calculates the MMS exact solution

!Use Select_Precision    ! Needed to set precision
Use MMS_Constants       ! Needed for MMS constants
Use Set_Inputs          ! Needed for Pi

Implicit None

Real(kind=Prec) ::  MMS_Exact_Solution      ! Exaxt MMS solution
Real(kind=Prec) ::  x						! x-coordinate				!NOTE: these 3 variables are local
Real(kind=Prec) ::  y						! y-coordinate
Real(kind=Prec) ::  rLength					! reference length


  MMS_Exact_Solution = -2.0_Prec*rnu*alpha/rLength*tanh( (x + y)*alpha/rLength)

End Function

!******************************************************************************

Function MMS_Source_Fcn(rLength,x,y) ! This function calculates MMS source term

!Use Select_Precision
Use MMS_Constants
Use Set_Inputs

Implicit None

Real(kind=Prec) :: MMS_Source_Fcn        ! MMS source term
Real(kind=Prec) :: term1 = 0.0_Prec      ! Temporary variables (term1 - term4)    !NOTE: next 7 variables are local
Real(kind=Prec) :: term2 = 0.0_Prec
Real(kind=Prec) :: term3 = 0.0_Prec
Real(kind=Prec) :: term4 = 0.0_Prec
Real(kind=Prec) :: rLength               ! Reference length
Real(kind=Prec) :: x                     ! x-coordinate
Real(kind=Prec) :: y                     ! y-coordinate


  MMS_Source_Fcn = 0.0_Prec
End Function

!******************************************************************************

Subroutine Set_MMS_Boundary_Conditions

!Use Select_Precision
Use Geometry
Use Solution

Implicit None

  ! Set boundary conditions (top and bottom boundaries)
  do i = 2, imax-1
    u(i,1) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(i,1) + y(i,1))*alpha/rLength)
    u(i,jmax) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(i,jmax) + y(i,jmax))*alpha/rLength)
  enddo
  ! Set boundary conditions (left and right boundaries, including corners)
  do j = 1, jmax
    u(1,j) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(1,j) + y(1,j))*alpha/rLength)
    u(imax,j) = -2.0_Prec*rnu*alpha/rLength*tanh( (x(imax,j) + y(imax,j))*alpha/rLength)
  enddo

End Subroutine

!******************************************************************************

Subroutine Grid_Convergence

!Use Select_Precision
Use Geometry
Use Solution

Implicit None             

Real(kind=Prec) :: Disc_Error = 0.0_Prec      ! Absolute value of the discretization error 
Real(kind=Prec) :: L_Infinity_Norm = 0.0_Prec ! L_infinity norm of the discretization error
Real(kind=Prec) :: L_1_Norm = 0.0_Prec        ! L_1 norm of the discretization error
Real(kind=Prec) :: L_2_Norm = 0.0_Prec        ! L_2 norm of the discretization error
          
do j = 1, jmax
do i = 1, imax
  Disc_Error = abs( u(i,j)-u_MMS(i,j) )
  L_1_Norm = L_1_Norm + Disc_Error
  L_2_Norm = L_2_Norm + Disc_Error**2
  L_Infinity_Norm = amax1( L_Infinity_Norm, Disc_Error )
enddo
enddo

L_1_Norm = L_1_Norm / float( imax*jmax )
L_2_Norm = sqrt( L_2_Norm / float( imax*jmax ) )

write (*,*) 'imax, jmax: ', imax, jmax
write (*,*) '     L1_norm, L2_norm, Linfinity_norm: ',L_1_Norm,L_2_Norm,L_Infinity_Norm
write (60,*) 'imax, jmax: ', imax, jmax
write (60,*) '     L1_norm, L2_norm, Linfinity_norm: ',L_1_Norm,L_2_Norm,L_Infinity_Norm

End Subroutine

!******************************************************************************

Subroutine Debug_Output ! This subroutine outputs diagnostic information for debugging purposes

!Use Select_Precision
Use Geometry
Use Set_Inputs
Use Solution

Implicit None

!Call Initialize_Constants

!============== DIAGNOSTIC OUTPUT ==============
! Check Inputs
Write (*,*) nmax, ' max # of steps '
Write (*,*) iterout, ' number of time steps between output'
Write (*,*) cfl, ' maximum CFL for each i,j location'
Write (*,*) toler, ' tolerance for convergence'
Write (*,*) alpha, ' thermal diffusivity (m^2/s)'
Write (*,*) Pi, ' Set Pi'
Write (*,*) fsmall, ' small parameter'
Write (*,*) dt_max, ' maximum delta t'
Write (*,*) rtime, ' set initial time to zero'

! Check Precision (single vs. double)
Write (*,*) 'Variables will be saved in the form: ',Pi
Write (*,*)
Write (*,100) 'Precision',Kind(Pi),Precision(Pi),Range(Pi)
100 Format(1X,A,': kind = ',I2,', Precision = ',I2,', Range = ',I3)

! Write solution to standard output
do i = 1, imax
do j = 1, jmax
  write(*,*) 'i,j,x,y: ',i,j,x(i,j),y(i,j)
  write(*,*) '  u, u_Old: ',u(i,j),u_Old(i,j)
enddo
enddo

End Subroutine

! *************************************************************************
! INFORMATION ON ALL SUBROUTINES
! *************************************************************************
! No.	Subroutine				Moduled Inc.	    	Subroutine Inc.
! =========================================================================
!1	Write_Precision_Info	Select_Precision					NA
!2	Initialize_Constants	Select_Precision; Set_Inputs		NA
!3	Set_Geometry			Select_Precision; Geometry			NA
!4	Set_Initial_Conditions	Select_Precision; Set_Inputs
!							Geometry; Solution					NA
!5	Set_Boundary_Conditions	Select_Precision; Set_Inputs
!							Geometry; Solution					NA
!6	Output_Solution			Select_Precision; Set_Inputs
!							Geometry; Solution					NA
!7	Set_Local_Time_Step		Select_Precision; Set_Inputs
!							Geometry; Solution					NA
!8	Set_Global_Time_Step	Select_Precision; Set_Inputs
!							Geometry; Solution					NA
!9	Euler_Explicit_Iteration Select_Precision; Set_Inputs
!							Geometry; Solution					NA
!10	Set_MMS_Boundary_Cond.	Select_Precision; Geometry
!							Solution							NA
!11	Grid_Convergence		Select_Precision; Geometry
!							Solution							NA
!12	Debug_Output			Select_Precision; Geometry
!							Set_Inputs; Solution				Initialize_Constants
!
! *************************************************************************
! INFORMATION ON ALL FUNCTIONS
! *************************************************************************
! No.	Function				Moduled Inc.	    		Subroutine Inc.
! =========================================================================
!1	MMS_Exact_Solution		Select_Precision; Set_Inputs	NA	
!							MMS_Constants					
!2	MMS_Source_Fcn			Select_Precision; Set_Inputs	NA	
!							MMS_Constants					

! *************************************************************************

! *************************************************************************
! ***********************   MAIN PROGRAM   ********************************
! *************************************************************************

Program Main

Use Select_Precision
Use Set_Inputs
Use Geometry
Use MMS_Constants
Use Solution

Implicit None

!Real(kind=Prec) :: zero = 0.0_Prec              ! = 0.0
!Real(kind=Prec) :: third = 1.0_Prec/3.0_Prec    ! = 1./3.
Real(kind=Prec) :: Residual                     ! L2 norm of steady-state residual
Real(kind=Prec) :: Initial_Residual             ! Initial steady-state residual L2 norm
Real(kind=Prec) :: dt_min = 1.e99_Prec          ! Minimum time step in domain
Real(kind=Prec) :: MMS_Exact_Solution           ! Set type for function MMS_Exact_Solution
Real(kind=Prec) :: MMS_Source_Fcn               ! Set type for function MMS_Exact_Solution
Integer ::  n                                   ! Iteration number


alpha = RE/2._prec
rnu    = rlength*uref/Re

!  Set up output files (history and profile)
open(30,file='residuals.dat',status='unknown')
write(30,*) 'TITLE = "2D Burgers Equation Residual History"'
write(30,*) 'variables="Iteration""Time(s)""Res"'
      
open(40,file='2dveloc.dat',status='unknown')
write(40,*) 'TITLE = "2D Burgers Equatino Field Data"'
write(40,*) 'variables="x(m)""y(m)""u(m/s)""u_exact(m/s)""DE""TE"'

open(60,file='DiscError.txt',status='unknown')

Write (*,*)
Write (*,*) '********************************************************'
Write (*,*) '********************************************************'
Write (*,*)

Call Write_Precision_Info

Call Initialize_Constants

Call Set_Geometry

Call Set_Initial_Conditions

if(imms == 0) then
  Call Set_Boundary_Conditions
elseif(imms == 1) then
  do j = 1, jmax
  do i = 1, imax
    u_MMS(i,j) = MMS_Exact_Solution(rLength,x(i,j),y(i,j))
  enddo
  enddo
!  Call MMS_Exact_Solution

  do j = 1, jmax
  do i = 1, imax
    Source_MMS(i,j) = MMS_Source_Fcn(rLength,x(i,j),y(i,j))
  enddo
  enddo
!  Call MMS_Source

  Call Set_MMS_Boundary_Conditions
else
  write(*,*) 'Error: imms must equal 0 or 1 !!!'
  stop
endif

Call Output_Solution(0)

if(idebug == 1) then
  Call Debug_Output
elseif (idebug /= 0) then
  write(*,*) 'Error: idebug must equal 0 or 1 !!!'
  stop
endif

!=============      
!  Main Loop
!=============

  do n = 1, nmax

    Call Set_Local_Time_Step(dt_min)

    Call Euler_Explicit_Iteration

    if(imms == 0) then
      Call Set_Boundary_Conditions
    elseif(imms == 1) then
      Call Set_MMS_Boundary_Conditions
    else
      write(*,*) 'Error: imms must equal 0 or 1 !!!'
      stop
    endif

!   Update time (time value valid for global time stepping only)  
    rtime = rtime + dt_min

!   Output every "iterout" steps
    if( (mod(n,iterout).eq.0) ) Call Output_Solution(n)

!   Iterative Convergence
    Residual = 0.0_Prec
    do j = 2, jmax - 1
    do i = 2, imax - 1
      Residual = Residual + ( (u(i,j) - u_Old(i,j))/dt(i,j) )**2
    enddo
    enddo
    Residual=amax1(Residual,1.e-20_Prec)
    Residual = sqrt(Residual/dfloat( (imax-2)*(jmax-2) ))
    if(n.eq.1) Initial_Residual = amax1(Residual,1.e-20)
    Residual = Residual/Initial_Residual

    if( (mod(n,100).eq.0).or.n.eq.1 ) then
      write (30,*) n,rtime,Residual
      write(*,300) n,rtime,dt_min,Residual
 300  format(1X,i8,2(e11.3),3(e15.6))
    endif        

    if(Residual.le.toler) then
      write(30,*) n,rtime,Residual
      goto 200 ! Exit iteration loop
    endif
  
  cfl = min(cfl*1.001,0.005)
  
  end do

Write(*,*) 'Solution failed to converge in ',nmax,' iterations!!!'

Call Output_Solution(n)

Call Grid_Convergence
stop  

200	continue
 
Write(*,*) 'Solution converged in ',n,' iterations!!!' 	 

Call Output_Solution(n)

Call Grid_Convergence

Close(30)
Close(40)
Close(60)

OPEN(500,file='finished.dat',status='unknown')
Write(500,*) '1'
close(500)

OPEN(600,file='matlaboutput.dat',status='unknown')
do j = 1, jmax
do i = 1, imax
	Write(600,*) u(i,j)
enddo
enddo
close(600)


End
