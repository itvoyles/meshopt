module burgers1d_data
  use select_precision, only : prec
! Burgers' Equation:
!
!         du    d2u 
! L(u) = u-- - v--- = 0
!         dx    dx2
!
! Nondimensional 
! xi = xi(x)


  real(prec) :: Lref  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: Re    = 8._prec
  real(prec) :: alpha = 32._prec
  real(prec) :: nu

  real(prec) :: dx, dxi

  integer    :: imax = 129


end module burgers1d_data













!*******************************************************************************
!*******************************************************************************
!                                 MAIN PROGRAM
!*******************************************************************************
!*******************************************************************************

program main
  use select_precision
  use burgers1d_data
  use burgers1d_functions
  use derivative_transform
  use fd_derivative_calc
  use fd_derivatives
  use burgers1d_testroutines
  !use functional_integral_te_squared
  use functional_integral_te_4th
  !use functional_integral_te_1st
  use calc_grid
  use functional
  implicit none

  integer :: i,j,k,l,m,n
  integer :: order = 2
  real(prec) :: pi
  real(prec), dimension(:), allocatable :: TE, TEex,dTE_dxvec,x,dJdx,x_xi,x_xixi
  character(50) :: input_file,grid_in,grid_info


  integer :: io_error
  integer :: imax_default
  integer :: tempi
  
  character(50)  :: grid_default = 'burgers1d.grd'
  character(100) :: temp
  character(20)  :: append_replace = 'append'
  
  logical :: connected
  
  integer :: fid_input = 101
  integer :: fid_newinput = 102
  integer :: fid_grid = 103
  integer :: fid_griddata = 107
  integer :: fid_functional = 105
  integer :: fid_history = 106

  real(prec) :: scaling, fnctnl



  !replace with input file -----------------------------------------------------
  !write(*,'(A)',advance='no') 'Input file name:  '
  read(*,*) input_file
  
  !create new input file if input file name is 'new' ***************************
  if (input_file .eq. 'test') then
    call run_burgers1d_tests
    call test_J_eq_integral_TE_squared
    stop
    
  
  elseif (input_file .eq. 'new') then
    write(*,'(A)',advance = 'no') 'New input file name:  '
    read(*,*) input_file
    open(fid_newinput,file=trim(input_file),status='unknown',iostat=io_error)
    
    write(*,'(A)',advance = 'no') 'Input number of nodes:  '
    read(*,*) imax_default
    
    write(*,'(A)',advance = 'no') 'Input Reynods number:  '
    read(*,*) Re    
    
    write(*,'(A)',advance = 'no') 'Input grid file name:  '
    read(*,'(A)') grid_default
    
    write(*,'(A)',advance = 'no') 'Input grid data file name:  '
    read(*,'(A)') grid_info
   
    write(*,'(A)',advance = 'no') '        append or replace?  '
    read(*,'(A)') append_replace
    
    write(*,'(A)',advance = 'no') 'Creating new input file...'
    write(fid_newinput,'(A)') trim(input_file)
    write(fid_newinput,'(I4,A)') imax_default, '     !number of grid nodes'
    write(fid_newinput,'(f6.0,A)') Re, '     !Reynolds number'
    write(fid_newinput,'(A)') trim(grid_default)// '     !burgers grid'
    write(fid_newinput,'(A)') trim(grid_info)// '     !burgers grid info'
    write(fid_newinput,'(A)') trim(append_replace)//&
     '     !"append" or "replace" grid data for grid adaption history'
    write(*,'(A)') 'Done'
    
    close(fid_newinput)
  endif
  !*****************************************************************************

  !read input file *************************************************************
  open(fid_input,file=trim(input_file),status='old',iostat=io_error)
  if (io_error .ne. 0) then
    write(*,'(A)') 'Error! Input file not found. To create new input file type'//&
                   ' "new" for the input file name.' 
    stop
  endif
  
  read(fid_input,*) !skip first line
  
  read(fid_input,*) imax !<--------input
  read(fid_input,*) re !<--------input
  
  !create default grid overridden later if grid input **************************
  allocate(x(imax))
  
  !evenly distributed grid [-4,4]
  !dx = Lref/real(imax-1,prec)
  !do i = 1,imax
  !  x(i) = -Lref/2._prec+(i-1)*dx
  !enddo

  !truncation error based grid generation
  nu    = Lref*uref/Re
  alpha = RE/2._prec
  call calc_initial_grid(x,imax,alpha,nu,Lref)

  
  !!evenly distributed grid (-1,1) with endpoints moved to -4 and 4
  !do i = 1,imax
  !  x(i) = -2._prec+(i-1)*4._prec/real(imax-1,prec)
  !enddo
  !x(1) = -Lref/2._prec
  !x(imax) = Lref/2._prec
  
  
  read(fid_input,'(A100)')temp !<--------input
  i = scan(temp,' ')
  grid_in = temp(1:i)
  
  !open grid file and create new grid if doesn't exist on query ----------------
  inquire(file=trim(grid_in), exist=connected)
  if (connected.eqv..false.) then
    write(*,'(A)', advance='no') 'Error reading grid file. '//&
                                 'Create new grid '//trim(grid_in)//' [y/n]?  '
    read(*,*) temp
    
    if (trim(temp) .eq. 'y') then
      open(fid_grid,file=trim(grid_in),status='unknown',iostat=io_error)
      write(fid_grid,*) 1
      write(fid_grid,*) imax, 1, 1

      write(fid_grid,'(4e23.14)') (x(i),i=1,imax)
      write(fid_grid,'(4e23.14)') (0._prec*real(i,prec),i=1,imax)
      write(fid_grid,'(4e23.14)') (0._prec*real(i,prec),i=1,imax)
      
      !clear variables
      close(fid_grid)
      
      write(*,'(A)') 'New grid file created.'
      
    else
      write(*,'(A)') 'Stopping program execution'
      stop
      
    endif
  
  endif
  
  
  read(fid_input,'(A100)') temp !<--------input
  i = scan(temp,' ')
  grid_info = temp(1:i)
  
  read(fid_input,'(A100)') temp !<--------input
  i = scan(temp,' ')
  append_replace = temp(1:i)
  
  
  !reads input to append or replace grid data for adaption history
  read(fid_input,'(A100)',iostat=io_error) temp !<--------input
  if (io_error .eq. 0) then
  i = scan(temp,' ')
    if ((temp(1:i) .eq. 'append') .or. (temp(1:i) .eq. 'replace')) then
      append_replace = temp(1:i)
    else
      append_replace = 'replace'
    endif
  endif
  
  !setup grid info file --------------------------------------------------------
  inquire(file=trim(grid_info), exist=connected)
  open(fid_history,file=trim(grid_info),status='unknown',iostat=io_error)
  if ((append_replace.eq.'replace').or.(connected.eqv..false.)) then
    write(fid_history,'(A)') 'Title="'//trim(grid_in)//' info"'
    write(fid_history,'(A)') 'VARIABLES="x""Xxi""Xxixi""dJdx""TE"'
  
  elseif (append_replace.eq.'append') then
    io_error = 0
    do while (io_error .eq. 0)
      read(fid_history,*,iostat=io_error) temp
    enddo
  endif
  
  write(fid_history,'(A)') 'zone T="grid"'
  write(fid_history,*) 'i=',imax
  
  
  !read grid file --------------------------------------------------------------
  open(fid_grid,file=trim(grid_in),status='unknown',iostat=io_error)
  read(fid_grid,*)
  read(fid_grid,*) tempi
  if (tempi .ne. imax) then
    write(*,'(A)') 'Error! Grid size does not match input file size.'
    stop
  endif
  
  allocate( x_xi(imax),x_xixi(imax),&
            te(imax),teex(imax),dJdx(imax),dTE_dxvec(imax) )
  
  read(fid_grid,*) (x(i),i=1,imax)
  read(fid_grid,*) 
  read(fid_grid,*) 
     

  dxi = 1._prec
  pi = acos(-1._prec)
  
  !-----------------------------------------------------------------------------
 





  !calculate grid metrics ------------------------------------------------------
  call dnfdxn(x_xi,x,dxi,1,order,imax)
  call dnfdxn(x_xixi,x,dxi,2,order,imax)

  !calculates the exact truncation from the centered, finite-difference residual
  TEex = Burgers_fd_cont_resid(x,imax,nu,alpha,Lref)

  !Derived truncation error
  TE = Burgers_TE(x,x_xi,x_xixi,dxi,nu,alpha,Lref)
  

  
  !calculate funcational sensitivity
  call calc_dJdx(djdx,x,imax,nu,alpha,Lref)
  
  fnctnl = functional_J(x,TE,imax)

  !write grid info to tecplot format  
  open(fid_griddata,file='griddata.dat',status='unknown')
  write(fid_griddata,*) imax
  write(fid_griddata,'(6e23.14)') (x(i),djdx(i),x_xi(i),-1._prec,0._prec,1._prec,i=1,imax)

  !write functional value to file
  open(fid_functional,file='functional.dat',status='unknown')
  write(fid_functional,*) fnctnl
  close(fid_functional)
  
  do i = 1,imax
    write(fid_history,'(5e23.14)') x(i),x_xi(i),x_xixi(i),djdx(i),TE(i)
  enddo

end program main
