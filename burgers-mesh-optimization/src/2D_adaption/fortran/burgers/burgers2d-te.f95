module burgers2d_data
  use select_precision, only : prec
! Burgers' Equation:
!
!         du    du    d2u    d2u
! L(u) = u-- + u-- - v--- - v--- = 0
!         dx    dy    dx2    dy2
!
! Nondimensional 
! xi = xi(x)


  real(prec) :: Lref  = 8._prec
  real(prec) :: uref  = 2._prec
  real(prec) :: Re    = 16._prec
  real(prec) :: alpha 
  real(prec) :: nu




end module burgers2d_data



program main
  use select_precision, only : prec
  use burgers2d_data
  use burgers2d_functions
  use functional
  use ghetto_adapt
  implicit none

  integer:: imax = 9
  integer:: jmax = 9
  real(prec), allocatable, dimension(:,:) :: x,y,xnew,ynew
  real(prec), allocatable, dimension(:,:) :: djdx,djdy,TE,u
  real(prec), allocatable, dimension(:,:,:)::djdxmat
  real(prec) :: Jfunc
  
  real(prec) :: dx,dy
  integer :: i,j, io_error



  logical :: connected
  
  character(100) :: temp
  
  integer :: fid_grid = 101
  character(50)  :: grid_in = 'grid.grd'
  
  integer :: fid_input = 102
  character(50)  :: input_file
  
  integer :: fid_gridinfo = 103
  character(50)  :: grid_info = 'gridhistory.dat'
  character(10)  :: append_replace ='append'


  !read input file *************************************************************
  read(*,*) input_file
  
  inquire(file=trim(input_file),exist=connected)
  open(fid_input,file=trim(input_file),status='unknown')
  if (connected) then
    read(fid_input,*)
    read(fid_input,*) imax
    read(fid_input,*) jmax
    read(fid_input,*) RE
    
    alpha = RE/2._prec
    nu    = Lref*uref/Re

    read(fid_input,'(A100)') temp !<--------input grid
    i = scan(temp,' ')
    grid_in = temp(1:i)
    
    read(fid_input,'(A100)') temp !<--------input grid info
    i = scan(temp,' ')
    grid_info = temp(1:i)
  
    read(fid_input,'(A100)') temp !<--------input grid info option
    i = scan(temp,' ')
    append_replace = temp(1:i)
    
  else
    write(*,'(A)',advance='no') 'Input file doesn''t exist. Write new [y/n]?  '
    read(*,*) temp
    
    if (trim(temp).eq.'y') then
      write(fid_input,'(A)') trim(input_file)
      write(fid_input,'(I5)') imax
      write(fid_input,'(I5)') jmax
      write(fid_input,'(f6.0)') RE
      write(fid_input,'(A)') trim(grid_in)
      write(fid_input,'(A)') trim(grid_info)
      write(fid_input,'(A)') trim(append_replace)
    endif
    
    stop 
  
  endif 
  close(fid_input)
!*******************************************************************************
!******************************************************************************* 




 




!Read grid file ****************************************************************
  inquire(file=trim(grid_in),exist=connected)
  open(fid_grid,file=trim(grid_in),status='unknown')
  if (connected) then
    read(fid_grid,*)
    read(fid_grid,*) imax,jmax
    
    allocate(x(imax,jmax),y(imax,jmax))
             
    read(fid_grid,*) ((x(i,j),i=1,imax),j=1,jmax)
    read(fid_grid,*) ((y(i,j),i=1,imax),j=1,jmax)
  
  else
    write(*,'(A)',advance='no') 'Grid file doesn''t exist. Write new [y/n]?  '
    read(*,*) temp
    
    !create evenly distributed grid
    if (trim(temp).eq.'y') then
      allocate(x(imax,jmax),y(imax,jmax))

      dx = Lref/real(imax-1,prec)
      dy = Lref/real(jmax-1,prec)
      do j = 1,jmax
       do i = 1,imax
         x(i,j) = -Lref/2._prec+real(i-1,prec)*dx
         y(i,j) = -Lref/2._prec+real(j-1,prec)*dy
       enddo
      enddo

    write(fid_grid,*) 1
    write(fid_grid,*) imax,jmax,1
    write(fid_grid,'(4e23.15)') ((x(i,j),i=1,imax),j=1,jmax)        !write x nodes
    write(fid_grid,'(4e23.15)') ((y(i,j),i=1,imax),j=1,jmax)        !write y nodes
    write(fid_grid,'(4e23.15)') ((0._prec*y(i,j),i=1,imax),j=1,jmax)!write z nodes (0)
    close(fid_grid)

    else
    stop
    
    endif
    
    
    
  endif
  close(fid_grid)
!*******************************************************************************
!*******************************************************************************




  

  



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
  allocate(djdx(imax,jmax),djdy(imax,jmax),TE(imax,jmax),u(imax,jmax)&
          ,djdxmat(imax,jmax,2),xnew(imax,jmax),ynew(imax,jmax))


!setup grid info file **********************************************************
  inquire(file=trim(grid_info), exist=connected)
  open(fid_gridinfo,file=trim(grid_info),status='unknown')
  if ((append_replace.eq.'replace').or.(connected.eqv..false.)) then
    write(fid_gridinfo,'(A)') 'Title="'//trim(grid_in)//' info"'
    write(fid_gridinfo,'(A)') 'VARIABLES="x""y""u""TE""J""dJ/dx""dJ/dy"'
  
  !write grid history
  TE = truncation_error(x,y,alpha,nu,Lref,imax,jmax)
  Jfunc = calc_J(x,y,nu,alpha,Lref,imax,jmax)
  djdxmat = calc_djdx(x,y,nu,alpha,Lref,imax,jmax)
  djdx = djdxmat(:,:,1)
  djdy = djdxmat(:,:,2)
  u = uxn(x,y,alpha,nu,Lref,0)


  
  write(fid_gridinfo,'(A)') 'zone T="Orginal grid"'
  write(fid_gridinfo,*) 'i=',imax
  write(fid_gridinfo,*) 'j=',jmax
  
  write(fid_gridinfo,'(7e23.15)') ((x(i,j),y(i,j),u(i,j),&
                                    TE(i,j),Jfunc,djdx(i,j),djdy(i,j),&
                                    i=1,imax),j=1,jmax)
  
  elseif (append_replace.eq.'append') then
    io_error = 0
    do while (io_error .eq. 0)
      read(fid_gridinfo,*,iostat=io_error) temp
    enddo
  endif

!*******************************************************************************
!*******************************************************************************

  





!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
  
  !call burgers2_opt(xnew,ynew,x,y,alpha,nu,Lref,imax,jmax,500)

  TE = truncation_error(x,y,alpha,nu,Lref,imax,jmax)
  Jfunc = calc_J(x,y,nu,alpha,Lref,imax,jmax)
  djdxmat = calc_djdx(x,y,nu,alpha,Lref,imax,jmax)
  djdx = djdxmat(:,:,1)
  djdy = djdxmat(:,:,2)
  u = uxn(x,y,alpha,nu,Lref,0)


  !write grid history
  write(fid_gridinfo,'(A)') 'zone T="grid"'
  write(fid_gridinfo,*) 'i=',imax
  write(fid_gridinfo,*) 'j=',jmax
  
  write(fid_gridinfo,'(7e23.15)') ((x(i,j),y(i,j),u(i,j),&
                                    TE(i,j),Jfunc,djdx(i,j),djdy(i,j),&
                                    i=1,imax),j=1,jmax)

  close(fid_gridinfo)
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************  
  
  

                                   
  
  
  
  !write sensitivities to file
  !open(101,file='djdx.dat',status='unknown')
  !write(101,*)imax, jmax
  !write(101,'(4e23.15)') ((x(i,j),y(i,j),djdx(i,j),djdy(i,j),&
  !                                 i=1,imax),j=1,jmax)
                                   
  !close(101)
  
  
  !write functional to file
  !open(101,file='functional.dat',status='unknown')
  !write(101,*) Jfunc
  !close(101)
  
  
  !write grid file
  !open(fid_grid,file=trim(grid_in),status='unknown')
  !write(fid_grid,*) 1
  !write(fid_grid,*) imax, jmax, 1
  !write(fid_grid,'(4e23.15)') ((xnew(i,j),i=1,imax),j=1,jmax)
  !write(fid_grid,'(4e23.15)') ((ynew(i,j),i=1,imax),j=1,jmax)
  !write(fid_grid,'(4e23.15)') ((0._prec*xnew(i,j),i=1,imax),j=1,jmax)    
  
  !close(fid_grid)
  
end program main

