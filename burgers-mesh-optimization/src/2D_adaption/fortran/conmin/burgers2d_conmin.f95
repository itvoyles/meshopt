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
  real(prec) :: Re    = 8._prec
  real(prec) :: alpha 
  real(prec) :: nu

  integer, parameter :: imax = 9
  integer, parameter :: jmax = 9
  integer :: fid_history = 101
  
  integer, parameter :: nl_const = (imax-2)*(jmax-2)
  integer, parameter :: l_const = imax*(jmax-1)+jmax*(imax-1)

  real(prec) :: cp_limit = 0.0872_prec  !0.0872=5(175)deg
  real(prec) :: lin_const_lim = 1e-7
  real(prec), dimension(imax,jmax) :: TEcart
  real(prec) :: TEcartmax, TEcartL2
  
contains
function mat_to_vec(x)	
  implicit none
  
  real(prec), intent(in), dimension(imax,jmax) :: x
  real(prec), dimension(imax*jmax) :: mat_to_vec
  
  integer :: i,j,cnt
  
    !reshape grid
    cnt = 1
    do j = 1,jmax
    do i = 1,imax
      mat_to_vec(cnt) = x(i,j)
   	 	cnt = cnt + 1
    enddo
 	  enddo
 	  
end function

function vec_to_mat(x)
  implicit none
  
  real(prec), intent(in), dimension(imax*jmax) :: x
  real(prec), dimension(imax,jmax) :: vec_to_mat
  
  integer :: i,j,cnt
  
    !reshape grid
    cnt = 1
    do j = 1,jmax
    do i = 1,imax
      vec_to_mat(i,j) = x(cnt)
   	 	cnt = cnt + 1
    enddo
 	  enddo
 	  
end function

end module burgers2d_data


module constraints
  use select_precision, only : prec
  use burgers2d_data
  implicit none
  
  integer, dimension((imax-2)*(jmax-2)) :: xi,yi,xip1,yip1,xjp1,yjp1
  integer, dimension((imax-2)*(jmax-2)) :: xim1,yim1,xjm1,yjm1

  
contains
subroutine setup_indicies
  implicit none
  
  integer, dimension(imax,jmax) :: ii,jj
  integer :: i,j,cnt,n
  
  n = imax*jmax
  
  cnt = 0
  do j = 1,jmax
  do i = 1,imax
    cnt = cnt + 1
    ii(i,j) = cnt
    jj(i,j) = cnt+n
  enddo
  enddo
  
  
  cnt = 0
  do j = 2,jmax-1
  do i = 2,imax-1
    cnt = cnt+1
    xi(cnt) = ii(i,j)
    yi(cnt) = jj(i,j)
    
    xip1(cnt) = xi(cnt)+1
    yip1(cnt) = yi(cnt)+1
    
    xim1(cnt) = xi(cnt)-1
    yim1(cnt) = yi(cnt)-1
    
    xjp1(cnt) = xi(cnt)+imax
    yjp1(cnt) = yi(cnt)+imax
    
    xjm1(cnt) = xi(cnt)-imax
    yjm1(cnt) = yi(cnt)-imax

  enddo
  enddo
  
  


end subroutine setup_indicies

function c1(X,i)
  use burgers2d_data, only : nl_const
  implicit none
  
  integer, intent(in) :: i
  real(prec), intent(in), dimension(i) :: X
  
  real(prec), dimension(nl_const) :: c1
  !real(prec), dimension(nl_const) :: vnorm
  
  !vnorm = sqrt( (X(xip1)-X(xi))**2+(X(yip1)-X(yi))**2 )&
  !      *sqrt( (X(xjp1)-X(xi))**2+(X(yjp1)-X(yi))**2 )
  
  c1 = -(X(xip1)*X(yjp1)-X(xjp1)*X(yip1)-X(xip1)*X(yi)&
        -X(xi)*X(yjp1)+X(xjp1)*X(yi)+X(xi)*X(yip1))+cp_limit

end function c1

function c2(X,i)
  use burgers2d_data, only : nl_const
  implicit none
  
  integer, intent(in) :: i
  real(prec), intent(in), dimension(i) :: X
  
  real(prec), dimension(nl_const) :: c2
  !real(prec), dimension(nl_const) :: vnorm
  
  !vnorm = sqrt( (X(xim1)-X(xi))**2+(X(yim1)-X(yi))**2 )&
  !       *sqrt( (X(xjm1)-X(xi))**2+(X(yjm1)-X(yi))**2 )
         
  c2 = -(X(xim1)*X(yjm1)-X(xjm1)*X(yim1)-X(xim1)*X(yi)&
        -X(xi)*X(yjm1)+X(xjm1)*X(yi)+X(xi)*X(yim1))+cp_limit

end function c2



end module constraints


module conmin_dim
use burgers2d_data

integer, parameter :: n1=2*jmax*imax+2
integer, parameter :: n2=2*nl_const+l_const+2*n1
integer, parameter :: n3=n2
integer, parameter :: n4=n2
integer, parameter :: n5=2*n4
integer, parameter :: n6=1
integer, parameter :: n7=1
integer, parameter :: n8=1
integer, parameter :: n9=1

end module

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
module common_grad
  use select_precision,only : prec
  use conmin_dim
  implicit none
  
  integer, dimension(n2) :: isc
  integer, dimension(n3) :: ic
  real(prec), dimension(n1) ::df
  real(prec), dimension(n1,n3):: A
  
  common /grad/ isc,ic,df,A

end module common_grad



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
module common_cnmn1
  use select_precision, only : prec
  use conmin_dim
  implicit none

  real(prec) :: delfun,dabfun,fdch,fdchm,ct,ctmin,ctl,ctlmin, &
                alphax,abobj1,theta,obj,phi

  integer :: ndv,ncon,nside,iprint,nfdg,itmax,itrm,icndir,igoto,nac,info,infog,&
             iter,linobj,nscal,  &
             nlim
             


  common /cnmn1/ delfun, dabfun, fdch, fdchm, ct, ctmin, ctl, ctlmin, &
                 alphax, abobj1, theta, obj, ndv, ncon, nside, iprint, &
                 nfdg, nscal, linobj, itmax, itrm, icndir, igoto, nac, &
                 info,infog, iter

end module common_cnmn1

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
module common_variable
  use select_precision, only : prec
  use conmin_dim
  implicit none

  real(prec) :: aobj
  real(prec), dimension(n1) :: DV
  real(prec), dimension(n2) :: g

  common /varable/ aobj, DV, g

end module common_variable

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
module common_andata
  implicit none
  
  integer :: loopcnt
  
end module common_andata


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
module data_input
  use select_precision, only : prec
  use conmin_dim
  implicit none
  
  real(prec), dimension(n1) :: vlb, vub
  real(prec), dimension(n2) :: g2
  real(prec), dimension(n1) :: s
  real(prec), dimension(n3,n3) :: b
  real(prec), dimension(n1) :: scal
  real(prec), dimension(n2) :: g1
  real(prec), dimension(n4) :: c
  integer, dimension(n5) :: ms1
end module data_input


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
module namelist_conpar
  use common_cnmn1
  use common_variable
  use common_andata
  use data_input
  use common_grad
  implicit none
  
  namelist /conpar/ infog, info, nfdg,iprint,ndv,itmax,ncon,nside,&
                    icndir,nscal,fdch,fdchm,ct,ctmin,ctlmin,theta,&
                    phi,delfun,dabfun,linobj,itrm,DV,vlb,vub,&
                    alphax,abobj1,ctl,isc,scal

end module namelist_conpar


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
program burgers2d_mincon
  use select_precision, only : prec
  use burgers2d_data
  use burgers2d_functions
  use springsystem2d_force
  implicit none


  real(prec), dimension(imax,jmax) :: x1,y1
  real(prec), dimension(imax,jmax) :: x2,y2
  real(prec), dimension(imax,jmax) :: TE
  real(prec), dimension(imax*jmax*2) :: Xvec,Xout
  real(prec), dimension(imax*jmax*2) :: Force,Forceout
  real(prec), dimension(imax,jmax,2) :: X
  logical :: newgrid
  
  integer :: i,j,k,cnt,ioerror
  
  nu  = Lref*uref/Re
  alpha = RE/2._prec



    !create cartesian grid
    do i = 1,imax
      x1(i,:) = -Lref/2._prec+real(i-1,prec)/real(imax-1,prec)*Lref
    enddo
    do i = 1,jmax
      y1(:,i) = -Lref/2._prec+real(i-1,prec)/real(jmax-1,prec)*Lref
    enddo
    TEcart = Truncation_error(x1,y1,alpha,nu,Lref,imax,jmax)
    TEcartmax = maxval(abs(TEcart))
    TEcartL2 = sqrt( sum(TEcart**2)/(imax*jmax) )



  !check for existing file
  newgrid = .true.
  open(121,file='test.grd',status='old',iostat=ioerror)
  if (ioerror.eq.0) then
    read(121,*)
    read(121,*) i,j
    if ((i==imax) .and. (j==jmax))then
      read(121,*) ((x1(i,j),i=1,imax),j=1,jmax)
      read(121,*) ((y1(i,j),i=1,imax),j=1,jmax)
      newgrid=.false.
    endif
    
  endif
  close(121)
  
  !setup spring system
  call springsystem_setup(imax,jmax)
  
  
  !setup history file --------------------------------------------
  open(fid_history,file='adaption-history.dat',status='unknown')
  if (newgrid) then !create new adaption history grid  
    write(fid_history,'(A)') 'Title="Adaption history"'
    write(fid_history,'(A)') 'VARIABLES="x""y""TE"'

    !write first evenly space grid data
    write(fid_history,'(A,I6,A)') 'zone T="',I,'"'
    write(fid_history,*) 'i=',imax
    write(fid_history,*) 'j=',jmax
    TE = Truncation_error(x1,y1,alpha,nu,Lref,imax,jmax)
  
    do j = 1,jmax
    do i = 1,imax
      write(fid_history,'(3e23.15)') x1(i,j),y1(i,j),TE(i,j)
    enddo
    enddo

    !start grid other than evenly spaced and write to *.grd file
    !call calc_initial_grid(x1,imax,alpha,nu,lref)
    open(121,file='initial.grd',status='unknown')
    write(121,*) 1
    write(121,*) imax, jmax, 1
    write(121,100) ((x1(i,j),i=1,imax),j=1,jmax)
    write(121,100) ((y1(i,j),i=1,imax),j=1,jmax)
    write(121,100) ((0._prec*x1(i,j),i=1,imax),j=1,jmax)
    close(121)


  else !move to append to old adaption history grid
    do while (ioerror.eq.0) 
      read(fid_history,*,iostat=ioerror)
    enddo
  endif

  !Calculate forces inspring system
  Force=reversespringsystem(x1,y1,imax,jmax)

  !call conmin
  call run_conmin(Forceout,Force,imax*jmax)
  
  X = forwardspringsystem(Forceout,imax,jmax)
  
  !set grid
  x2 = X(:,:,1)
  y2 = X(:,:,2)
 

  open(121,file='test.grd',status='unknown')
  write(121,*) 1
  write(121,*) imax, jmax, 1
  write(121,100) ((x2(i,j),i=1,imax),j=1,jmax)
  write(121,100) ((y2(i,j),i=1,imax),j=1,jmax)
  write(121,100) ((0._prec*x2(i,j),i=1,imax),j=1,jmax)
  close(121)
  100 format (4e23.15)

end program



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
subroutine run_conmin(DVout,DVin,ii)
  use select_precision, only : prec
  use common_cnmn1
  use common_variable
  use common_andata
  use namelist_conpar
  use common_grad
  use data_input
  use conmin_dim
  use burgers2d_data
  use burgers2d_functions
  use constraints
  use springsystem2d_force
  implicit none

  integer, intent(in) :: ii  
  real(prec), intent(in), dimension(2*ii) :: DVin
  real(prec), intent(out), dimension(2*ii) :: DVout
  real(prec), dimension(imax,jmax) :: TE
  real(prec), dimension(imax,jmax,2) :: X
  real(prec) :: TEmax, TEL2
  integer :: i,iii,jjj,n,writecheck,printcheck
   
  n = imax*jmax
   
  !INITIALIZE
  infog=0
  info=0
  nfdg=1
  iprint=1
  ndv=2*imax*jmax
  itmax=1000
  ncon=l_const
  nside=0
  icndir=0
  nscal=100
  fdch=0.01_prec
  fdchm=0.01_prec
  ct=0.001_prec
  ctmin=0.0000001_prec
  ctl=0.000001_prec
  ctlmin=0.0000001_prec
  theta=0.0_prec
  phi=5.0_prec
  delfun=0.000000001_prec
  dabfun=0.000000000001_prec
  linobj=0.0_prec
  itrm=0
  alphax=0.5_prec
  abobj1=0.0_prec



  DV(1:2*n) = DVin
  vlb(1:2*n) = -Lref/2._prec
  vub(1:2*n) = Lref/2._prec
  scal(1:2*n) = 0._prec



  
  do i = 1,l_const
  isc=1
  enddo
  !read the parameters from conmin
  !read(5,conpar) !use default values
  !write(6,conpar)
  nlim=itmax*(ndv+5)
  
  !non-iterative part of analysis
  igoto = 0
  writecheck = 0
  printcheck = 0
  call setup_indicies

  
  !iterative part of analysis
  Do I=1,nlim
    loopcnt=I

    !call optimization routine conmin
    call conmin(DV,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)

    if (igoto.eq.0) loopcnt=-999
    
    !analysis module
    
    call analys
    !print tolerance
    if ((mod(Iter,1)==0).and.(printcheck.ne.iter)) then
    X = forwardspringsystem(DV,imax,jmax)
    TE =  Truncation_error(X(:,:,1),&
                                                  X(:,:,2),&
                                                  alpha,nu,Lref,imax,jmax)
                                                  
    TEmax= maxval(abs(TE))
    TEL2 = sqrt( sum(TE**2)/(imax*jmax))
    
      write(*,'(I5,I5,2e12.4)') I,Iter,TEcartmax/temax,TEcartL2/TEL2
      printcheck=iter
    endif
    
    obj=aobj
    if (igoto.eq.0) goto 1100


  !write adaption history
  if ((mod(iter,1)==0).and.(writecheck.ne.iter)) then
  write(fid_history,'(A,I6,A)') 'zone T="',Iter,'"'
  write(fid_history,*) 'i=',imax
  write(fid_history,*) 'j=',jmax

  X = forwardspringsystem(DV,imax,jmax)
  TE =  Truncation_error(X(:,:,1),&
                         X(:,:,2),&
                         alpha,nu,Lref,imax,jmax)
                                      
  
  do jjj = 1,jmax
  do iii = 1,imax
  write(fid_history,'(3e23.15)') X(iii,jjj,1),X(iii,jjj,2),TE(iii,jjj)
  enddo
  enddo
  
  writecheck=iter
  endif
  !end write adaption history

  enddo
  
  
  1100 continue
 
  write(fid_history,'(A,I6,A)') 'zone T="',Iter,'"'
  write(fid_history,*) 'i=',imax
  write(fid_history,*) 'j=',jmax
  X = forwardspringsystem(DV,imax,jmax)
  TE =  Truncation_error(X(:,:,1),&
                         X(:,:,2),&
                         alpha,nu,Lref,imax,jmax)
  
  do jjj = 1,jmax
  do iii = 1,imax
  write(fid_history,'(3e23.15)') X(iii,jjj,1),X(iii,jjj,2),TE(iii,jjj)
  enddo
  enddo
  
  close(fid_history)
  
  DVout = DV(1:ndv)
  
end subroutine

!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
subroutine analys
  use common_variable
  use common_grad
  use common_cnmn1
  use burgers2d_data
  use functional
  use constraints
  use springsystem2d_force
  implicit none
  
  real(prec) :: eps
  real(prec), dimension(imax,jmax) :: x2,y2
  real(prec), dimension(imax,jmax,2) :: X,dftemp
  integer :: i,j,k,cnt,cnt2
  real(prec) :: Vnorm
  eps = lin_const_lim

  !reshape grid
  X = forwardspringsystem(DV,imax,jmax)
  x2 = X(:,:,1)
  y2 = X(:,:,2)
    
  


  if ((info.ge.1).and.(info.le.2)) then
 
  
    aobj = calc_J(x2,y2,nu,alpha,Lref,imax,jmax)


    !linear consraint values
    cnt = 1
    do j = 1,jmax
    do i = 1,imax-1
      g(cnt) = x2(i,j)-x2(i+1,j)+eps
      cnt = cnt + 1
    enddo
    enddo

    do j = 1,jmax-1
    do i = 1,imax
      g(cnt) = y2(i,j)-y2(i,j+1)+eps
      cnt = cnt + 1
    enddo
    enddo

    !nonlinear constraint values
    g(cnt:cnt+nl_const-1) = 0._prec!c1(X,ndv)
    g(cnt+nl_const:cnt+2*nl_const-1) = 0._prec!c2(X,ndv)


  endif

  if (info.eq.2) then

    !consraint values
    !Set linear constraint gradients
    nac = 0
    cnt = 0
    do i = 1,jmax*(imax-1)
      cnt = cnt + 1
      if (g(i).ge.ct) then
        nac = nac + 1
        ic(nac) = i    
        A(i,nac) = 1._prec
        A(i+1,nac) = -1._prec
      endif
    enddo
    
    do j = 1,imax*(jmax-1)+jmax*(imax-1)
      cnt = cnt+1
      if (g(j).ge.ct) then
        nac = nac + 1
        ic(nac) = j    
        A(j,nac) = 1._prec
        A(j+imax,nac) = -1._prec
      endif
    enddo
    
    !Set non-linear constraint gradients
!--    do i = 1,nl_const
!--      if (g(i).ge.ct) then
!--        nac = nac + 1
!--        ic(nac) = i
!--        vnorm = sqrt( (X(xip1(i))-X(xi(i)))**2+(X(yip1(i))-X(yi(i)))**2 )&
!--               *sqrt( (X(xjp1(i))-X(xi(i)))**2+(X(yjp1(i))-X(yi(i)))**2 )
!--        A(xi(i),nac)   =  ( X(yjp1(i))-X(yip1(i)) )/vnorm
!--        A(yi(i),nac)   =  ( X(xip1(i))-X(xjp1(i)) )/vnorm
!--        A(xip1(i),nac) =  (-X(yjp1(i))+X(yi(i))   )/vnorm
!--        A(yip1(i),nac) =  ( X(xjp1(i))-X(xi(i))   )/vnorm
!--        A(xjp1(i),nac) =  ( X(yip1(i))-X(yi(i))   )/vnorm
!--        A(yjp1(i),nac) =  (-X(xip1(i))+X(xi(i))   )/vnorm  
!--      endif    
!--    enddo
!--    
!--    do i = 1,nl_const
!--      if (g(i).ge.ct) then
!--        nac = nac + 1
!--        ic(nac) = i
!--        
!--        vnorm = sqrt( (X(xim1(i))-X(xi(i)))**2+(X(yim1(i))-X(yi(i)))**2 )&
!--               *sqrt( (X(xjm1(i))-X(xi(i)))**2+(X(yjm1(i))-X(yi(i)))**2 )
!--               
!--        A(xi(i),nac)   =  ( X(yjm1(i))-X(yim1(i)) )/vnorm
!--        A(yi(i),nac)   =  ( X(xim1(i))-X(xjm1(i)) )/vnorm
!--        A(xim1(i),nac) =  (-X(yjm1(i))+X(yi(i))   )/vnorm
!--        A(yim1(i),nac) =  ( X(xjm1(i))-X(xi(i))   )/vnorm
!--        A(xjm1(i),nac) =  ( X(yim1(i))-X(yi(i))   )/vnorm
!--        A(yjm1(i),nac) =  (-X(xim1(i))+X(xi(i))   )/vnorm     
!--      endif    
!--    enddo

    !gradients
    dftemp =  calc_dJdx(x2,y2,nu,alpha,Lref,imax,jmax)


    df(1:imax*jmax) = mat_to_vec(dftemp(:,:,1))
    df(imax*jmax+1:imax*jmax*2) = mat_to_vec(dftemp(:,:,2))


  endif




end subroutine




