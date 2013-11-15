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

module conmin_dim

integer, parameter :: n1=6
integer, parameter :: n2=11
integer, parameter :: n3=11
integer, parameter :: n4=11
integer, parameter :: n5=22
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
  real(prec), dimension(n1) :: x
  real(prec), dimension(n2) :: g

  common /varable/ aobj, x, g

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
                    phi,delfun,dabfun,linobj,itrm,x,vlb,vub,&
                    alphax,abobj1,ctl,isc,scal

end module namelist_conpar



!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
program EXAMPLE1
  use select_precision, only : prec
  use common_cnmn1
  use common_variable
  use common_andata
  use namelist_conpar
  use common_grad
  use data_input
  use conmin_dim
  implicit none


  integer :: i
   
  !INITIALIZE
  infog=0
  info=0
  nfdg=1
  iprint=2
  ndv=4
  itmax=40
  ncon=3
  nside=0
  icndir=0
  nscal=0
  fdch=0.0_prec
  fdchm=0.0_prec
  ct=0.0_prec
  ctmin=0.0_prec
  ctl=0.0_prec
  ctlmin=0.0_prec
  theta=0.0_prec
  phi=0.0_prec
  delfun=0.001_prec
  dabfun=0.0000000001_prec
  linobj=0.0_prec
  itrm=0
  alphax=0.0_prec
  abobj1=0.0_prec



  do i = 1,ndv
  x(i) = 1.0_prec
  vlb(i) = -9999._prec
  vub(i) = 9999._prec
  scal(i) = 0._prec
  enddo

  
  do i = 1,ncon
  isc=0
  enddo
  !read the parameters from conmin
  !read(5,conpar) !use default values
  !write(6,conpar)
  nlim=itmax*(ndv+5)
  
  !non-iterative part of analysis
  igoto = 0

  !iterative part of analysis
  Do I=1,nlim
    loopcnt=I

    !call optimization routine conmin
    call conmin(x,vlb,vub,g,scal,df,a,s,g1,g2,b,c,isc,ic,ms1,n1,n2,n3,n4,n5)

    if (igoto.eq.0) loopcnt=-999
    
    !analysis module
    
    call analys
    obj=aobj
    if (igoto.eq.0) goto 1100

  enddo
  
  1100 continue
  
  
end program

subroutine analys
  use common_variable
  use common_grad
  use common_cnmn1
  implicit none
  
  if (info.lt.2) then
 
  aobj = x(1)**2-5._prec*x(1)+x(2)**2-5._prec*x(2)+2._prec*x(3)**2-21._prec*x(3)&
         +x(4)**2+7.0_prec*x(4)+50._prec
  
  g(1) = x(1)**2+x(1)+x(2)**2-x(2)+x(3)**2+x(3)+x(4)**2-x(4)-8.0_prec
  g(2) = x(1)**2-x(1)+2._prec*x(2)**2+x(3)**2+2._prec*x(4)**2-x(4)-10.0_prec
  g(3) = 2._prec*x(1)**2+2._prec*x(1)+x(2)**2-x(2)+x(3)**2-x(4)-5.0_prec
  

  else
  
  !gradient information

  df(1) = 2.0_prec*x(1)-5.0_prec
  df(2) = 2.0_prec*x(2)-5.0_prec
  df(3) = 4.0_prec*x(3)-21._prec
  df(4) = 2.0_prec*x(4)+7.0_prec

  !gradients of active and violated contraints
  nac = 0
  
  if (g(1).ge.ct) then
    nac = 1
    ic = 1
    A(1,1) = 2._prec*x(1)+1._prec
    A(2,1) = 2._prec*x(2)-1._prec
    A(3,1) = 2._prec*x(3)+1._prec
    A(4,1) = 2._prec*x(4)-1._prec
  endif
  
  if (g(2).ge.ct) then
    nac = nac+1
    ic(nac) = 2
    A(1,nac) = 2._prec*x(1)-1.0_prec
    A(2,nac) = 4._prec*x(2)
    A(3,nac) = 2._prec*x(3)
    A(4,nac) = 4._prec*x(4)-1.0_prec
  endif
  
  if (g(3).ge.ct) then
    nac = nac+1
    ic(nac) = 3
    A(1,nac) = 4._prec*x(1)+2._prec
    A(2,nac) = 2._prec*x(2)-1._prec
    A(3,nac) = 2._prec*x(3)
    A(4,nac) = -1._prec
  endif


  endif


end subroutine
