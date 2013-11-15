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


!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
!*******************************************************************************
module common_cnmn1
  use select_precision, only : prec
  implicit none

  real(prec) :: delfun,dabfun,fdch,fdchm,ct,ctmin,ctl,ctlmin, &
                alphax,abobj1,theta,obj,phi

  integer :: ndv,ncon,nside,iprint,nfdg,itmax,itrm,icndir,igoto,nac,info,infog,&
             iter,linobj,nscal,ms1,  &
             n1, n2, n3, n4, n5,nlim

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
  implicit none

  real(prec) :: aobj
  real(prec), dimension(6) :: x
  real(prec), dimension(11):: g

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
  implicit none
  
  real(prec), dimension(6) :: s,vlb,vub,scal,df
  real(prec), dimension(11) :: g1,g2,c,ic
  integer, dimension(11) :: isc
  real(prec), dimension(11,11) :: b
  real(prec), dimension(6,11) :: A


end module data_input



module namelist_conpar
  use common_cnmn1
  use common_variable
  use common_andata
  use data_input
  implicit none
  
  namelist /conpar/ infog, info, nfdg,iprint,ndv,itmax,ncon,nside,&
                    icndir,nscal,fdch,fdchm,ct,ctmin,ctlmin,theta,&
                    phi,delfun,dabfun,linobj,itrm,x,vlb,vub,&
                    n1,n2,n3,n4,n5,alphax,abobj1,ctl,isc,scal

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
  use data_input
  implicit none


  integer :: i
   
  !INITIALIZE
  infog=0
  info=0
  nfdg=0
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
  delfun=0.0_prec
  dabfun=0.0_prec
  linobj=0.0_prec
  itrm=0
  n1=6
  n2=11
  n3=11
  n4=11
  n5=22
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
  write(6,conpar)
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
  implicit none
  
 
  aobj = x(1)**2-5._prec*x(1)+x(2)**2-5._prec*x(2)+2._prec*x(3)**2-21._prec*x(3)&
         +x(4)**2+7.0_prec*x(4)+50._prec
  
  g(1) = x(1)**2+x(1)+x(2)**2-x(2)+x(3)**2+x(3)+x(4)**2-x(4)-8.0_prec
  g(2) = x(1)**2-x(1)+2._prec*x(2)**2+x(3)**2+2._prec*x(4)**2-x(4)-10.0_prec
  g(3) = 2._prec*x(1)**2+2._prec*x(1)+x(2)**2-x(2)+x(3)**2-x(4)-5.0_prec
  


end subroutine
