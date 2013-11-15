!module to calculate generic finite difference derivatives for a given vector
! dnfdxn(deriv,f,dx,dorder,order,n)
!   INPUTS:
!          f => vector of size n (real)
!         dx => scalar with grid spacing (real)
!     dorder => derivative order (integer)
!      order => order of accuracy (integer)
!          n => size of function vector (integer)
!
!  OUTPUTS:
!      deriv => vector of derivatives of size n
!
!   ERRORS:
!      If n is less than the required number of points to calculate a derivative
!      then program stop.
!*******************************************************************************
module fd_derivative_calc

contains

subroutine dnfdxn(deriv,f,dx,dorder,order,n)
  use select_precision, only : prec
  use misc_func
! calculate an arbitrary derivative (dorder) of arbitrary order of accuracy (order)
! deriv : derivative vector of size n
! dx: grid spacing
! f : value vector of size n
! dorder : derivative order
! order  : order of accuracy

  implicit none
  integer, intent(in)    :: n      !number of elements in f and deriv
  integer, intent(in)    :: dorder !order of derivative
  integer, intent(in)    :: order  !order of accuracy of derivative
  real(prec), intent(in) :: dx     !grid spacing
  real(prec),dimension(n),intent(in) :: f
  real(prec),dimension(n),intent(out):: deriv

  integer:: ntrunc=2 !number of additional terms to include for order check

  !misc for determining coefficients of the derivatives
  integer:: nterms
  real(prec) :: sten_width
  integer :: i,j, cnt, len_cur_sten
  real(prec), allocatable, dimension(:,:) :: stencil,A,trunc,coef
  logical, allocatable, dimension(:) :: Ine0, Ieq0
  real(prec), allocatable, dimension(:) :: cur_stencil,mat,b,final_order
  real(prec), allocatable, dimension(:) :: final_trunc, trunc_terms,multipliers
  integer, allocatable,dimension(:):: ones,pwr
  
  !misc for calculating the derivatives
  integer :: Icenter, nshift,isten
  integer ,allocatable,dimension(:,:) :: phy_stencil
  real(prec),allocatable,dimension(:)   :: cur_coef



nterms = (dorder-1) + (order-1)+1  !number of terms in stencil


sten_width = real(nterms,prec)/2.0_prec; ! number of nodes in stencil minus 1
if (sten_width+1.0_prec .gt. n) then
    write(*,'(A)')'Error! Not enough points to calculate derivative!'
    stop
endif


!create all possible stencils from backward differences to forward differences
allocate(stencil(nterms+1,nterms+1))
do i=1,nterms+1
  stencil(1,i)=real(i-1,prec)
enddo

do j = 1,nterms 
    stencil(j+1,:) = stencil(j,:)-1.0_prec;
enddo

!solve Taylor Series 
allocate(Ine0(nterms+1),&
         Ieq0(nterms+1),&
         ones(nterms+1),&
         pwr(nterms+ntrunc),&
         mat(nterms+ntrunc),&
         coef(nterms+1,nterms+1),&
         trunc_terms(nterms+1),&
         final_order(nterms+1),&
         final_trunc(nterms+1))

do i = 1,nterms+1
   !det variable defaults
   Ine0=.false.
   Ieq0=.false.
   ones = 0

   !find parts of current stencil that does not equal 0 (i.e. excluding current
   !grid node
   where (stencil(i,:).ne.0)
    Ine0=.true.
    ones=1
   end where

   !get length of the current stencil
   len_cur_sten=sum(ones)

   !extract current stencil from stencil matrix
   allocate(cur_stencil(len_cur_sten),&
            A(len_cur_sten,len_cur_sten),&
            trunc(nterms+ntrunc-len_cur_sten,len_cur_sten),&
            b(len_cur_sten),&
            multipliers(len_cur_sten))

   cnt = 0
   do j = 1,nterms+1
     if (Ine0(j)) then
       cnt=cnt+1
       cur_stencil(cnt)=stencil(i,j)
     endif
   enddo

   !create vector of exponents
   do j = 1,nterms+ntrunc
     pwr(j) = j
   enddo
      
   !calculate TS expansion in matrix form (A) and additional terms for order
   ! checking (trunc)
   do j = 1,len_cur_sten
      mat = cur_stencil(j)**pwr/real(factoriali(pwr),prec)
      A(:,j) = mat(1:len_cur_sten);
      trunc(:,j) = mat(len_cur_sten+1:nterms+ntrunc);
   enddo
   
   !set up b matrix
   b = 0._prec;
   b(dorder) = 1;

   !solve using linsolve.o package
   call gaussian_elimination(multipliers,A,b,len_cur_sten)

   !reconstruct the stencil to include the ith node
   cnt=0
   do j = 1,nterms+1
     if (Ine0(j)) then
       cnt=cnt+1
       coef(i,j)=multipliers(cnt)
     endif
   enddo

   where (Ine0.eqv..false.)
     coef(i,:) = -sum(coef(i,:),Ine0)
   end where

   !reduce truncation error terms using the previously found coefficients
   do j = 1,nterms+ntrunc-len_cur_sten
     trunc_terms(j)=sum(trunc(j,:)*multipliers)
   enddo
   
   !determine first truncation error term that has a non-zero coefficient
   do j = 1,ntrunc
     if (trunc_terms(j).ne.0._prec) then
       final_order(i)=pwr(len_cur_sten+j)-dorder
       final_trunc(i)=trunc_terms(j)
       exit
     endif
   enddo


deallocate(cur_stencil,A,b,trunc,multipliers) 
enddo

!divide by grid spacing
 coef = coef/dx**dorder
 

!Calculate derivatives **********************************************************
allocate(phy_stencil(nterms+1,nterms+1),cur_coef(nterms+1))
Icenter = ceiling(real(nterms+1,prec)/2._prec);

do i = 1,n
   phy_stencil = int(stencil) + i;

   if (minval(phy_stencil(Icenter,:)) .lt. 1) then
       nshift = 1-minval(phy_stencil(Icenter,:))
       Isten = Icenter-nshift
   elseif (maxval(phy_stencil(Icenter,:)) .gt. n) then
       nshift = maxval(phy_stencil(Icenter,:) - n)
       Isten = Icenter + nshift
   else
       Isten = Icenter
   endif
      
   cur_coef = coef(Isten,:)

 
   deriv(i) = sum(cur_coef*f(phy_stencil(Isten,:)))
    
enddo

end subroutine dnfdxn

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


end module fd_derivative_calc
