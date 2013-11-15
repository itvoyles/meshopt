module functional
  use functional_integral_te_squared
  use select_precision, only : prec
  
  real(prec) :: Jfunctional
  
contains

function calc_J(x,y,v,a,l,imax,jmax)
  use select_precision, only : prec
  use burgers2d_functions,  only : truncation_error	
  implicit none
  
  integer, intent(in) :: imax,jmax
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  real(prec), intent(in) :: v,a,l
  real(prec) :: calc_J
  
  real(prec), dimension(imax,jmax) :: TE
  
  TE = Truncation_error(x,y,a,v,l,imax,jmax)
  calc_J = functional_J(x,y,TE,imax,jmax)



end function calc_J


function calc_djdx(x,y,v,a,L,imax,jmax)
  use select_precision, only : prec
  use burgers2d_functions
  implicit none

  integer, intent(in) :: imax,jmax
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  real(prec), intent(in) :: v,a,L
  real(prec), dimension(imax,jmax,2) :: calc_djdx,calc_djdx2
  integer, dimension(4) :: ex
  integer :: len
  
  real(prec) :: err_lim !error tolerance (normalized)
  
  integer :: j,i,iloc,jloc
  real(prec) :: dx,err
  integer :: il,ih,jl,jh

  Jfunctional = calc_J(x,y,v,a,l,imax,jmax)
  
  err_lim = 0.0001_prec
  
  calc_djdx = 0._prec
  do j = 1,jmax
  do i = 1,imax
  
  dx = (maxval(x(:,j))-minval(x(:,j)))/100000._prec

  !find index for submatrix to calculate djdx (submatrix = 9x9)
  len = 5
    if (i<1+len) then
      il = max(i-len,1)
      ih = il+2*len
      iloc = i
      ex(1) = 0
      ex(2) = 2
    elseif (i>imax-len) then
      ih = min(i+len,imax)
      il = ih-2*len
      iloc = imax-il+1
      ex(1) = 2
      ex(2) = 0
   else
      il = i-len
      ih = i+len
      iloc = 1+len
      ex(1) = 1
      ex(2) = 1
   endif

    if (j<1+len) then
      jl = max(j-len,1)
      jh = jl+2*len
      jloc = j
      ex(3) = 0
      ex(4) = 2
    elseif (j>jmax-len) then
      jh = min(j+len,jmax)
      jl = jh-2*len
      jloc=jmax-jl+1
      ex(3) = 2
      ex(4) = 0
   else 
     jl = j-len
     jh = j+len
     jloc = 1+len
   endif

  if (il<1) then
    il = 1
    ex(1) = 0
    iloc = i
  endif
  if (ih > imax) then
    ih = imax
    ex(2) = 0
    iloc = i
  endif
  if (jh > jmax) then
    jh = jmax
    ex(4) = 0
    jloc = j
  endif
  if (jl<1) then
    jl=1
    ex(3) = 0
    jloc = j
  endif


  !TEloc = TE(il:ih,jl:jh)
  !iloc = i
  !jloc = j
  calc_djdx(i,j,1) = djdx_loc1(x(il:ih,jl:jh),y(il:ih,jl:jh),iloc,jloc,dx,&
                              v,a,L,ih-il+1,jh-jl+1,ex)
  calc_djdx(i,j,2) = djdy_loc1(x(il:ih,jl:jh),y(il:ih,jl:jh),iloc,jloc,dx,&
                              v,a,L,ih-il+1,jh-jl+1,ex)
                              
  !calc_djdx(i,j,1) = djdx_loc(x,y,i,j,dx,&
   !                           err_lim,v,a,L,imax,jmax)
  !calc_djdx(i,j,2) = djdy_loc(x,y,i,j,dx,&
  !                            err_lim,v,a,L,imax,jmax)
                              
    enddo
  enddo



  calc_djdx(1,1:jmax,1) = 0._prec
  calc_djdx(imax,1:jmax,1) = 0._prec

  calc_djdx(1:imax,1,2) = 0._prec
  calc_djdx(1:imax,jmax,2) = 0._prec

  !remove NaN
  where (calc_djdx.ne.calc_djdx)
  calc_djdx = 0._prec
  end where


  !write(*,'(9e12.4)') calc_djdx(1:9,1:9)

end function calc_djdx


function djdx_loc(x,y,ii,jj,dxin,err_lim,v,a,L,imax,jmax,exin)
!given a mesh, this function calculates the mesh sensitivity of the functional
!J by perturbing the mesh by dx. Richardson extrapolation is then used to 
!estimate the error. If the error is below the defined tolerance (err_lim) then
!djdx is returned, if not the required refinement factor is calculated to reach 
!this limit and the derivative is calculated again.
  use select_precision, only : prec
  use burgers2d_functions
  implicit none
  
  integer, intent(in) :: ii,jj,imax,jmax
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  integer, dimension(4), intent(in), optional :: exin
  integer, dimension(4) :: ex
  real(prec), intent(in) ::dxin,err_lim,v,a,L
  real(prec) :: djdx_loc
  
  real(prec) :: err,dx
  real(prec), dimension(imax,jmax) :: xip1,xim1
  real(prec) :: Jip1,Jim1
  real(prec),dimension(2) :: dJdx
  real(prec) :: pf = 2._prec
  real(prec) :: r
  integer :: niter
  integer :: il,ih,jl,jh

  ex = (/ 0,0,0,0 /)
  if (present(exin)) ex = exin

  err = 1._prec
  dx = dxin
  niter = 0



  xip1 = x
  xim1 = x


  do while (abs(err)>err_lim)
  niter = niter + 1

  !calc 1
  xip1(ii,jj) = x(ii,jj)+dx
  Jip1 = Jmod(xip1,y,imax,jmax,a,v,l,ex)

  xim1(ii,jj) = x(ii,jj)-dx
  Jim1 = Jmod(xim1,y,imax,jmax,a,v,l,ex)

  dJdx(1) = (Jip1-Jim1)/(2._prec*dx)


  !calc 2
  dx = dx/2

  xip1(ii,jj) = x(ii,jj)+dx
  Jip1 = Jmod(xip1,y,imax,jmax,a,v,l,ex)

  xim1(ii,jj) = x(ii,jj)-dx
  Jim1 = Jmod(xim1,y,imax,jmax,a,v,l,ex)

  dJdx(2) = (Jip1-Jim1)/(2._prec*dx)

  err = (dJdx(1) - dJdx(2))/(2._prec**pf-1._prec)
  if (dJdx(2)/=0._prec) err = err/dJdx(2)


  !print*, ii,jj,Jip1,Jim1,djdx(1),djdx(2),err,abs(err/err_lim),r,dx
  !print*, (abs(err)>err_lim)


  r = (abs(err/err_lim))**(1._prec/pf)
  dx = dx/r


  enddo

  djdx_loc = djdx(2)
  !print*, ii,jj,niter,abs(err)

end function djdx_loc




function djdy_loc(x,y,ii,jj,dyin,err_lim,v,a,L,imax,jmax,exin)
!given a mesh, this function calculates the mesh sensitivity of the functional
!J by perturbing the mesh by dy. Richardson extrapolation is then used to 
!estimate the error. If the error is below the defined tolerance (err_lim) then
!djdy is returned, if not the required refinement factor is calculated to reach 
!this limit and the derivative is calculated again.
  use select_precision, only : prec
  use burgers2d_functions
  implicit none
  
  integer, intent(in) :: ii,jj,imax,jmax
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  integer, dimension(4), intent(in), optional :: exin
  integer, dimension(4) :: ex
  real(prec), intent(in) ::dyin,err_lim,v,a,L
  real(prec) :: djdy_loc
  
  real(prec) :: err,dy
  real(prec), dimension(imax,jmax) :: yip1,yim1
  real(prec) :: Jip1,Jim1
  real(prec),dimension(2) :: djdy
  real(prec) :: pf = 2._prec
  real(prec) :: r
  integer :: niter


  ex = (/ 0,0,0,0 /)
  if (present(exin)) ex = exin


  err = 1._prec
  dy = dyin
  niter = 0

  yip1 = y
  yim1 = y



  do while (abs(err)>err_lim)
  niter = niter + 1

  !calc 1
  yip1(ii,jj) = y(ii,jj)+dy
  Jip1 = Jmod(x,yip1,imax,jmax,a,v,l,ex)
	
  yim1(ii,jj) = y(ii,jj)-dy
  Jim1 = Jmod(x,yim1,imax,jmax,a,v,l,ex)


  djdy(1) = (Jip1-Jim1)/(2._prec*dy)


  !calc 2
  dy = dy/2._prec

  yip1(ii,jj) = y(ii,jj)+dy
  Jip1 = Jmod(x,yip1,imax,jmax,a,v,l,ex)

  yim1(ii,jj) = y(ii,jj)-dy
  Jim1 = Jmod(x,yim1,imax,jmax,a,v,l,ex)

  djdy(2) = (Jip1-Jim1)/(2._prec*dy)

  err = (djdy(1) - djdy(2))/(2._prec**pf-1._prec)
  if (dJdy(2)/=0._prec) err = err/dJdy(2)
  
  r = (abs(err/err_lim))**(1._prec/pf)
  dy = dy/r

  !print*, Jip1,Jim1,djdy(1),djdy(2),err,abs(err/err_lim),r,dy
  !print*, (abs(err)>err_lim)
  enddo

  djdy_loc = djdy(2)

  !print*, ii,jj,niter,abs(err)

end function djdy_loc

function djdx_loc1(x,y,ii,jj,dx,v,a,L,imax,jmax,exin)
!given a mesh, this function calculates the mesh sensitivity of the functional
!J by perturbing the mesh by dx. Richardson extrapolation is then used to 
!estimate the error. If the error is below the defined tolerance (err_lim) then
!djdx is returned, if not the required refinement factor is calculated to reach 
!this limit and the derivative is calculated again.
  use select_precision, only : prec
  use burgers2d_functions
  implicit none
  
  integer, intent(in) :: ii,jj,imax,jmax
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  integer, dimension(4), intent(in), optional :: exin
  integer, dimension(4) :: ex
  real(prec), intent(in) ::dx,v,a,L
  real(prec) :: djdx_loc1
  

  real(prec), dimension(imax,jmax) :: xip1
  real(prec) :: Jip1,Jim1
  integer :: il,ih,jl,jh

  ex = (/ 0,0,0,0 /)
  if (present(exin)) ex = exin




  xip1 = x
  xip1(ii,jj) = x(ii,jj)+dx
  Jip1 = Jmod(xip1,y,imax,jmax,a,v,l,ex)

  dJdx_loc1 = (Jip1-Jfunctional)/dx

end function djdx_loc1

function djdy_loc1(x,y,ii,jj,dy,v,a,L,imax,jmax,exin)
!given a mesh, this function calculates the mesh sensitivity of the functional
!J by perturbing the mesh by dx. Richardson extrapolation is then used to 
!estimate the error. If the error is below the defined tolerance (err_lim) then
!djdx is returned, if not the required refinement factor is calculated to reach 
!this limit and the derivative is calculated again.
  use select_precision, only : prec
  use burgers2d_functions
  implicit none
  
  integer, intent(in) :: ii,jj,imax,jmax
  real(prec), dimension(imax,jmax), intent(in) :: x,y
  integer, dimension(4), intent(in), optional :: exin
  integer, dimension(4) :: ex
  real(prec), intent(in) ::dy,v,a,L
  real(prec) :: djdy_loc1
  

  real(prec), dimension(imax,jmax) :: yip1
  real(prec) :: Jip1
  integer :: il,ih,jl,jh

  ex = (/ 0,0,0,0 /)
  if (present(exin)) ex = exin



  yip1 = y
  yip1(ii,jj) = y(ii,jj)+dy
  Jip1 = Jmod(x,yip1,imax,jmax,a,v,l,ex)

  dJdy_loc1 = (Jip1-Jfunctional)/dy

end function djdy_loc1



!used for calculating dJ/dx. Excludes boundary. This is not meant to calculate J only the change in J!
!The input matrix is assumed to be 9x9, the maximum stencil that a change in x affects TE plus a boundary
!required to calculate TE.
function Jmod(x,y,imax,jmax,a,v,l,exin)
  use select_precision, only : prec
  use burgers2d_functions
  implicit none

  integer, intent(in) :: imax,jmax
  real(prec),intent(in) :: a, v, l
  integer, intent(in), dimension(4),optional :: exin
  integer, dimension(4) :: ex
  real(prec), dimension(imax,jmax),intent(in) :: x,y
  real(prec)                       :: Jmod
  real(prec), dimension(imax,jmax) :: TEtemp

  ex = (/ 0,0,0,0 /)
  if (present(exin)) ex = exin

   Jmod = functional_J(x(1+ex(1):imax-ex(2),1+ex(3):jmax-ex(4)),&
                       y(1+ex(1):imax-ex(2),1+ex(3):jmax-ex(4)),&
                       Truncation_error(x(1+ex(1):imax-ex(2),1+ex(3):imax-ex(4)),&
                                        y(1+ex(1):imax-ex(2),1+ex(3):imax-ex(4)),&
                                        a,v,l,imax-ex(1)-ex(2),jmax-ex(3)-ex(4)),&
                       imax-ex(1)-ex(2),jmax-ex(3)-ex(4))


end function Jmod

end module functional
