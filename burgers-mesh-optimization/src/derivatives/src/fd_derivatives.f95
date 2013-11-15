!Module of finite difference derivatives
module fd_derivatives
  use select_precision, only : prec


contains
  elemental function fd_dudx(uim1,uip1,dx)
    implicit none
    real(prec), intent(in) :: uim1,uip1,dx
    real(prec)             :: fd_dudx
    fd_dudx = (uip1-uim1)/(2._prec*dx)
  end function

  elemental function fd_d2udx2(uim1,ui,uip1,dx)
    implicit none
    real(prec), intent(in) :: uim1,ui,uip1,dx
    real(prec)             :: fd_d2udx2
    fd_d2udx2 = (uip1-2._prec*ui+uim1)/(dx**2)
  end function

  elemental function fd_d3udx3(uim2,uim1,uip1,uip2,dx)
    implicit none
    real(prec), intent(in) :: uim2,uim1,uip1,uip2,dx
    real(prec)             :: fd_d3udx3
    fd_d3udx3 = (uip2-2._prec*uip1+2._prec*uim1-uim2)/(2._prec*dx**3)
  end function

  elemental function fd_d4udx4(uim2,uim1,ui,uip1,uip2,dx)
    implicit none
    real(prec), intent(in) :: uim2,uim1,ui,uip1,uip2,dx
    real(prec)             :: fd_d4udx4
    fd_d4udx4 = (uip2-4._prec*uip1+6._prec*ui-4._prec*uim1+uim2)/(dx**4)
  end function

  elemental function fd_d2udxy(uip1jp1,uip1jm1,uim1jp1,uim1jm1,dxi,deta)
    implicit none
    real(prec), intent(in) :: uip1jp1,uip1jm1,uim1jp1,uim1jm1,dxi,deta
    real(prec)             :: fd_d2udxy

    fd_d2udxy = (uip1jp1-uip1jm1-uim1jp1+uim1jm1)/(4._prec*dxi*deta)


  end function





end module fd_derivatives
