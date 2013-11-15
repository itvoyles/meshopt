!Module is used to transform solution derivatives from computational space to
!physcial space
module derivative_transform
  use select_precision, only : prec
  implicit none

contains
!*******************************************************************************
  elemental function dudx_from_dudxi(xi_x,eta_x,u_xi,u_eta)
    implicit none
    real(prec), intent(in) :: xi_x,u_xi,eta_x,u_eta
    real(prec) :: dudx_from_dudxi

    dudx_from_dudxi = xi_x*u_xi+eta_x*u_eta
  end function dudx_from_dudxi
!*******************************************************************************
  elemental function d2udx2_from_d2udxi2(xi_x,eta_x,xi_xx,eta_xx,&
                                         u_xi,u_eta,u_xixi,u_etaeta,u_xieta)
    implicit none
    real(prec), intent(in) :: xi_x,eta_x,xi_xx,eta_xx,&
                              u_xi,u_eta,u_xixi,u_xieta,u_etaeta
    real(prec) :: d2udx2_from_d2udxi2

    d2udx2_from_d2udxi2 = xi_x**2*u_xixi+2*xi_x*eta_x*u_xieta+eta_x**2*u_etaeta&
                        + u_xi*xi_xx+u_eta*eta_xx
    
  end function d2udx2_from_d2udxi2
!*******************************************************************************
  elemental function dudy_from_dudeta(xi_y,eta_y,u_xi,u_eta)
    implicit none
    real(prec), intent(in) :: xi_y,u_xi,eta_y,u_eta
    real(prec) :: dudy_from_dudeta

    dudy_from_dudeta = xi_y*u_xi+eta_y*u_eta
    
  end function dudy_from_dudeta

!*******************************************************************************
  elemental function d2udy2_from_d2udeta2(xi_y,eta_y,xi_yy,eta_yy,&
                                          u_xi,u_eta,u_xixi,u_etaeta,u_xieta)
    implicit none
    real(prec), intent(in) :: xi_y,eta_y,xi_yy,eta_yy,&
                              u_xi,u_eta,u_xixi,u_xieta,u_etaeta
    real(prec) :: d2udy2_from_d2udeta2

    d2udy2_from_d2udeta2 = xi_y**2*u_xixi+2*xi_y*eta_y*u_xieta+eta_y**2*u_etaeta&
                        + u_xi*xi_yy+u_eta*eta_yy
    
  end function d2udy2_from_d2udeta2  
  
end module derivative_transform
