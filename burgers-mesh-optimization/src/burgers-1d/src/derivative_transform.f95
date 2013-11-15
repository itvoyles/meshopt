!Module is used to transform solution derivatives from computational space to
!physcial space
module derivative_transform
  use select_precision, only : prec
  implicit none

contains
!*******************************************************************************
  elemental function dudx_from_dudxi(dudxi,dxdxi)
    implicit none
    real(prec), intent(in) :: dudxi,dxdxi
    real(prec) :: dudx_from_dudxi

    dudx_from_dudxi = dudxi/dxdxi
  end function dudx_from_dudxi
!*******************************************************************************
  elemental function d2udx2_from_d2udxi2(dudxi,d2udxi2,dxdxi,d2xdxi2)
    implicit none
    real(prec), intent(in) :: dudxi,d2udxi2,dxdxi,d2xdxi2
    real(prec) :: d2udx2_from_d2udxi2

    d2udx2_from_d2udxi2 = d2udxi2/dxdxi**2-dudxi*(d2xdxi2/dxdxi**3)
  end function d2udx2_from_d2udxi2
!*******************************************************************************
  elemental function d3udx3_from_d3udxi3(dudxi,d2udxi2,d3udxi3,&
                                         dxdxi,d2xdxi2,d3xdxi3)
    implicit none
    real(prec), intent(in) :: dudxi,d2udxi2,d3udxi3,dxdxi,d2xdxi2,d3xdxi3
    real(prec) :: d3udx3_from_d3udxi3
    real(prec) :: dudx,d2udx2
    
    dudx   =dudx_from_dudxi(dudxi,dxdxi)
    d2udx2 =d2udx2_from_d2udxi2(dudxi,d2udxi2,dxdxi,d2xdxi2)

    d3udx3_from_d3udxi3 =   d3udxi3/dxdxi**3                                   &
                          - 3._prec*d2udx2*d2xdxi2/dxdxi                       &
                          - dudx*d3xdxi3/dxdxi**3

  end function d3udx3_from_d3udxi3

!*******************************************************************************
  elemental function d4udx4_from_d4udxi4(dudxi,d2udxi2,d3udxi3,d4udxi4,        &
                                         dxdxi,d2xdxi2,d3xdxi3,d4xdxi4)
    implicit none
    real(prec), intent(in) :: dudxi,d2udxi2,d3udxi3,d4udxi4
    real(prec), intent(in) :: dxdxi,d2xdxi2,d3xdxi3,d4xdxi4
    real(prec) :: d4udx4_from_d4udxi4
    real(prec) :: dudx,d2udx2,d3udx3
    
    dudx   =dudx_from_dudxi(dudxi,dxdxi)
    d2udx2 =d2udx2_from_d2udxi2(dudxi,d2udxi2,dxdxi,d2xdxi2)
    d3udx3 =d3udx3_from_d3udxi3(dudxi,d2udxi2,d3udxi3,dxdxi,d2xdxi2,d3xdxi3)

    d4udx4_from_d4udxi4 =   d4udxi4/dxdxi**4                                   &
                          - 6._prec*d3udx3*d2xdxi2/dxdxi**2                    &
                          - 3._prec*d2udx2*d2xdxi2**2/dxdxi**4                 &
                          - 4._prec*d2udx2*d3xdxi3/dxdxi**3                    &
                          - dudx*d4xdxi4/dxdxi**4

  end function d4udx4_from_d4udxi4
!*******************************************************************************
end module derivative_transform
