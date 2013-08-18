module crossProd_mod
  use dc_types
  use vector3d_mod
  implicit none
  private
  
  interface operator(.cross.)
     module procedure crossProd
  end interface operator(.cross.)

  public operator(.cross.)

contains
function crossProd(v1, v2) result(cross)
  type(vector3d), intent(in) :: v1, v2
  type(vector3d) :: cross
  
  cross%v_(1) = v1%v_(2)*v2%v_(3) - v1%v_(3)*v2%v_(2)
  cross%v_(2) = v1%v_(3)*v2%v_(1) - v1%v_(1)*v2%v_(3)
  cross%v_(3) = v1%v_(1)*v2%v_(2) - v1%v_(2)*v2%v_(1)
 
end function crossProd

end module crossProd_mod
