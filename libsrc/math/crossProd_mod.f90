module crossProd_mod

  use dc_types
  use Vector2d_mod
  use vector3d_mod

  implicit none
  private
  
  interface operator(.cross.)
     module procedure crossProd_vector2d
     module procedure crossProd_vector3d
  end interface operator(.cross.)

  public operator(.cross.)

contains
function crossProd_vector3d(v1, v2) result(cross)
  type(vector3d), intent(in) :: v1, v2
  type(vector3d) :: cross
  
  cross%v_(1) = v1%v_(2)*v2%v_(3) - v1%v_(3)*v2%v_(2)
  cross%v_(2) = v1%v_(3)*v2%v_(1) - v1%v_(1)*v2%v_(3)
  cross%v_(3) = v1%v_(1)*v2%v_(2) - v1%v_(2)*v2%v_(1)
 
end function crossProd_vector3d

function crossProd_vector2d(v1, v2) result(cross)
  type(vector2d), intent(in) :: v1, v2
  type(vector3d) :: cross
  

  cross%v_(1) = 0d0
  cross%v_(2) = 0d0
  cross%v_(3) = v1%v_(1)*v2%v_(2) - v1%v_(2)*v2%v_(1)
 
end function crossProd_vector2d

end module crossProd_mod
