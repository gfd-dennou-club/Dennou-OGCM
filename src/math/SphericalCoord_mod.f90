module SphericalCoord_mod
  use dc_types
  use VectorSpace_mod
  implicit none
  private

  public :: SphToCartPos, CartToSphPos
  public :: RadToDegUnit, DegToRadUnit

  real(DP), save  :: PI = acos(-1d0)

contains
function SphToCartPos(lambda, phi) result(cartPos)
  real(DP), intent(in) :: lambda, phi
  type(Vector3d) :: cartPos
  
  cartPos%v_(1:3) = (/ &
       & cos(phi)*cos(lambda), cos(phi)*sin(lambda), sin(phi) /)

end function SphToCartPos

function CartToSphPos(cartPos) result(sphPos)
  type(Vector3d), intent(in) :: cartPos
  type(Vector2d) :: sphPos
  
  real(DP) :: tmp

  tmp = sqrt(cartPos%v_(1)**2 + cartPos%v_(2)**2)

  if ( cartPos%v_(2) >= 0d0 ) then
     sphPos%v_(:) = (/ acos(cartPos%v_(1)/tmp), atan(cartPos%v_(3)/tmp) /)
  else
     sphPos%v_(:) = (/ -acos(cartPos%v_(1)/tmp), atan(cartPos%v_(3)/tmp) /)
  end if

end function CartToSphPos

function RadToDegUnit(radPos) result(degPos)
  type(Vector2d), intent(in) :: radPos
  type(Vector2d) :: degPos

  degPos = (180d0/PI) * radPos

end function RadToDegUnit

function DegToRadUnit(degPos) result(radPos)
  type(Vector2d), intent(in) :: degPos
  type(Vector2d) :: radPos

  radPos = (PI/180d0) * degPos

end function DegToRadUnit

end module SphericalCoord_mod
