module SphericalCoord_mod
  use dc_types
  use VectorSpace_mod
  implicit none
  private

  public :: SphToCartPos, CartToSphPos
  public :: RadToDegUnit, DegToRadUnit

  public :: geodesicArcLength, sphericalTriArea

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

function geodesicArcLength(p1, p2) result(dist)
  type(Vector3d), intent(in) :: p1, p2
  real(DP) :: dist

  dist = &
    & l2norm(p1) * &
    & 2d0 * asin( 0.5* l2norm(normalizedVec(p1) - normalizedVec(p2)) )

end function geodesicArcLength

function sphericalTriArea(p1, p2, p3) result(area)
  type(Vector3d), intent(in) :: p1, p2, p3
  real(DP) :: area

  real(DP) :: a_2, b_2, c_2, s_2
  type(Vector3d) :: nP1, nP2, nP3

  nP1 = normalizedVec(p1)
  nP2 = normalizedVec(p2)
  nP3 = normalizedVec(p3)

  a_2 = 2d0/2d0 * asin( l2norm(nP1 - nP2) )
  b_2 = 2d0/2d0 * asin( l2norm(nP2 - nP3) )
  c_2 = 2d0/2d0 * asin( l2norm(nP3 - nP1) )
  s_2 = a_2 + b_2 + c_2

  area = l2norm(p1) * 4d0 * atan( &
    &  sqrt( abs(tan(s_2)*tan(s_2 - a_2)*tan(s_2 - b_2)*tan(s_2 - c_2) ) )  &
    & )
  
end function sphericalTriArea

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
