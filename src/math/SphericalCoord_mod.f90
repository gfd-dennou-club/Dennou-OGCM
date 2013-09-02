module SphericalCoord_mod
  use dc_types
  use VectorSpace_mod
  implicit none
  private

  interface SphToCartPos
     module procedure SphToCartPos1
     module procedure SphToCartPos2
  end interface SphToCartPos

  interface SphToCartVec
     module procedure SphToCartVec1
     module procedure SphToCartVec2
  end interface SphToCartVec

  public :: SphToCartPos, CartToSphPos
  public :: SphToCartVec
  public :: RadToDegUnit, DegToRadUnit

  public :: geodesicArcLength, sphericalTriArea

  real(DP), save  :: PI = acos(-1d0)

contains
function SphToCartPos1(lambda, phi, r) result(cartPos)
  real(DP), intent(in) :: lambda, phi, r
  type(Vector3d) :: cartPos
  
  cartPos%v_(1:3) = (/ &
       & r*cos(phi)*cos(lambda), r*cos(phi)*sin(lambda), r*sin(phi) /)

end function SphToCartPos1

function SphToCartPos2(lonlatVec) result(cartPos)
  type(vector3d), intent(in) :: lonlatVec
  type(Vector3d) :: cartPos
  
  cartPos = SphToCartPos1(lonlatVec%v_(1), lonlatVec%v_(2), lonlatVec%v_(3))

end function SphToCartPos2


function SphToCartVec1(u, v, w, cartPos) result(cartVec)
  real(DP), intent(in) :: u, v, w
  type(Vector3d), intent(in) :: cartPos
  type(Vector3d) :: cartVec
  
  type(Vector3d) :: sphPos

  sphPos = cartToSphPos(cartPos)
  cartVec%v_(1) = w*cos(sphPos%v_(1))*cos(sphPos%v_(2)) &
       & - v*cos(sphPos%v_(1))*sin(sphPos%v_(2)) - u*sin(sphPos%v_(1))
  cartVec%v_(2) = w*sin(sphPos%v_(1))*cos(sphPos%v_(2)) &
       & - v*sin(sphPos%v_(1))*sin(sphPos%v_(2)) + u*cos(sphPos%v_(3))
  cartVec%v_(3) = w*sin(sphPos%v_(2)) + v*cos(sphPos%v_(2))

end function SphToCartVec1

function SphToCartVec2(sphVec, cartPos) result(cartVec)
  type(vector3d), intent(in) :: sphVec
  type(Vector3d), intent(in) :: cartPos
  type(Vector3d) :: cartVec
  
  cartVec = SphToCartVec1(sphVec%v_(1), sphVec%v_(2), sphVec%v_(3), cartPos)

end function SphToCartVec2

function CartToSphPos(cartPos) result(sphPos)
  type(Vector3d), intent(in) :: cartPos
  type(Vector3d) :: sphPos
  
  real(DP) :: tmp, r

  r = l2norm(cartPos)
  tmp = sqrt(cartPos%v_(1)**2 + cartPos%v_(2)**2)

  if ( cartPos%v_(2) >= 0d0 ) then
     sphPos%v_(:) = (/ acos(cartPos%v_(1)/tmp), atan(cartPos%v_(3)/tmp), r /)
  else
     sphPos%v_(:) = (/ -acos(cartPos%v_(1)/tmp), atan(cartPos%v_(3)/tmp), r /)
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

  area = l2norm(p1)**2 * 4d0 * atan( &
    &  sqrt( abs(tan(s_2)*tan(s_2 - a_2)*tan(s_2 - b_2)*tan(s_2 - c_2) ) )  &
    & )
  
end function sphericalTriArea

function RadToDegUnit(radPos) result(degPos)
  type(Vector3d), intent(in) :: radPos
  type(Vector3d) :: degPos

  degPos%v_(1:2) = (180d0/PI) * radPos%v_(1:2)
  degPos%v_(3) = radPos%v_(3)

end function RadToDegUnit

function DegToRadUnit(degPos) result(radPos)
  type(Vector3d), intent(in) :: degPos
  type(Vector3d) :: radPos

  radPos%v_(1:2) = (PI/180d0) * degPos%v_(1:2)
  radPos%v_(3) = degPos%v_(3)

end function DegToRadUnit

end module SphericalCoord_mod
