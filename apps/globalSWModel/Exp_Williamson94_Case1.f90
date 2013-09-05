module Exp_Williamson94_Case1

  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  implicit none
  private

  public :: setIniCond_ExpWS94Case1

contains
subroutine setIniCond_ExpWS94Case1()
  

  real(DP) :: u0
  real(DP), parameter :: alpha = 0.5*acos(-1d0)
  integer :: faceNum, faceId, cellNum, cellId, lfaceId
  type(Vector3d) :: geo_vel, geoPos, faceNormal
  integer :: maxId(1)
  real(DP) :: r
  type(Vector3d) :: centerPos

  real(DP), parameter :: cosBellRadius = radius/3d0
  real(DP), parameter :: h0 = 1000d0

  u0 = 2d0*PI*radius / (3600d0*24d0*12d0)
write(*,*) "U0=", u0

  cellNum = getCellListSize(plMesh)
  faceNum = getFaceListSize(plMesh)

  do faceId=1, faceNum
     geoPos = cartToSphPos( fvmInfo%s_faceCenter.At.faceId )
     geo_vel = (/ u0*(cos(geoPos%v_(2))), 0d0, 0d0 /)
     faceNormal = normalizedVec( fvmInfo%s_faceAreaVec.At.faceId )
     s_normalVel%data%v_(faceId) = SphToCartVec( geo_vel, fvmInfo%s_faceCenter.At.faceId) .dot. faceNormal

  end do

  centerPos = SphToCartPos(3d0*PI/2d0, 0d0, radius)
  do cellId=1, cellNum
     r = geodesicArcLength(centerPos, plmesh%cellPosList(cellId))
     if(r < cosBellRadius) then
        v_height%data%v_(cellId) = 0.5d0*h0*(1d0 + cos(PI*r/cosBellRadius))
     else
        v_height%data%v_(cellId) = 0d0
     end if
  end do

end subroutine setIniCond_ExpWS94Case1

end module Exp_Williamson94_Case1
