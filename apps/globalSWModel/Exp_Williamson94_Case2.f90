module Exp_Williamson94_Case2

  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  implicit none
  private

  public :: setIniCond_ExpWS94Case2, callBack_EndCurrentTimeStep
  
  type(volScalarField) :: v_height0

contains
subroutine setIniCond_ExpWS94Case2()
  

  real(DP) :: u0
  real(DP), parameter :: alpha = 0.5*acos(-1d0)
  integer :: faceNum, faceId, cellNum, cellId, lfaceId
  type(Vector3d) :: geo_vel, geoPos, faceNormal
  integer :: maxId(1)
  real(DP) :: r
  type(Vector3d) :: cartPos

  real(DP), parameter :: h0 = 2.94d04 / Grav ! 

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

  do cellId=1, cellNum
     geoPos = cartToSphPos( plMesh%CellPosList(cellId) )
     v_height%data%v_(cellId) =  h0 &
          &                    - (Radius*Omega*u0 + 0.5d0*u0**2) * sin(geoPos%v_(2))**2 / Grav
  end do

  call GeometricField_Init(v_height0, plMesh, "v_height0")
  call DeepCopy( v_height0, v_height )

end subroutine setIniCond_ExpWS94Case2

subroutine callBack_EndCurrentTimeStep(tstep, v_heightN, v_normalVelN)
  integer, intent(in) :: tstep
  type(volScalarField), intent(in) :: v_heightN
  type(surfaceScalarField), intent(in) :: v_normalVelN

  real(DP) :: L2ErrorNorm, LinfErrorNorm

if( mod( tstep*DelTime, 100 ) /= 0 ) return;
 
  L2ErrorNorm = sqrt( sum( (v_heightN%data%v_ - v_height0%data%v_)**2 * fvmInfo%v_CellVol%data%v_ ) ) &
       & / sqrt( sum( v_height0%data%v_**2 * fvmInfo%v_CellVol%data%v_ ) )

  LinfErrorNorm = maxVal( abs(v_heightN%data%v_ - v_height0%data%v_) )/ maxval(v_height0%data%v_)
  
write(*,*) "l2 error norm=", L2ErrorNorm, ": Linf error norm=", LinfErrorNorm

end subroutine callBack_EndCurrentTimeStep

end module Exp_Williamson94_Case2
