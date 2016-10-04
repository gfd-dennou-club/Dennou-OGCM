module Exp_Williamson94_Case5

  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod

  implicit none
  private

  public :: setIniCond_ExpWS94Case5, callBack_EndCurrentTimeStep
  
  type(volScalarField), save :: v_h0

contains
subroutine setIniCond_ExpWS94Case5()
  

  real(DP) :: u0
  real(DP), parameter :: alpha = 0.5*acos(-1d0)
  integer :: faceNum, faceId, cellNum, cellId, lfaceId
  type(Vector3d) :: geo_vel, geoPos, geoPos1, geoPos2, faceNormal
  integer :: maxId(1)
  real(DP) :: r
  type(Vector3d) :: cartPos

  real(DP), parameter :: h0 = 5.96d04 / Grav ! 
  real(DP), parameter :: b0 = 2000d0
  real(DP), parameter :: R0 = PI/9d0
  real(DP), parameter :: lambdaC = -PI/2d0
  real(DP), parameter :: thetaC = PI/6d0
 
  u0 = 2d0*PI*radius / (3600d0*24d0*12d0)
write(*,*) "U0=", u0

  cellNum = getCellListSize(plMesh)
  faceNum = getFaceListSize(plMesh)

  do faceId=1, faceNum
     geoPos = cartToSphPos( At(fvmInfo%s_faceCenter,faceId) )
     geo_vel = (/ u0*(cos(geoPos%v_(2))), 0d0, 0d0 /)
     faceNormal = normalizedVec( At(fvmInfo%s_faceAreaVec,faceId) )
     s_normalVel%data%v_(1,faceId) = SphToCartVec( geo_vel, At(fvmInfo%s_faceCenter,faceId) ).dot.faceNormal

  end do

  do cellId=1, cellNum
     geoPos = cartToSphPos( plMesh%CellPosList(cellId) )

     r = sqrt( min(R0**2, (geoPos%v_(1)-lambdaC)**2 + (geoPos%v_(2)-thetaC)**2) )
     v_hb%data%v_(1,cellId) = b0*(1d0 - r/R0)

     v_h%data%v_(1,cellId) =  h0 &
          &                    - (Radius*Omega*u0 + 0.5d0*u0**2) * sin(geoPos%v_(2))**2 / Grav &
          &                    - At(v_hb,CellId)

  end do

  call GeometricField_Init(v_h0, plMesh, "v_h0")
  call DeepCopy( v_h0, v_h )

end subroutine setIniCond_ExpWS94Case5

subroutine callBack_EndCurrentTimeStep(tstep, v_hN, v_normalVelN)
  integer, intent(in) :: tstep
  type(volScalarField), intent(in) :: v_hN
  type(surfaceScalarField), intent(in) :: v_normalVelN

  real(DP) :: L2ErrorNorm, LinfErrorNorm


end subroutine callBack_EndCurrentTimeStep

end module Exp_Williamson94_Case5
