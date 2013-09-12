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
  
  type(volScalarField) :: v_h0

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
     geoPos = cartToSphPos( fvmInfo%s_faceCenter.At.faceId )
     !geoPos1 = cartToSphPos( plMesh%PointList(fvmInfo%Face_PointId(1,faceId)) )
     !geoPos2 = cartToSphPos( plMesh%PointList(fvmInfo%Face_PointId(2,faceId)) )
     !s_normalVel%data%v_(faceId) = - Radius*u0*( sin(geoPos1%v_(2)) - sin(geoPos2%v_(2)) )/ l2norm(fvmInfo%s_faceAreaVec.At.FaceId)
     geo_vel = (/ u0*(cos(geoPos%v_(2))), 0d0, 0d0 /)
     faceNormal = normalizedVec( fvmInfo%s_faceAreaVec.At.faceId )
     s_normalVel%data%v_(faceId) = SphToCartVec( geo_vel, fvmInfo%s_faceCenter.At.faceId) .dot. faceNormal

  end do

  do cellId=1, cellNum
     geoPos = cartToSphPos( plMesh%CellPosList(cellId) )

     r = sqrt( min(R0**2, (geoPos%v_(1)-lambdaC)**2 + (geoPos%v_(2)-thetaC)**2) )
     v_hb%data%v_(cellId) = b0*(1d0 - r/R0)

     v_h%data%v_(cellId) =  h0 &
          &                    - (Radius*Omega*u0 + 0.5d0*u0**2) * sin(geoPos%v_(2))**2 / Grav &
          &                    - (v_hb.At.CellId)

  end do

  call GeometricField_Init(v_h0, plMesh, "v_h0")
  call DeepCopy( v_h0, v_h )

end subroutine setIniCond_ExpWS94Case5

subroutine callBack_EndCurrentTimeStep(tstep, v_hN, v_normalVelN)
  integer, intent(in) :: tstep
  type(volScalarField), intent(in) :: v_hN
  type(surfaceScalarField), intent(in) :: v_normalVelN

  real(DP) :: L2ErrorNorm, LinfErrorNorm

!!$if( mod( tstep*DelTime, 5*DelTime ) /= 0 ) return;
!!$ 
!!$  L2ErrorNorm = sqrt( sum( (v_hN%data%v_ - v_h0%data%v_)**2 * fvmInfo%v_CellVol%data%v_ ) ) &
!!$       & / sqrt( sum( v_h0%data%v_**2 * fvmInfo%v_CellVol%data%v_ ) )
!!$
!!$  LinfErrorNorm = maxVal( abs(v_hN%data%v_ - v_h0%data%v_) )/ maxval(v_h0%data%v_)
!!$  
!!$write(*,*) "l2 error norm=", L2ErrorNorm, ": Linf error norm=", LinfErrorNorm

!if(mod( tstep*DelTime, outputIntrVal) == 0 ) call OutputData(tstep)

end subroutine callBack_EndCurrentTimeStep

subroutine OutputData( tstep )

  use vtkDataWriter_mod
  use netcdfDataWriter_mod
  use dc_string
  use SphericalCoord_mod

  integer, intent(in) :: tstep
  character(STRING) :: dataFileName
  type(vtkDataWriter) :: writer

type(PointScalarField) :: p_zeta
type(volScalarField) :: v_hError

call GeometricField_Init(v_hError, plMesh, "v_hError")

!call GeometricField_Init(p_zeta, plMesh, "p_zeta", "relative vorticity", "s-1")
!p_zeta = curl(s_normalVel)
!v_div = div(s_normalVel)
v_hError = v_h - v_h0

  dataFileName = CPrintf("errorcheck-%08d.vtk", i=(/ tstep /))
  call vtkDataWriter_Init(writer, dataFileName, plMesh)
  call vtkDataWriter_Regist(writer, (/ v_hError /))
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)

!call GeometricField_Final(p_zeta)
call GeometricField_Final(v_hError)

end subroutine OutputData


end module Exp_Williamson94_Case5
