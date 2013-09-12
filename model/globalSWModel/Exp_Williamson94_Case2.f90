module Exp_Williamson94_Case2

  use VectorSpace_mod
  use SphericalCoord_mod
  use GeometricField_mod

  use SimParameters_mod
  use GridSet_mod
  use VariableSet_mod, only: &
       & v_h, s_normalVel

  use OutputData_mod, only: &
       & OutputDataAnalysisInfo


  implicit none
  private

  public :: setIniCond_ExpWS94Case2, callBack_EndCurrentTimeStep
  
  type(volScalarField) :: v_h0

contains
subroutine setIniCond_ExpWS94Case2()
  

  real(DP) :: u0
  real(DP), parameter :: alpha = 0.5*acos(-1d0)
  integer :: faceNum, faceId, cellNum, cellId, lfaceId
  type(Vector3d) :: geo_vel, geoPos, geoPos1, geoPos2, faceNormal
  integer :: maxId(1)
  real(DP) :: r
  type(Vector3d) :: cartPos

  real(DP), parameter :: h0 = 2.94d04 / Grav ! 

  !call CheckGridQuality()

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
     v_h%data%v_(cellId) =  h0 &
          &                    - (Radius*Omega*u0 + 0.5d0*u0**2) * sin(geoPos%v_(2))**2 / Grav
  end do

  call GeometricField_Init(v_h0, plMesh, "v_h0")
  call DeepCopy( v_h0, v_h )

end subroutine setIniCond_ExpWS94Case2

subroutine callBack_EndCurrentTimeStep(tstep, v_hN, v_normalVelN)
  integer, intent(in) :: tstep
  type(volScalarField), intent(in) :: v_hN
  type(surfaceScalarField), intent(in) :: v_normalVelN

  real(DP) :: L2ErrorNorm, LinfErrorNorm

  if( mod( tstep*DelTime, outputIntrVal ) /= 0 ) return;
 
  L2ErrorNorm = sqrt( sum( (v_hN%data%v_ - v_h0%data%v_)**2 * fvmInfo%v_CellVol%data%v_ ) ) &
       & / sqrt( sum( v_h0%data%v_**2 * fvmInfo%v_CellVol%data%v_ ) )

  LinfErrorNorm = maxVal( abs(v_hN%data%v_ - v_h0%data%v_) )/ maxval(v_h0%data%v_)
  
  call MessageNotify("M", "globalSWM", "l2 error norm=%f : Linf error norm=%f", d=(/ L2ErrorNorm, LinfErrorNorm /) )
  call OutputDataAnalysisInfo(tstep, L2ErrorNorm, LinfErrorNorm)

end subroutine callBack_EndCurrentTimeStep
!!$
!!$subroutine OutputData( tstep )
!!$
!!$  use vtkDataWriter_mod
!!$  use netcdfDataWriter_mod
!!$  use dc_string
!!$  use SphericalCoord_mod
!!$
!!$  integer, intent(in) :: tstep
!!$  character(STRING) :: dataFileName
!!$  type(vtkDataWriter) :: writer
!!$
!!$type(PointScalarField) :: p_zeta
!!$type(volScalarField) :: v_hError
!!$
!!$call GeometricField_Init(v_hError, plMesh, "v_hError")
!!$
!!$!call GeometricField_Init(p_zeta, plMesh, "p_zeta", "relative vorticity", "s-1")
!!$!p_zeta = curl(s_normalVel)
!!$!v_div = div(s_normalVel)
!!$v_hError = v_h - v_h0
!!$
!!$  dataFileName = CPrintf("errorcheck-%08d.vtk", i=(/ tstep /))
!!$  call vtkDataWriter_Init(writer, dataFileName, plMesh)
!!$  call vtkDataWriter_Regist(writer, (/ v_hError /))
!!$  call vtkDataWriter_write(writer)
!!$  call vtkDataWriter_Final(writer)
!!$
!!$!call GeometricField_Final(p_zeta)
!!$call GeometricField_Final(v_hError)
!!$
!!$end subroutine OutputData


subroutine checkGridQuality()

  use SphericalCoord_mod

  integer :: cellNum, faceNum, ptNum, i, j
  integer :: lfaceNum, neighCellId, faceId, ptId
  
  real(DP) :: dists(6)
  real(DP) :: sigMin, sigAvg, qMin, qAvg
  type(volScalarField) :: sigma
  type(pointScalarField) :: q

  cellNum = getCellListSize(plMesh)
  faceNum = getFaceListSize(plMesh)
  ptNUm = getPointListSize(plMesh)

  call GeometricField_Init(sigma, plMesh, "sigma")
  call GeometricField_Init(q, plMesh, "q")

  do i=1, cellNum
     lfaceNum = plMesh%cellList(i)%faceNum

     do j=1, lfaceNum
        faceId = fvmInfo%Cell_FaceId(j,i)
        if( fvmInfo%Face_CellId(1,faceId) == i) then
           neighCellId = fvmInfo%Face_CellId(2,faceId)
        else
           neighCellId = fvmInfo%Face_CellId(1,faceId)
        end if
!        dists(j) = geodesicArcLength( plMesh%CellPosList(i), plMesh%CellPosList(neighCellId) )
        dists(j) = l2norm(plMesh%CellPosList(i) - plMesh%CellPosList(neighCellId) )
     end do

     sigma%data%v_(i) = minval(dists(1:lfaceNum)) / maxval(dists(1:lfaceNum))
!write(*,*) i, sum(dists(1:lfaceNum))/lfaceNum
  end do

  do ptId=1, ptNum
     q%data%v_(ptId) = evalTriQuality(plMesh%CellPosList(fvmInfo%Point_CellId(1:3, ptId)))
  end do


  sigMin = minval(sigma%data%v_(:) )
  sigAvg = sum(sigma%data%v_(:))/dble(cellNum)
  call MessageNotify("M", "globalSWM::CheckGridQuality", &
       & "minimum of local Uniformity for VoronoiMesh=%f", d=(/ sigMin /) )
  call MessageNotify("M", "globalSWM::CheckGridQuality", &
       & "average of local Uniformity for VoronoiMesh=%f", d=(/ sigAvg /) )

  qMin = minval(q%data%v_(:) )
  qAvg = sum(q%data%v_(:))/dble(ptNum)
  call MessageNotify("M", "globalSWM::CheckGridQuality", &
       & "minimum of local Uniformity for dual DelaunayTriMesh=%f", d=(/ qMin /) )
  call MessageNotify("M", "globalSWM::CheckGridQuality", &
       & "average of local Uniformity for dual DelaunayTriMesh=%f", d=(/ qAvg /) )

  call GeometricField_Final(sigma)
  call GeometricField_Final(q)

stop
contains

function evalTriQuality(pts) result(ret)
  type(Vector3d), intent(in) :: pts(3)
  real(DP) :: ret

  real(DP) :: a, b, c, s

  a = geodesicArcLength(pts(1), pts(2))
  b = geodesicArcLength(pts(2), pts(3))
  c = geodesicArcLength(pts(3), pts(1))
  s = a + b + c
  ret = (s-2d0*a)*(s-2d0*b)*(s-2d0*c)/(a*b*c)

end function evalTriQuality

end subroutine checkGridQuality


end module Exp_Williamson94_Case2
