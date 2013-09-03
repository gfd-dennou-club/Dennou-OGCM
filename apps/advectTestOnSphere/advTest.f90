program advTest

  use dc_types
  use dc_message
  use VectorSpace_mod
  use PolyMesh_mod
  use HexTriIcMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use fvCalculus_mod

  implicit none

  type(HexTriIcMesh) :: htiMesh
  type(fvMeshInfo) :: fvInfo
  type(PolyMesh) :: plMesh
  type(volScalarField) :: v_CellVol
  character(STRING) :: gridFilePath = "grid/grid-glevel5.nc"

  real(DP), parameter :: radius  = 6.37122d06
  real(DP), parameter :: Omega = 7.292d-05
  real(DP), parameter :: Grav = 9.80616
  real(DP), parameter :: meanDepth = 2000d0

  integer, parameter :: dt       = 1200
  integer, parameter :: endTime  = 12*24*3600
  integer, parameter :: endTStep = endTime / dt 
  integer, parameter :: outputIntrVal = 43200
  real(DP), parameter :: PI = acos(-1d0)

  integer :: tStep

  type(volScalarField) :: v_height, v_div
  type(surfaceScalarField) :: s_normalVel
  
  !
  !

  call MessageNotify( 'M', "advectionTest", "Set up grid. Load grid data from '%a' ..", ca=(/ gridFilePath /) ) 
  call Setup_grid()

  ! Initialize some modules for finite volume method
  call MessageNotify( 'M', "advectionTest", "Initialize some modules for finite volume method..")
  call fvMeshInfo_Init(fvInfo, plMesh)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)
  call fvCalculus_Init(fvInfo)



  !
  call GeometricField_Init(v_height, plmesh, "height", "height", "m")
  call GeometricField_Init(v_div, plmesh, "div", "div", "1")
  call GeometricField_Init(s_normalVel, plmesh, "normalVel", "velocity normal to surface", "ms-1")

  !
  call MessageNotify( 'M', "advectionTest", "Set initial condition..")
  call Set_InitialCondition(v_height, s_normalVel)
  
  call OutputData(0)

  do tStep=1, endTStep
     if( mod(tStep, endTStep/20) == 0 ) then
        call MessageNotify( 'M', "advectionTest", "t=%d [s] (endTime=%d [s]) ..", i=(/ tStep*dt, endTime /) )
     end if


     if( mod(tStep*dt, outputIntrVal) == 0 ) call OutputData(tStep)

     call Solve_AdvetEquation

  end do


  ! Finalize
  !
  call GeometricField_Final(v_height)
  call GeometricField_Final(v_div)
  call GeometricField_Final(s_normalVel)

  call fvCalculus_Final()
  call fvMeshInfo_Final(fvInfo)
  call HexTriIcMesh_Final(htiMesh)
  call PolyMesh_Final(plMesh)

contains

subroutine Setup_grid()

  use netcdfDataReader_mod

  type(netcdfDataReader) :: ncReader

  call netcdfDataReader_Init(ncReader, gridFilePath, plMesh)
  call netcdfDataReader_Final(ncReader)
  call HexTriIcMesh_Init(htiMesh, plMesh, radius)


end subroutine Setup_grid

subroutine Solve_AdvetEquation()

  type(volScalarField) :: v_mdhdt, v_height0
  type(surfaceScalarField) :: s_mdvdt, s_normalVel0

  integer :: i
  real(DP) :: dt_

  dt_ = dble(dt)
  call GeometricField_Init(v_mdhdt, plmesh, "v_Tend_height", "tendency of height", "ms-1")
  call GeometricField_Init(v_height0, plMesh, "v_height0", "height of current tstep", "m")
  call GeometricField_Init(s_mdvdt, plmesh, "s_Tend_normalVel", "tendency of velocity", "ms-2")
  call GeometricField_Init(s_normalVel0, plMesh, "s_normalVel0", "normal velocity of current tstep", "ms-1")


  call DeepCopy(v_height0,  v_height)
  call DeepCopy(s_normalVel0, s_normalVel)

  call calcTendency(v_mdhdt, s_mdvdt, v_height, s_normalVel )
  v_height = v_height0 - (dt_/6d0) * v_mdhdt

  call calcTendency(v_mdhdt, s_mdvdt, v_height0 + (-0.5d0*dt_)*v_mdhdt, s_normalVel )
  v_height = v_height - (2d0*dt_/6d0) * v_mdhdt

  call calcTendency(v_mdhdt, s_mdvdt, v_height0 + (-0.5d0*dt_)*v_mdhdt, s_normalVel )
  v_height = v_height - (2d0*dt_/6d0) * v_mdhdt

  call calcTendency(v_mdhdt, s_mdvdt, v_height0 + (-dt_)*v_mdhdt, s_normalVel )
  v_height = v_height - (dt_/6d0) * v_mdhdt


  call GeometricField_Final(v_mdhdt)
  call GeometricField_Final(s_mdvdt)
  call GeometricField_Final(v_height0)
  call GeometricField_Final(s_normalVel0)


end subroutine Solve_AdvetEquation

subroutine calcTendency(dmhdt, dmvdt, h, v)
  type(volScalarField), intent(inout) :: dmhdt
  type(surfaceScalarField), intent(inout) :: dmvdt
  type(volScalarField), intent(in) :: h
  type(surfaceScalarField), intent(in) :: v

  dmhdt = div(h, v)

end subroutine calcTendency



subroutine Set_InitialCondition(v_height, s_normalVel)

  use SphericalCoord_mod

  type(volScalarField), intent(inout) :: v_height
  type(surfaceScalarField), intent(inout) :: s_normalVel

  real(DP) :: u0
  real(DP), parameter :: alpha = 0.5*acos(-1d0)
  integer :: faceNum, faceId, cellNum, cellId, lfaceId
  type(Vector3d) :: geo_vel, geoPos, faceNormal, vxs(2)
  real(DP) :: distor, lmean
  integer :: maxId(1)
  real(DP), allocatable :: r
  type(Vector3d) :: centerPos

  real(DP), parameter :: cosBellRadius = radius/3d0
  real(DP), parameter :: h0 = 1000d0

  u0 = 2d0*PI*radius / (3600d0*24d0*12d0)
write(*,*) "U0=", u0

  cellNum = getCellListSize(plMesh)
  faceNum = getFaceListSize(plMesh)

  do faceId=1, faceNum
     geoPos = cartToSphPos( fvInfo%s_faceCenter.At.faceId )
     geo_vel = (/ u0*(cos(geoPos%v_(2))), 0d0, 0d0 /)
     faceNormal = normalizedVec( fvInfo%s_faceAreaVec.At.faceId )
     s_normalVel%data%v_(faceId) = SphToCartVec( geo_vel, fvInfo%s_faceCenter.At.faceId) .dot. faceNormal

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


end subroutine Set_InitialCondition

subroutine OutputData( tstep )

  use vtkDataWriter_mod
  use netcdfDataWriter_mod
  use dc_string
  use SphericalCoord_mod

  integer, intent(in) :: tstep
  character(STRING) :: dataFileName
  type(vtkDataWriter) :: writer

integer :: maxId(1)
maxId = maxloc(v_height%data%v_)
write(*,*) "max height=", v_height.At.maxId(1), ":", RadToDegUnit(CartToSphPos(plmesh%cellPosList(maxId(1))))
  dataFileName = CPrintf("data-%08d.vtk", i=(/ tstep /))
  call vtkDataWriter_Init(writer, dataFileName, plMesh)
  call vtkDataWriter_Regist(writer, (/ v_height, fvInfo%v_CellVol, v_div /))
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)


end subroutine OutputData

end program advTest
