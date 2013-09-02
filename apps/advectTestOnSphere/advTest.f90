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
  integer, parameter :: dt       = 900
  integer, parameter :: endTime  = 12*24*3600
  integer, parameter :: endTStep = endTime / dt 
  integer, parameter :: outputIntrVal = 43200

  integer :: tStep

  type(volScalarField) :: v_height
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
  call GeometricField_Init(s_normalVel, plmesh, "normalVel", "velocity normal to surface", "ms-1")

  !
  call MessageNotify( 'M', "advectionTest", "Set initial condition..")
  call Set_InitialCondition(v_height, s_normalVel)
  
  call OutputData(0)

  do tStep=1, endTStep
     if( mod(tStep, endTStep/10) == 0 ) then
        call MessageNotify( 'M', "advectionTest", "t=%d [s] (endTime=%d [s]) ..", i=(/ tStep*dt, endTime /) )
     end if


     if( mod(tStep*dt, outputIntrVal) == 0 ) call OutputData(tStep)

     call Solve_AdvetEquation

  end do


  ! Finalize
  !
  call GeometricField_Final(v_height)
  call GeometricField_Final(s_normalVel)

  call fvCalculus_Final()
  call fvMeshInfo_Final(fvInfo)
  call HexTriIcMesh_Final(htiMesh)
  call PolyMesh_Final(plMesh)

contains

subroutine Setup_grid()

  use netcdfDataReader_mod

  type(netcdfDataReader) :: ncReader
  integer :: i

  call netcdfDataReader_Init(ncReader, gridFilePath, plMesh)
  call netcdfDataReader_Final(ncReader)
  call HexTriIcMesh_Init(htiMesh, plMesh, radius)


!do i=1, getPointListSize(plMesh)
!   write(*,*) i, plmesh%pointList(i)
!end do

end subroutine Setup_grid

subroutine Solve_AdvetEquation()

  type(volScalarField) :: v_Tend_height(4)
  type(surfaceScalarField) :: s_Tend_normalVel

  integer :: i
  do i=1, 4
     call GeometricField_Init(v_Tend_height(i), plmesh, "v_Tend_height", "tendency of height", "ms-1")
  end do

  v_Tend_height(1) = div( v_height, s_normalVel )
  v_Tend_height(2) =  div( v_height - (0.5d0*dble(dt))*v_Tend_height(1), s_normalVel )
  v_Tend_height(3) =  div( v_height - (0.5d0*dble(dt))*v_Tend_height(2), s_normalVel )
  v_Tend_height(4) =  div( v_height - dble(dt)*v_Tend_height(3), s_normalVel )
  
  v_height = v_height - dt/6d0 * ( v_Tend_height(1) + 2d0*v_Tend_height(2) + 2d0*v_Tend_height(3) + v_Tend_height(4))
 

 do i=1, 4
     call GeometricField_Final(v_Tend_height(i))
  end do

end subroutine Solve_AdvetEquation

subroutine Set_InitialCondition(v_height, s_normalVel)

  use SphericalCoord_mod

  type(volScalarField), intent(inout) :: v_height
  type(surfaceScalarField), intent(inout) :: s_normalVel

  real(DP) :: u0
  integer :: faceNum, faceId, cellNum, cellId
  type(Vector3d) :: geo_vel, geo_facePos, faceNormal

  u0 = 2d0*acos(-1d0)*radius / (3600d0*24d0*12d0)
write(*,*) "U0=", u0

  cellNum = getCellListSize(plMesh)
  faceNum = getFaceListSize(plMesh)

  do faceId=1, faceNum
     geo_facePos = cartToSphPos( fvInfo%s_faceCenter.At.faceId )
     geo_vel = (/ u0*geo_facePos%v_(1), 0d0, 0d0 /)
     faceNormal = normalizedVec( fvInfo%s_faceAreaVec.At.faceId )
     s_normalVel%data%v_(faceId) = SphToCartVec( geo_vel, fvInfo%s_faceCenter.At.faceId) .dot. faceNormal
  end do

  do cellId=1, cellNum
     v_height%data%v_(cellId) = 1d0
  end do
  
end subroutine Set_InitialCondition

subroutine OutputData( tstep )

  use vtkDataWriter_mod
  use netcdfDataWriter_mod
  use dc_string

  integer, intent(in) :: tstep
  character(STRING) :: dataFileName
  type(vtkDataWriter) :: writer

  dataFileName = CPrintf("data-%08d.vtk", i=(/ tstep /))
  call vtkDataWriter_Init(writer, dataFileName, plMesh)
  call vtkDataWriter_Regist(writer, (/ v_height, fvInfo%v_CellVol /))
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)


end subroutine OutputData

end program advTest
