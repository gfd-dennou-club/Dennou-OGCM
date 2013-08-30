program advTest

  use dc_types
  use VectorSpace_mod
  use HexTriIcMesh_mod
  use fvMeshInfo_mod
  use fvCalculus_mod
  use vtkDataWriter_mod
  use netcdfDataWriter_mod

  implicit none

  type(HexTriIcMesh) :: htiMesh
  type(fvMeshInfo) :: fvInfo
  integer, parameter :: glevel = 4
  type(vtkDataWriter) :: writer
  type(netcdfDataWriter) :: ncwriter

  write(*,*) 'glevel:', glevel
  call HexTriIcMesh_Init(htiMesh, glevel)
  call HexTriIcMesh_generate(htiMesh)

  call fvMeshInfo_Init(fvInfo, htiMesh%mesh)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)

  !
  !

  call fvCalculus_Init(fvInfo)
  call fvCalculus_Final()

  !
  !
  call vtkDataWriter_Init(writer, "data.vtk", htiMesh%mesh)
  call vtkDataWriter_Regist(writer, (/ fvInfo%v_CellVol /))
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)

  !
  call netcdfDataWriter_Init(ncwriter, "data.nc", fvInfo)
  call netcdfDataWriter_Regist(ncwriter, (/ fvInfo%v_CellVol /))
  call netcdfDataWriter_write(ncwriter, fvInfo%v_CellVol)
  call netcdfDataWriter_Final(ncwriter)


  !
  !
  call fvMeshInfo_Final(fvInfo)
  call HexTriIcMesh_Final(htiMesh)

contains
end program advTest
