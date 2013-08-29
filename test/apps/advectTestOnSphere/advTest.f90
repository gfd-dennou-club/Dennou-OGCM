program advTest

  use dc_types
  use VectorSpace_mod
  use HexTriIcMesh_mod
  use fvMeshInfo_mod
  use fvCalculus_mod
  use vtkDataWriter_mod

  implicit none

  type(HexTriIcMesh) :: htiMesh
  type(fvMeshInfo) :: fvInfo
  integer, parameter :: glevel = 4
  type(vtkDataWriter) :: writer

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
!  call vtkDataWriter_Regist(writer, volScalarFields=(/ fvInfo%v_CellVol /)) 
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)


  !
  !
  call fvMeshInfo_Final(fvInfo)
  call HexTriIcMesh_Final(htiMesh)

contains
end program advTest
