program advTest
  use dc_types
  use VectorSpace_mod
  use PolyMesh_mod
  use HexTriIcMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use fvCalculus_mod
  use vtkDataWriter_mod
  use netcdfDataWriter_mod
  use netcdfDataReader_mod

  implicit none

  type(HexTriIcMesh) :: htiMesh
  type(fvMeshInfo) :: fvInfo
  integer, parameter :: glevel = 4
  type(vtkDataWriter) :: writer
  type(netcdfDataReader) :: ncReader
  type(PolyMesh) :: plMesh
  type(volScalarField) :: v_CellVol

  call netcdfDataReader_Init(ncReader, "gridData.nc", plMesh)
  call netcdfDataReader_Final(ncReader)
  call HexTriIcMesh_Init(htiMesh, plMesh)


  call fvMeshInfo_Init(fvInfo, plMesh)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)

  !
  !

  call fvCalculus_Init(fvInfo)
  call fvCalculus_Final()


  call vtkDataWriter_Init(writer, "data.nc", plMesh)
  call vtkDataWriter_Regist(writer, (/ fvInfo%v_CellVol /))
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)



  !
  !
  call fvMeshInfo_Final(fvInfo)
  call HexTriIcMesh_Final(htiMesh)
  call PolyMesh_Final(plMesh)

contains


end program advTest
