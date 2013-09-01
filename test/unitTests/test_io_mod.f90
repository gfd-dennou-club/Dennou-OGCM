program test_io_mod

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

  type(HexTriIcMesh) :: htiMesh, htiMesh_r

  type(fvMeshInfo) :: fvInfo, fvInfo_r
  integer, parameter :: glevel = 4
  integer, parameter :: maxItrNum = 10
  type(vtkDataWriter) :: writer
  type(netcdfDataWriter) :: ncwriter
  type(netcdfDataReader) :: ncReader
  type(volScalarField) :: v_CellVol
  type(PolyMesh) :: pMesh


  ! Setup

  write(*,*) 'hexgonal grid generation.. :glevel:', glevel
  call HexTriIcMesh_Init(htiMesh)
  call HexTriIcMesh_generate(htiMesh, glevel, maxItrNum)

  call fvMeshInfo_Init(fvInfo, htiMesh%mesh)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)


  !
  !
  
  !
  write(*,*) "vtkDataWriter test.."
  call vtkDataWriter_Init(writer, "data.vtk", htiMesh%mesh)
  call vtkDataWriter_Regist(writer, (/ fvInfo%v_CellVol /))
  call vtkDataWriter_write(writer)
  call vtkDataWriter_Final(writer)

  !
  write(*,*) "* netcdfDataWriter test.."
  call netcdfDataWriter_Init(ncwriter, "data.nc", htiMesh%mesh)
  call netcdfDataWriter_Regist(ncwriter, (/ fvInfo%v_CellVol /))
  call netcdfDataWriter_write(ncwriter, fvInfo%v_CellVol)
  call netcdfDataWriter_Final(ncwriter)

  write(*,*) "* netcdfDataReader test.."
  call netcdfDataReader_Init(ncReader, "data.nc", pMesh)
  call netcdfDataReader_get(ncReader, "v_CellVol", v_CellVol)
  call netcdfDataReader_Final(ncReader)

  write(*,*) "Check the data read from netcdf file.."
  call dataCheck()

  call HexTriIcMesh_Init(htiMesh_r, pMesh)
!  call HexTriIcMesh_configfvMeshInfo(htiMesh_r, fvInfo_r)
  call HexTriIcMesh_Final(htiMesh_r)

  ! Finalization
  call PolyMesh_Final(pmesh)
  call GeometricField_Final(v_CellVol)
  call HexTriIcMesh_Final(htiMesh)

contains
subroutine setup()

end subroutine setup

subroutine dataCheck()
  integer :: i

  do i=1, getCellListSize(htiMesh%mesh)
     if( (v_CellVol.At.i) /= (fvInfo%v_CellVol.At.i) ) stop
  end do
end subroutine dataCheck

end program test_io_mod
