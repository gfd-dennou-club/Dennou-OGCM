program gridGen

  use dc_types
  use VectorSpace_mod
  use PolyMesh_mod
  use HexTriIcMesh_mod
  use GeometricField_mod
  use fvMeshInfo_mod
  use netcdfDataWriter_mod
  use vtkDataWriter_mod


  implicit none

  type(HexTriIcMesh) :: htiMesh
  type(fvMeshInfo) :: fvInfo
  type(netcdfDataWriter) :: ncwriter
  type(vtkDataWriter) :: vtkWriter

  integer, parameter :: glevel = 4
  integer, parameter :: maxItrNum = 2800
  logical, parameter :: outputVtkData = .true.
  character(STRING), parameter :: fileName = 'grid-glevel4'

  ! Setup
  write(*,*) 'hexgonal grid generation.. :glevel:', glevel
  call HexTriIcMesh_Init(htiMesh)
  call HexTriIcMesh_generate(htiMesh, glevel, maxItrNum)

  call fvMeshInfo_Init(fvInfo, htiMesh%mesh)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvInfo)

  !
  write(*,*) "* output netcdf .."
  call netcdfDataWriter_Init(ncwriter, trim(fileName)//".nc", htiMesh%mesh)
  call netcdfDataWriter_Regist(ncwriter, (/ fvInfo%v_CellVol /))
  call netcdfDataWriter_write(ncwriter, fvInfo%v_CellVol)
  call netcdfDataWriter_Final(ncwriter)

  if( outputVtkData ) then
     call vtkDataWriter_Init(vtkwriter, trim(fileName)//".vtk", htiMesh%mesh)
     call vtkDataWriter_Regist(vtkwriter, (/ fvInfo%v_CellVol /))
     call vtkDataWriter_write(vtkwriter)
     call vtkDataWriter_Final(vtkwriter)
  end if

  ! Finalization
  call HexTriIcMesh_Final(htiMesh)
  

end program gridGen
