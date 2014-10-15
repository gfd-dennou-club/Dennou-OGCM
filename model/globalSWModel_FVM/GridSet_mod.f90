module GridSet_mod
  use dc_types
  use dc_message

  use PolyMesh_mod
  use fvMeshInfo_mod
  use HexTriIcMesh_mod

  use SimParameters_mod, only: gridFilePath, Radius

  implicit none

  type(PolyMesh) :: plMesh  
  type(fvMeshInfo) :: fvmInfo
  type(HexTriIcMesh) :: htiMesh


contains
subroutine GridSet_Init()

  ! Load grid data of HexTriIcMesh from netcdf file. 
  call load_gridData()
  call HexTriIcMesh_Init(htiMesh, plMesh, Radius)

  ! Initialize some modules for finite volume method
  call MessageNotify( 'M', "globalSWM", "Initialize some modules for finite volume method..")
  call fvMeshInfo_Init(fvmInfo, plMesh, dualMeshFlag=.true.)
  call HexTriIcMesh_configfvMeshInfo(htiMesh, fvmInfo)

end subroutine GridSet_Init

subroutine GridSet_Final()

  call fvMeshInfo_Final(fvmInfo)
  call HexTriIcMesh_Final(htiMesh)
  call PolyMesh_Final(plMesh)

end subroutine GridSet_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine load_gridData()

  use netcdfDataReader_mod

  type(netcdfDataReader) :: ncReader

  call MessageNotify( 'M', "globalSWM", "Set up grid. Load grid data from '%a' ..", ca=(/ gridFilePath /) ) 

  call netcdfDataReader_Init(ncReader, gridFilePath, plMesh)
  call netcdfDataReader_Final(ncReader)


end subroutine Load_gridData

end module GridSet_mod
