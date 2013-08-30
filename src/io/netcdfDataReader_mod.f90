module netcdfDataReader_mod
  
  use VectorSpace_mod
  use PolyMesh_mod
  use fvMeshInfo_mod
  use netcdfDataHelper_mod

  use dc_types
  use gtool_history
  
  use GeometricField_mod

  implicit none
  private

  type, public :: netcdfDataReader
    type(fvMeshInfo), pointer :: mesh => null()
    character(STRING) :: filePath
    
    type(Mesh2_ncInfo) :: mesh2Info

  end type netcdfDataReader

  interface netcdfDataReader_get
     module procedure netcdfDataReader_readvScalarData
  end interface netcdfDataReader_get

  public :: netcdfDataReader_Init, netcdfDataReader_Final
  public :: netcdfDataReader_get


contains
subroutine netcdfDataReader_Init(reader, fileName, mesh)
  type(netcdfDataReader), intent(inout) ::reader
  character(*), intent(in) :: fileName
  type(fvMeshInfo), intent(in), target :: mesh

  reader%mesh => mesh
  reader%filePath = fileName
  call Mesh2_ncInfo_Init(reader%mesh2Info, reader%mesh%mesh)
  
  call read_meshData(reader)

end subroutine netcdfDataReader_Init

subroutine netcdfDataReader_readvScalarData(reader, fieldName, field)
  type(netcdfDataReader), intent(in) ::reader
  character(*), intent(in) :: fieldName
  type(volScalarField), intent(inout) :: field

  call GeometricField_Init(field, reader%mesh%mesh, fieldName) 
  call HistoryGet(reader%filePath, fieldName, field%data%v_(:))
  call HistoryGetAttr(reader%filePath, fieldName, 'long_name', field%long_name)
  call HistoryGetAttr(reader%filePath, fieldName, 'units', field%units)

end subroutine netcdfDataReader_readvScalarData

subroutine netcdfDataReader_Final(reader)
  type(netcdfDataReader), intent(inout) :: reader

end subroutine netcdfDataReader_Final

!
!
subroutine read_meshData(reader)
  use SphericalCoord_mod

  type(netcdfDataReader), intent(inout) :: reader


  real(DP), pointer :: vlon(:) => null()
  real(DP), pointer :: vlat(:) => null()
  integer :: cellSize, i

  call HistoryGetPointer(reader%filePath, reader%mesh2Info%face_x%element_name, vlon)
  call HistoryGetPointer(reader%filePath, reader%mesh2Info%face_y%element_name, vlat)
  
  cellSize = size(vlon)
!!$  
!!$  call PolyMeshInit(reader%mesh%mesh, 0, 0, cellSize)
!!$  do i=1, cellSize
!!$!     call PolyMesh_SetCell(
!!$  end do
!!$ 
!!$  call fvMeshInfo_Init(reader%mesh, reader%mesh%mesh)

  deallocate(vlon, vlat)

end subroutine read_meshData

end module netcdfDataReader_mod
