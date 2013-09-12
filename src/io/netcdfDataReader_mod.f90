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
    type(PolyMesh), pointer :: mesh => null()
    character(STRING) :: filePath
    
    type(Mesh2_ncInfo) :: mesh2Info

  end type netcdfDataReader

  interface netcdfDataReader_get
     module procedure netcdfDataReader_readvScalarData
     module procedure netcdfDataReader_readpScalarData
  end interface netcdfDataReader_get

  public :: netcdfDataReader_Init, netcdfDataReader_Final
  public :: netcdfDataReader_get


contains
subroutine netcdfDataReader_Init(reader, fileName, mesh)
  type(netcdfDataReader), intent(inout) ::reader
  character(*), intent(in) :: fileName
  type(PolyMesh), intent(in), target :: mesh

  reader%mesh => mesh
  reader%filePath = fileName
  call Mesh2_ncInfo_Init(reader%mesh2Info, reader%mesh)
  
  call read_meshData(reader)

end subroutine netcdfDataReader_Init

subroutine netcdfDataReader_readvScalarData(reader, fieldName, field, range)
  type(netcdfDataReader), intent(in) ::reader
  character(*), intent(in) :: fieldName
  type(volScalarField), intent(inout) :: field
  character(*), intent(in), optional :: range

  call GeometricField_Init(field, reader%mesh, fieldName) 
  call HistoryGet(reader%filePath, fieldName, field%data%v_(:), range=range)
  call HistoryGetAttr(reader%filePath, fieldName, 'long_name', field%long_name)
  call HistoryGetAttr(reader%filePath, fieldName, 'units', field%units)

end subroutine netcdfDataReader_readvScalarData

subroutine netcdfDataReader_readpScalarData(reader, fieldName, field, range)
  type(netcdfDataReader), intent(in) ::reader
  character(*), intent(in) :: fieldName
  type(pointScalarField), intent(inout) :: field
  character(*), intent(in), optional :: range

  call GeometricField_Init(field, reader%mesh, fieldName) 
  call HistoryGet(reader%filePath, fieldName, field%data%v_(:), range=range)
  call HistoryGetAttr(reader%filePath, fieldName, 'long_name', field%long_name)
  call HistoryGetAttr(reader%filePath, fieldName, 'units', field%units)

end subroutine netcdfDataReader_readpScalarData

subroutine netcdfDataReader_Final(reader)
  type(netcdfDataReader), intent(inout) :: reader

end subroutine netcdfDataReader_Final

!
!
subroutine read_meshData(reader)
  use SphericalCoord_mod

  type(netcdfDataReader), intent(inout), target :: reader


  real(DP), pointer :: vlon(:) => null()
  real(DP), pointer :: vlat(:) => null()
  real(DP), pointer :: plon(:) => null()
  real(DP), pointer :: plat(:) => null()

  integer :: cellNum, ptNum, faceNum, i
  type(Cell) :: tmpCell
  type(Vector3d) :: cellPos, ptPos
  type(Vector3d) :: latlon
  type(Mesh2_ncInfo), pointer :: meshInfo
  integer, pointer :: cell_points(:,:) => null()
  integer, pointer :: cell_faces(:,:) => null()
  integer, pointer :: face_points(:,:) => null()
  integer, pointer :: face_links(:,:) => null()


  meshInfo => reader%mesh2Info

  call HistoryGetPointer(reader%filePath, meshInfo%node_x%element_name, plon)
  call HistoryGetPointer(reader%filePath, meshInfo%node_y%element_name, plat)
  call HistoryGetPointer(reader%filePath, meshInfo%face_x%element_name, vlon)
  call HistoryGetPointer(reader%filePath, meshInfo%face_y%element_name, vlat)
  !
  call HistoryGetPointer(reader%filePath, meshInfo%face_nodes%element_name, cell_points)
  call HistoryGetPointer(reader%filePath, meshInfo%face_edges%element_name, cell_faces)
  call HistoryGetPointer(reader%filePath, meshInfo%edge_nodes%element_name, face_points)
  call HistoryGetPointer(reader%filePath, meshInfo%face_links%element_name, face_links)

  
  cellNum = size(vlon)
  ptNum = size(plon)
  faceNum = size(face_points, 2)
  
  call PolyMesh_Init(reader%mesh, ptNum, faceNum, cellNum)
  do i=1, cellNum
     latlon = (/ vlon(i), vlat(i), 1d0 /)
     cellPos = SphToCartPos(DegToRadUnit(latlon))

     call PolyMesh_SetCell(reader%mesh, i, tmpCell, cellPos)
  end do

  do i=1, ptNum
     latlon = (/ plon(i), plat(i), 1d0 /)
     ptPos = SphToCartPos(DegToRadUnit(latlon))
     call PolyMesh_SetPoint(reader%mesh, i, ptPos)
  end do 


  call PolyMesh_setConnectivity(reader%mesh, cell_points, cell_faces, face_points, face_links)

  deallocate(vlon, vlat)
  deallocate(plon, plat)
  deallocate(cell_faces, cell_points, face_points, face_links)

end subroutine read_meshData

end module netcdfDataReader_mod
