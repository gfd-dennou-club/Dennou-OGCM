module netcdfDataWriter_mod

  use VectorSpace_mod
  use PolyMesh_mod
  use fvMeshInfo_mod
  use netcdfDataHelper_mod

  use dc_types
  use gtool_history
  
  use GeometricField_mod

  implicit none
  private

  type, public :: netcdfDataWriter
    type(fvMeshInfo), pointer :: mesh => null()
    character(STRING) :: filePath
    
    type(VolScalarField), pointer :: ref_VScalarFList(:) => null()
    type(VolVectorField), pointer :: ref_VVectorFList(:) => null()

    type(Mesh2_ncInfo) :: mesh2Info

  end type netcdfDataWriter

  interface netcdfDataWriter_write
     module procedure write_volScalarField
  end interface netcdfDataWriter_write

  public :: netcdfDataWriter_Init, netcdfDataWriter_Final
  public :: netcdfDataWriter_Regist
  public :: netcdfDataWriter_Write

contains
subroutine netcdfDataWriter_Init(writer, fileName, mesh)
  type(netcdfDataWriter), intent(inout) :: writer
  character(*), intent(in) :: fileName
  type(fvMeshInfo), intent(in), target :: mesh

  writer%mesh => mesh
  writer%filePath = fileName
  call Mesh2_ncInfo_Init(writer%mesh2Info, writer%mesh%mesh)
  call netcdfDataWriter_writeMetaData(writer)
  call netcdfDataWriter_writeGridData(writer)

end subroutine netcdfDataWriter_Init


subroutine netcdfDataWriter_Regist(writer, &
  & volScalarFields, volVectorFields )

  type(netcdfDataWriter), intent(inout), target :: writer
  type(volScalarField), target, optional :: volScalarFields(:)
  type(volVectorField), target, optional :: volVectorFields(:)

  integer :: i
  type(Mesh2_ncInfo), pointer :: meshInfo
  
  meshInfo => writer%mesh2Info
  
  if(present(volScalarFields)) then
    writer%ref_VScalarFList => volScalarFields
    do i=1, size(volScalarFields)
       call HistoryAddVariable( &
            & varname=volScalarFields(i)%name, &
            & dims=(/ meshInfo%face%element_name, 't' /), &
            & longname=volScalarFields(i)%long_name, & 
            & units=volScalarFields(i)%units,  xtype='double')
    end do
  end if

  if(present(volVectorFields)) then
    writer%ref_VVectorFList => volVectorFields
  end if

end subroutine netcdfDataWriter_Regist

subroutine write_volScalarField(writer, v_scalar)
  type(netcdfDataWriter), intent(inout) :: writer
  type(volScalarField), intent(in) :: v_scalar

  call HistoryPut(v_scalar%name, v_scalar%data%v_(:))

end subroutine write_volScalarField

subroutine netcdfDataWriter_Final(writer)
  type(netcdfDataWriter), intent(inout) :: writer

  call HistoryClose()

  writer%mesh => null()

end subroutine netcdfDataWriter_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine netcdfDataWriter_writeMetaData(writer)
  type(netcdfDataWriter), intent(inout), target :: writer

  integer :: pointNum, cellNum
  type(Mesh2_ncInfo), pointer :: meshInfo
  
  meshInfo => writer%mesh2Info
  cellNum = getCellListSize(writer%mesh%mesh)
  pointNum = getPointListSize(writer%mesh%mesh)

  call HistoryCreate( &
       & file=writer%filePath, title=" ", &
       & source='Sample program of gtool_history/gtool5',   &
       & institution='GFD_Dennou Club davis project',       &
       & dims=(/ 'nMesh2_point', 'nMesh2_cell ', 't           '/), &
       & dimsizes=(/ pointNum, cellNum, 0/), &
       & longnames=(/'point', 'cell ', 'time ' /),       &
       & units=(/'1', '1', 's'/),   &
       & origin=real(0), interval=real(1) )

end subroutine netcdfDataWriter_writeMetaData

subroutine netcdfDataWriter_writeGridData(writer)
  use SphericalCoord_mod

  type(netcdfDataWriter), intent(inout), target :: writer

  type(Mesh2_ncInfo), pointer :: meshInfo
  type(PolyMesh), pointer :: mesh
  type(volScalarField) :: v_lon, v_lat
  integer :: cellNum, cellId
  type(Vector2d) :: geoPos
  integer :: i

  mesh => writer%mesh%mesh
  meshInfo => writer%mesh2Info
  cellNum = getCellListSize(mesh)

  call GeometricField_Init(v_lon, mesh, "v_lon")
  call GeometricField_Init(v_lat, mesh, "v_lat")
  do i=1, cellNum
     geoPos = RadToDegUnit( cartToSphPos( writer%mesh%v_cellVec.At.i ) )
     v_lon%data%v_(i) = geoPos%v_(1)
     v_lat%data%v_(i) = geoPos%v_(2)
  end do

  call HistoryAddVariable( &
       & varname=meshInfo%face_x%element_name, &
       & dims=(/ meshInfo%face_x%dimension_element%element_name /), &
       & longname=meshInfo%face_x%long_name, & 
       & units=meshInfo%face_x%units,  xtype='double')

  call HistoryPut(meshInfo%face_x%element_name, v_lon%data%v_(:))

  call HistoryAddVariable( &
       & varname=meshInfo%face_y%element_name, &
       & dims=(/ meshInfo%face_y%dimension_element%element_name /), &
       & longname=meshInfo%face_y%long_name, & 
       & units=meshInfo%face_y%units,  xtype='double')

  call HistoryPut(meshInfo%face_y%element_name, v_lat%data%v_(:))

  call GeometricField_Final(v_lon)
  call GeometricField_Final(v_lat)

end subroutine netcdfDataWriter_writeGridData


end module netcdfDataWriter_mod
