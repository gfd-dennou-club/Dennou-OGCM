module netcdfDataWriter_mod

  use VectorSpace_mod
  use PolyMesh_mod
  use fvMeshInfo_mod
  use netcdfDataHelper_mod

  use dc_types
  use gtool_history
  
  use GeometricField_mod
  use netcdf

  implicit none
  private

  type, public :: netcdfDataWriter
    type(PolyMesh), pointer :: mesh => null()
    character(STRING) :: filePath
    type(Mesh2_ncInfo) :: mesh2Info
    integer :: ncId
    integer :: rec_dimId
    integer :: recCounter = 1
  end type netcdfDataWriter

  interface netcdfDataWriter_write
     module procedure write_volScalarField
     module procedure write_ptScalarField
  end interface netcdfDataWriter_write

  public :: netcdfDataWriter_Init, netcdfDataWriter_Final
  public :: netcdfDataWriter_Regist
  public :: netcdfDataWriter_Write
  public :: netcdfDataWriter_WriteGlobalAttr
  public :: netcdfDataWriter_AdvanceTimeStep


  character(TOKEN), parameter :: RECODE_NAME = "time"

contains
subroutine netcdfDataWriter_Init(writer, fileName, mesh)
  type(netcdfDataWriter), intent(inout) :: writer
  character(*), intent(in) :: fileName
  type(PolyMesh), intent(in), target :: mesh

  writer%mesh => mesh
  writer%filePath = fileName
  call Mesh2_ncInfo_Init(writer%mesh2Info, writer%mesh)


  !
  call check_nf90_status( &
       & nf90_create(writer%filePath, NF90_CLOBBER, writer%ncID) )


end subroutine netcdfDataWriter_Init


subroutine netcdfDataWriter_Regist(writer, &
  & volScalarFields, pointScalarFields )

  type(netcdfDataWriter), intent(inout), target :: writer
  type(volScalarField), target, optional :: volScalarFields(:)
  type(pointScalarField), target, optional :: pointScalarFields(:)

  integer :: i, vScalarFListSize, varId
  type(Mesh2_ncInfo), pointer :: meshInfo
  

  call netcdfDataWriter_writeGridMetaData(writer)

  meshInfo => writer%mesh2Info
  
  if(present(volScalarFields)) then

     do i=1, size(volScalarFields)
        varId = netcdfDataWriter_defFieldData( writer, &
             & (/ writer%mesh2Info%layers%dimId, writer%mesh2Info%face%dimId, writer%rec_dimId /), &
             & volScalarFields(i)%name,  volScalarFields(i)%long_name, volScalarFields(i)%units &
             & )
    end do
  end if

  if(present(pointScalarFields)) then
     do i=1, size(pointScalarFields)
        varId = netcdfDataWriter_defFieldData( writer, & 
             & (/ writer%mesh2Info%layers%dimId, writer%mesh2Info%node%dimId, writer%rec_dimId /), &
             & pointScalarFields(i)%name,  pointScalarFields(i)%long_name, pointScalarFields(i)%units &
             & )
    end do
  end if


  call check_nf90_status( nf90_enddef( writer%ncID ) )
  call netcdfDataWriter_writeGridData(writer)
  !call netcdfDataWriter_AdvanceTimeStep(writer, writer%initTime)

end subroutine netcdfDataWriter_Regist

subroutine write_volScalarField(writer, v_scalar)
  type(netcdfDataWriter), intent(inout) :: writer
  type(volScalarField), intent(in) :: v_scalar

  integer :: varId


  call check_nf90_status( &
      & nf90_inq_varid(writer%ncID, v_scalar%name, varID) &
      & )

  call check_nf90_status( nf90_put_var( &
       & writer%ncID, varid, v_scalar%data%v_(1:v_scalar%vLayerNum,:), &
       &  start=(/ 1, 1, writer%recCounter /), &
       & count=(/ v_scalar%vLayerNum, size(v_scalar%data%v_,2), 1 /) ) &
       & )

end subroutine write_volScalarField

subroutine write_ptScalarField(writer, p_scalar)
  type(netcdfDataWriter), intent(inout) :: writer
  type(pointScalarField), intent(in) :: p_scalar

  integer :: varId


  call check_nf90_status( &
      & nf90_inq_varid(writer%ncID, p_scalar%name, varID) &
      & )

  call check_nf90_status( nf90_put_var( &
       & writer%ncID, varid, p_scalar%data%v_(1:p_scalar%vLayerNum,:), &
       &  start=(/ 1, 1, writer%recCounter /), &
       & count=(/ p_scalar%vLayerNum, size(p_scalar%data%v_,2), 1 /) ) &
       & )

end subroutine write_ptScalarField

subroutine netcdfDataWriter_Final(writer)
  type(netcdfDataWriter), intent(inout) :: writer

  call check_nf90_status( nf90_close(writer%ncID) )

  writer%mesh => null()

end subroutine netcdfDataWriter_Final

subroutine netcdfDataWriter_AdvanceTimeStep( writer, newStepTime )
  type(netcdfDataWriter), intent(inout) :: writer
  real(DP), intent(in) :: newStepTime

  integer :: varId

  call check_nf90_status( &
      & nf90_inq_varid(writer%ncID, RECODE_NAME, varID) &
      & )

  call check_nf90_status( nf90_put_var( &
       & writer%ncID, varid, (/ newStepTime /), &
       & start=(/ writer%recCounter /), count=(/ 1 /) ) &
 & )

  writer%recCounter = writer%recCounter + 1
  
end subroutine netcdfDataWriter_AdvanceTimeStep

subroutine netcdfDataWriter_WriteGlobalAttr(writer, attrName, attrVal)
  type(netcdfDataWriter), intent(inout) :: writer
  character(*), intent(in) :: attrName
  real(DP), intent(in) :: attrVal

  call check_nf90_status( &
    & nf90_put_att(writer%ncID, NF90_GLOBAL, attrName, attrVal), &
    message='write a global attribute.'  &
    & )
end subroutine netcdfDataWriter_WriteGlobalAttr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine netcdfDataWriter_writeGridMetaData(writer)
  type(netcdfDataWriter), intent(inout), target :: writer

  integer :: pointNum, cellNum
  type(Mesh2_ncInfo), pointer :: meshInfo
  integer :: timeVarId
  
  meshInfo => writer%mesh2Info
  cellNum = getCellListSize(writer%mesh)
  pointNum = getPointListSize(writer%mesh)


  ! NetCDF ファイルに次元要素 node を定義する.
  call ncdef_dimension(writer, meshInfo%node )

  ! NetCDF ファイルに座標要素 : node_x, node_y を定義する.
  call ncdef_mesh_coordinate(writer, meshInfo%node_x )
  call ncdef_mesh_coordinate(writer, meshInfo%node_y )


  ! NetCDF ファイルに次元要素 face を定義する.
  call ncdef_dimension(writer, writer%Mesh2Info%face )

  ! NetCDF ファイルに座標要素 : face_x, face_y を定義する.
  call ncdef_mesh_coordinate(writer, meshInfo%face_x )
  call ncdef_mesh_coordinate(writer, meshInfo%face_y )

  ! NetCDF ファイルに次元要素 edge を定義する.
  call ncdef_dimension(writer, writer%Mesh2Info%edge )

  ! NetCDF ファイルに座標要素 : edge_x, edge_y を定義する.
  call ncdef_mesh_coordinate(writer, meshInfo%edge_x )
  call ncdef_mesh_coordinate(writer, meshInfo%edge_y )

  ! NetCDF ファイルに次元要素 layers を定義する.
  call ncdef_dimension(writer, writer%Mesh2Info%layers )

  ! NetCDF ファイルに時間座標を定義する.
  call check_nf90_status( &
    & nf90_def_dim(writer%ncID, RECODE_NAME, NF90_UNLIMITED, writer%rec_dimid), &
    message='define dimension.'  &
    & )


  timeVarId = netcdfDataWriter_defFieldData( &
       & writer, (/ writer%rec_dimId /), &
       & RECODE_NAME,  "time", "second" &
       & )
  

  ! Set some information with topology
  ! NetCDF ファイルに次元要素 face を定義する.
  call ncdef_dimension(writer, meshInfo%Max_face_nodes )
  call ncdef_dimension(writer, meshInfo%Two )
  
  call ncdef_mesh_connectivity(writer, meshInfo%face_nodes)
  call ncdef_mesh_connectivity(writer, meshInfo%edge_nodes)
  call ncdef_mesh_connectivity(writer, meshInfo%face_edges)
  call ncdef_mesh_connectivity(writer, meshInfo%face_links)

end subroutine netcdfDataWriter_writeGridMetaData

subroutine netcdfDataWriter_writeGridData(writer)
  use SphericalCoord_mod

  type(netcdfDataWriter), intent(inout), target :: writer

  type(Mesh2_ncInfo), pointer :: meshInfo
  type(PolyMesh), pointer :: mesh
  type(volScalarField) :: v_lon, v_lat
  type(pointScalarField) :: p_lon, p_lat
  type(surfaceScalarField) :: s_lon, s_lat
  integer, allocatable :: cell_points(:,:), cell_faces(:,:), face_points(:,:), face_links(:,:)
  integer :: pointNum, cellNum, cellId
  integer :: ptNum, ptId

  type(Vector3d) :: geoPos
  integer :: i

  mesh => writer%mesh
  meshInfo => writer%mesh2Info
  cellNum = getCellListSize(mesh)
  pointNum = getPointListSize(mesh)

  call GeometricField_Init(v_lon, mesh, "v_lon", vlayerNum=1)
  call GeometricField_Init(v_lat, mesh, "v_lat", vlayerNum=1)
  call GeometricField_Init(p_lon, mesh, "p_lon", vlayerNum=1)
  call GeometricField_Init(p_lat, mesh, "p_lat", vlayerNum=1)
  call GeometricField_Init(s_lon, mesh, "s_lon", vlayerNum=1)
  call GeometricField_Init(s_lat, mesh, "s_lat", vlayerNum=1)


  do i=1, pointNum
     geoPos = RadToDegUnit( cartToSphPos( mesh%pointPosList(i) ) )
     p_lon%data%v_(1,i) = geoPos%v_(1)
     p_lat%data%v_(1,i) = geoPos%v_(2)
  end do

  call check_nf90_status( &
       & nf90_put_var( writer%ncID, meshInfo%node_x%varID, p_lon%data%v_(1,:) ) )

  call check_nf90_status( &
       & nf90_put_var( writer%ncID, meshInfo%node_y%varID, p_lat%data%v_(1,:) ) )

  do i=1, cellNum
     geoPos = RadToDegUnit( cartToSphPos( mesh%cellPosList(i) ) )
     v_lon%data%v_(1,i) = geoPos%v_(1)
     v_lat%data%v_(1,i) = geoPos%v_(2)
  end do

  call check_nf90_status( &
       & nf90_put_var( writer%ncID, meshInfo%face_x%varID, v_lon%data%v_(1,:) ) )

  call check_nf90_status( &
       & nf90_put_var( writer%ncID, meshInfo%face_y%varID, v_lat%data%v_(1,:) ) )

  call GeometricField_Final(v_lon)
  call GeometricField_Final(v_lat)
  call GeometricField_Final(s_lon)
  call GeometricField_Final(s_lat)
  call GeometricField_Final(p_lon)
  call GeometricField_Final(p_lat)

  !
  call PolyMesh_getConnectivity(mesh, cell_points, cell_faces, face_points, face_links)

  call check_nf90_status( &
       nf90_put_var( writer%ncID, meshInfo%face_nodes%varID, cell_points(:,:)) , message="define cell_points" )

  call check_nf90_status( &
       nf90_put_var( writer%ncID, meshInfo%edge_nodes%varID, face_points(1:meshInfo%edge_nodes%connected_dim_element%num,:)), &
       & message="define face_points" )

  call check_nf90_status( &
       nf90_put_var( writer%ncID, meshInfo%face_edges%varID, cell_faces(:,:)) , message="define cell_faces" )

  call check_nf90_status( &
       nf90_put_var( writer%ncID, meshInfo%face_links%varID, face_links(:,:)) , message="define face_links" )

end subroutine netcdfDataWriter_writeGridData


function netcdfDataWriter_defFieldData(writer, dimIds, varname, long_name, units ) result(varId)

  type(netcdfDataWriter), intent(inout) :: writer
  integer, intent(in) :: dimIds(:)
  character(*), intent(in) :: varname, long_name, units
  integer :: varid
  

  ! NetCDF ファイルに, 後に物理場のデータを書き込むことになる物理場の情報を定義する.
  call check_nf90_status( &
       & nf90_def_var(writer%ncID, varname, NF90_DOUBLE, dimIds, varid), &
       message='Define the field `' // trim(varname) // '`.' &
       & )

  ! 物理場の name を設定.
  call check_nf90_status( &
    & nf90_put_att( writer%ncID, varid, 'name', varname ) &
    & )

  ! 物理場の long_name を設定.
  call check_nf90_status( &
    & nf90_put_att( writer%ncID, varid, 'long_name', long_name ) &
    & )

  ! 物理場の単位の設定.
  call check_nf90_status( &
    & nf90_put_att( writer%ncID, varid, 'units', units ) &
    & ) 

end function netcdfDataWriter_defFieldData

!
!> \~japanese
!! \brief NetCDF ファイルのヘッダー部分に, 次元情報を定義する.
!!
!! @param[in]     self        構造型 IcGrid_ncWriter の変数.
!! @param[in,out] dim_element 次元情報を保持する構造型 Mesh_dim_element の変数.
!!
!! \~english
!! \brief Defines the informaion about the dimension in the header part of the NetCDF file.
!!
!! @param[in]     self        The variable of derived type IcGrid_ncWriter.
!! @param[in,out] dim_element The variable of derived type Mesh_dim_element containing the information about the dimension.
!!
subroutine ncdef_dimension(self, dim_element)

  ! 宣言文 ; Declaration statements
  !
  type(netcdfDataWriter), intent(in) :: self
  type(Mesh_dim_element), intent(inout) :: dim_element

  ! 実行文 ; Executable statments
  !

  call check_nf90_status( &
    & nf90_def_dim(self%ncID, dim_element%element_name, dim_element%num, dim_element%dimID), &
    & message='define dimension. ' // trim(dim_element%element_name) &
    & )

end subroutine ncdef_dimension

!
!> \~japanese
!! \brief NetCDF ファイルのヘッダー部分に, 座標情報を定義する.
!!
!! @param[in]     self                構造型 IcGrid_ncWriter の変数.
!! @param[in,out] coordinate_element 座標情報を保持する構造型 Mesh_dim_element の変数.
!!
!! \~english
!! \brief  Defines the informaion about the coordinate in the header part of the NetCDF file.
!!
!! @param[in]     self                The variable of derived type IcGrid_ncWriter.
!! @param[in,out] coordinate_element The variable of derived type Mesh_dim_element containing the informaion about the coordinate.
!!
subroutine ncdef_mesh_coordinate(self, coordinate_element)

  ! 宣言文 ; Declaration statements
  !
  type(netcdfDataWriter), intent(in) :: self
  type(Mesh_coord_element), intent(inout) :: coordinate_element

  ! 作業変数
  ! Work variables
  !

  ! 実行文 ; Executable statements
  !

  ! 座標変数を定義する.
  call check_nf90_status( &
    & nf90_def_var( &
       & self%ncID, &
       & coordinate_element%element_name, NF90_DOUBLE, &
       & coordinate_element%dimension_element%dimID, coordinate_element%varID &
       & ),  &
       & message = "define mesh coordinate: " // trim(coordinate_element%element_name) &
     & )

  ! 座標変数のメタ情報を加える.
  !

  call check_nf90_status( &
    & nf90_put_att( self%ncID, coordinate_element%varID, 'name', coordinate_element%standard_name ) &
    & )

  call check_nf90_status( &
    & nf90_put_att( self%ncID, coordinate_element%varID, 'long_name', coordinate_element%long_name ) &
    & )

  ! 座標の単位を設定.
  call check_nf90_status( &
    & nf90_put_att( self%ncID, coordinate_element%varID, 'units', coordinate_element%units ) &
    & )

end subroutine ncdef_mesh_coordinate

subroutine ncdef_mesh_connectivity(self, connectivity_element)

  ! 宣言文 ; Declaration statements
  !
  type(netcdfDataWriter), intent(in) :: self
  type(Mesh_connectivity_element), intent(inout) :: connectivity_element

  ! 作業変数
  ! Work variables
  !

  ! 実行文 ; Executable statements
  !


  ! 座標変数を定義する.
  call check_nf90_status( &
    & nf90_def_var( &
       & self%ncID, &
       & connectivity_element%element_name, NF90_INT, &
       & (/ connectivity_element%connected_dim_element%dimId, connectivity_element%own_dim_element%dimID /), &
       & connectivity_element%varID &
       & ), &
       & message = "define mesh connectivity: " // trim(connectivity_element%element_name) &
     & )

  ! 座標変数のメタ情報を加える.
  !

  call check_nf90_status( &
    & nf90_put_att( self%ncID, connectivity_element%varID, 'name', connectivity_element%element_name ) &
    & )

  call check_nf90_status( &
    & nf90_put_att( self%ncID, connectivity_element%varID, 'cf_role', connectivity_element%cf_role ) &
    & )

  call check_nf90_status( &
    & nf90_put_att( self%ncID, connectivity_element%varID, 'long_name', connectivity_element%long_name ) &
    & )

  call check_nf90_status( &
    & nf90_put_att( self%ncID, connectivity_element%varID, 'start_index', connectivity_element%start_index ) &
    & )


end subroutine ncdef_mesh_connectivity

end module netcdfDataWriter_mod
