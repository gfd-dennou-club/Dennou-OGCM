!> \~japanese
!! \brief IcGrid_ncReader/Writer_mod モジュールで用いられるデータ構造を管理するための構造型を提供するモジュール.
!> \~english
!! \brief This module provides some derived types to manage some data structure used in IcGrid_ncReader/Writer_mod module.
!!
!!
!! <br><br>
!! \b Copyright (C) GFD Dennou Club, 2011-2012. All rights reserved. <br>
!! \b license ?? <br>
!! \author Yuta Kawai
!!
module netcdfDataHelper_mod
  
  ! モジュール引用 ; Use statements
  !

  ! 種類型パラメタ 
  ! Kind type parameter
  !
  use dc_types, only: DP, &        ! 倍精度実数型. Double precision.
    &                 TOKEN, &
    &                 STRING

  ! メッセージ出力
  ! Message dump
  !
  use dc_message, only: &
    & MessageNotify

  ! NetCDF
  !
  !
  use netcdf

  !
  !
  use PolyMesh_mod

  ! 宣言文 ; Declaration statements
  !
  implicit none
  private

  ! 公開変数
  ! Public variables
  !

  !
  !> \~japanese
  !! \brief Mesh の次元情報を管理する構造体.
  !! \~english
  !! \brief The derived type to manage the dimension information about Mesh.
  !!
  type, public :: Mesh_dim_element

    !> \~japanese Mesh の次元の要素名.
    !! \~english
    character(TOKEN) :: element_name

    !> \~japanese Mesh の次元の要素数.
    !! \~english
    integer :: num

    !> \~japanese Mesh の次元データと対応する NetCDF の変数 ID.
    !! \~english
    integer :: dimID

  end type Mesh_dim_element

  !
  !> \~japanese
  !! \brief Mesh の座標情報を管理する構造体.
  !! \~english
  !! \brief The derived type to manage the coordinate information about Mesh.
  !!
  type, public :: Mesh_coord_element

    !> \~japanese Mesh の座標を作る次元要素(構造型 Mesh_dim_element の変数)へのポインタ.
    !! \~english
    type(Mesh_dim_element), pointer :: dimension_element

    !> \~japanese Mesh の座標の要素名.
    !! \~english
    character(TOKEN) :: element_name

    !> \~japanese Mesh の座標の正式名.
    !! \~english
    character(TOKEN) :: standard_name

    !> \~japanese
    !! \~english
    character(STRING) :: long_name

    !> \~japanese Mesh の座標の単位.
    !! \~english
    character(STRING) :: units

    !> \~japanese Mesh の座標データと対応する NetCDF における変数 ID.
    !! \~english
    integer :: varID

  end type Mesh_coord_element
  
  !> \~japanese
  !! \brief Mesh2 のノード, 面, エッジ, 面のリンクに関する情報を管理する構造体.
  !! \~english
  !! \brief The derived type to manage the information about node, face, edge and face link in Mesh2.
  !!
  type, public :: Mesh2_ncInfo

    !> \~japanese Mesh2 の次元要素 node を表す構造型 Mesh_dim_element の変数.
    !! \~english
    type(Mesh_dim_element) :: node

    !> \~japanese Mesh2 の次元要素 face を表す構造型 Mesh_dim_element の変数.
    !! \~english
    type(Mesh_dim_element) :: face

    !> \~japanese
    !! \~english
    type(Mesh_dim_element) :: edge

    !> \~japanese
    !! \~english
    type(Mesh_dim_element) :: face_links
    
    !> \~japanese Mesh2 の座標要素 node_x(次元要素が node である経度方向の座標)を表す構造型 Mesh_coord_element の変数.
    !! \~english
    type(Mesh_coord_element) :: node_x

    !> \~japanese Mesh2 の座標要素 node_y(次元要素が node である緯度方向の座標)を表す構造型 Mesh_coord_element の変数.
    !! \~english
    type(Mesh_coord_element) :: node_y

    !> \~japanese
    !! \~english
    type(Mesh_coord_element) :: face_x, face_y

    !> \~japanese
    !! \~english
    type(Mesh_coord_element) :: edge_x, edge_y
    
  end type Mesh2_ncInfo

  !> \~japanese
  !! \brief Mesh3 のノード, 面, エッジ, 面のリンクに関する情報を管理する構造体.
  !! \~english
  !! \brief The derived type to manage the information about node, face, edge and face link in Mesh2.
  !!
  type, public :: Mesh3_ncInfo

    !> \~japanese
    !! \~english
    type(Mesh_dim_element) :: node

    !> \~japanese
    !! \~english
    type(Mesh_dim_element) :: face

    !> \~japanese
    !! \~english
    type(Mesh_dim_element) :: edge

    !> \~japanese
    !! \~english
    type(Mesh_dim_element) :: face_links

    !> \~japanese
    !! \~english
    type(Mesh_dim_element) :: vlayers

    !> \~japanese
    !! \~english
    type(Mesh_coord_element) :: node_x, node_y

    !> \~japanese
    !! \~english
    type(Mesh_coord_element) :: face_x, face_y

    !> \~japanese
    !! \~english
    type(Mesh_coord_element) :: edge_x, edge_y

    !> \~japanese
    !! \~english
    type(Mesh_coord_element) :: layers

  end type Mesh3_ncInfo

  ! 公開手続き
  ! Public procedure
  !
  public :: Mesh_coord_element_Init, Mesh2_ncInfo_Init, check_nf90_status
  !public :: Mesh3_ncInfo_Init

contains

!
!> \~japanese
!! \brief 座標要素を表現する構造型 Mesh_coord_element の変数を初期化する.
!!
!! @param[in,out] self 構造型 Mesh_coord_element の変数.
!! @param[in,out] dimension_element 座標要素と対応付ける次元要素を表す構造型 Mesh_dim_element の変数.
!! @param[in] element_name 座標要素の名前.
!! @param[in] standard_name 座標要素の正式な名前.
!! @param[in] long_name 座標要素の長い名前.
!! @param[in] units 座標要素の単位.
!!
!! \english
!! \brief Initialize the variable of derived type Mesh_coord_element which represents the coordinate element.
!!
!! @param[in,out] self The variable of derived type Mesh_coord_element.
!! @param[in,out] dimension_element
!! @param[in] element_name
!! @param[in] standard_name
!! @param[in] long_name
!! @param[in] units
!!
!!
subroutine Mesh_coord_element_Init( &
  & self,                             &  ! (inout)
  & dimension_element,                &  ! (inout)
  & element_name, standard_name, long_name, units &  ! (in)
  & )

  ! 宣言文 ; Declaration statement
  !
  type(Mesh_coord_element), intent(inout) :: self
  type(Mesh_dim_element), intent(inout), target :: dimension_element
  character(*), intent(in) :: element_name
  character(*), intent(in) :: standard_name
  character(*), intent(in) :: long_name
  character(*), intent(in) :: units

  ! 実行文 ; Execuatble statement
  !
  self%dimension_element => dimension_element
  self%element_name = element_name
  self%standard_name = standard_name
  self%long_name = long_name
  self%units = units

end subroutine Mesh_coord_element_Init

!
!> \~japanese
!! \brief NetCDF の読み書きの際に, 球面上の水平方向の非構造格子データを管理する構造型 Mesh2_ncInfo の変数を初期化する.
!!
!! @param[in,out] self 構造型 Mesh2_ncInfo の変数.
!! @param[in,out] icgrid 正二十面格子データを管理する構造型 FVM_IcGrid の変数.
!!
!! \~english
!! \brief Initializes the variable of derived type Mesh2_ncInfo which manages the data of unsructured grid horizontally distributed on a sphere
!! when the data is written to or read from NetCDF file.
!!
!! @param[in,out] self The variable of derived type Mesh2_ncInfo.
!! @param[in,out] icgrid The variable of derived type FVM_IcGrid which manges icosahedral grid data.
!!
subroutine Mesh2_ncInfo_Init( &
  & self, mesh &  ! (inout)
  & )
  
  ! モジュール引用 ; Use statements
  !
  

  ! 宣言文 ; Declaration statements
  !
  type(Mesh2_ncInfo), intent(inout) :: self
  type(PolyMesh), intent(in) :: mesh

  ! 作業変数
  ! Work variables
  !
  character(TOKEN), parameter :: LAT = "latitude"
  character(TOKEN), parameter :: LON = "longitude"
  character(TOKEN), parameter :: LAT_UNITS = "degrees_north"
  character(TOKEN), parameter :: LON_UNITS = "degrees_east"

  ! 実行文 ; Executable statement
  !

  !

  !
  self%node%element_name = 'nMesh2_point'
  self%node%num = getPointListSize(mesh)
  self%face%element_name = 'nMesh2_cell '
  self%face%num = getCellListSize(mesh)


  ! Mesh2_node_x の初期化
  call Mesh_coord_element_Init( &
    & self%node_x, dimension_element=self%node, &
    & element_name='Mesh2_node_x', standard_name=LON, &
    & long_name='Latiude of 2D mesh point', units=LON_UNITS &
    & )

  ! Mesh2_node_y の初期化
  call  Mesh_coord_element_Init( &
    & self%node_y, dimension_element=self%node, &
    & element_name='Mesh2_node_y', standard_name=LAT, &
    & long_name='Latiude of 2D mesh point', units=LAT_UNITS &
    & )

  ! Mesh2_node_x の初期化
  call Mesh_coord_element_Init( &
    & self%face_x, dimension_element=self%face, &
    & element_name='Mesh2_face_x', standard_name=LON, &
    & long_name='Latiude of 2D mesh cell', units=LON_UNITS &
    & )

  ! Mesh2_node_y の初期化
  call  Mesh_coord_element_Init( &
    & self%face_y, dimension_element=self%face, &
    & element_name='Mesh2_face_y', standard_name=LAT, &
    & long_name='Latiude of 2D mesh cell', units=LAT_UNITS &
    & )

end subroutine Mesh2_ncInfo_Init

!!$subroutine Mesh3_ncInfo_Init( &
!!$  & self, &
!!$  & icgrid3D &
!!$  & )
!!$
!!$
!!$  ! モジュール引用 ; Use statement
!!$  !
!!$
!!$  use IcGrid2D_FVM_Manager, only: &
!!$    & IcGrid2D_FVM, &
!!$    & get_EffSize_Min, get_EffSize_Max, &
!!$    & RC_REGIONS_NUM
!!$
!!$  use IcGrid3D_FVM_Manager, only: &
!!$    & IcGrid3D_FVM, &
!!$    & get_icgrid2D, get_vertical_level_num, &
!!$    & get_IcGrid3D_VerticalCoord_name, &
!!$    & get_IcGrid3D_VerticalCoord_long_name, &
!!$    & get_IcGrid3D_VerticalCoord_units
!!$
!!$  ! 宣言文 ; Declaration statement
!!$  !
!!$  type(Mesh3_ncInfo), intent(inout) :: self
!!$  type(IcGrid3D_FVM), intent(inout) :: icgrid3D
!!$
!!$  ! 作業変数
!!$  ! Work variable
!!$  !
!!$  integer :: EMin, Emax
!!$  character(TOKEN), parameter :: LAT = "latitude"
!!$  character(TOKEN), parameter :: LON = "longitude"
!!$  character(TOKEN), parameter :: LAT_UNITS = "degrees_north"
!!$  character(TOKEN), parameter :: LON_UNITS = "degrees_east"
!!$
!!$  ! 実行文 ; Executable statement
!!$  !
!!$
!!$  !
!!$  EMin = get_EffSize_Min(get_icgrid2D(icgrid3D))
!!$  EMax = get_EffSize_Max(get_icgrid2D(icgrid3D))
!!$
!!$  !
!!$  !
!!$
!!$  !
!!$  self%node%element_name = 'nMesh3_node'
!!$  ! dimension: Mesh3_node 個数の設定
!!$  self%node%num = RC_REGIONS_NUM * ( EMax - EMin + 1)**2
!!$
!!$  ! Mesh3_node_x の初期化
!!$  call Mesh_coord_element_Init( &
!!$    & self%node_x, dimension_element=self%node, &
!!$    & element_name='Mesh3_node_x', standard_name=LON, &
!!$    & long_name='Latiude of 3D mesh node', units=LON_UNITS &
!!$    & )
!!$
!!$  ! Mesh3_node_y の初期化
!!$  call  Mesh_coord_element_Init( &
!!$    & self%node_y, dimension_element=self%node, &
!!$    & element_name='Mesh3_node_y', standard_name=LAT, &
!!$    & long_name='Latiude of 3D mesh node', units=LAT_UNITS &
!!$    & )
!!$
!!$  !
!!$  !
!!$
!!$  !
!!$  self%vlayers%element_name = 'nMesh3_layers'
!!$  ! dimension: Mesh3_node 個数の設定
!!$  self%vlayers%num = get_vertical_level_num(icgrid3D)
!!$
!!$  ! Mesh3_layers の初期化
!!$  call  Mesh_coord_element_Init( &
!!$    & self%layers, dimension_element=self%vlayers, &
!!$    & element_name='Mesh3_layers', &
!!$    & standard_name=get_IcGrid3D_VerticalCoord_name(icgrid3D), &
!!$    & long_name=get_IcGrid3D_VerticalCoord_long_name(icgrid3D), &
!!$    & units=get_IcGrid3D_VerticalCoord_units(icgrid3D) &
!!$    & )
!!$
!!$end subroutine Mesh3_ncInfo_Init

!
!> \~japanese
!! \brief NetCDF ライブラリが提供する関数の終了ステータスを調べる.
!!
!! netCDF ライブラリが提供する手続き実行後に返されるステータスを解析し,
!! エラーがあった場合はエラー出力を行い, プログラムを停止する.
!!
!! @param[in] status  netCDF ライブラリが提供する手続きが返すステータス
!! @param[in] message エラーがあった場合に出力するメッセージ
!!
!! \~english
!! \brief Checks the exit status of a function provided by NetCDF library.
!!
!! netCDF ライブラリが提供する手続き実行後に返されるステータスを解析し,
!! エラーがあった場合はエラー出力を行い, プログラムを停止する. 
!!
!! @param[in] status  netCDF ライブラリが提供する手続きが返すステータス
!! @param[in] message エラーがあった場合に出力するメッセージ
!!
subroutine check_nf90_status( &
  & status, message &  ! (in)
  & )
  
  ! 宣言文 ; Declaration statement
  !
  integer, intent(in) :: status 
  character(*), intent(in), optional :: message

  ! 作業変数
  ! Work varibles
  !
  
  character(STRING) :: output_str
  character(STRING) :: netCDF_error 
 
  ! 実行文 ; Execuatble statements
  !

  if ( status /= nf90_noerr ) then
    ! ステータスを解析し, エラーの原因を取得する. 
    netCDF_error = '(netCDF error status:' // trim(nf90_strerror(status)) // ')'
 
    if ( present(message) ) then
      output_str = trim(message) // netCDF_error
    else
      output_str = netCDF_error
    end if

    ! エラー出力を行い, プログラムを停止する. 
!    write(*,*) trim(output_str)
    call MessageNotify( 'E', 'IGMBaseLib IO:IcGrid_ncStream_helper:', trim(output_str) ) 
    stop

  end if

end subroutine check_nf90_status

end module NetcdfDataHelper_mod
