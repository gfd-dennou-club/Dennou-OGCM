!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module GridSet_mod

  ! モジュール引用; Use statements
  !

  use dc_types, only: &
       & STRING

  use dc_message, only: &
       & MessageNotify
       
  use PolyMesh_mod
  use fvMeshInfo_mod
  use HexTriIcMesh_mod

  !  use SimParameters_mod, only: gridFilePath, Radius

  ! 宣言文 ; Declaration statements  
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedures
  !
  public :: GridSet_Init, GridSet_Final

  ! 公開変数
  ! Public variables
  !
  type(PolyMesh), public, save :: plMesh  
  type(fvMeshInfo), public, save :: fvmInfo
  type(HexTriIcMesh), public, save :: htiMesh
  integer, public, save :: nCell
  integer, public, save :: nEdge
  integer, public, save :: nVertex
  integer, public, save :: nVzLyr
  integer, public, save :: nVrLyr
  integer, parameter, public :: VHaloSize = 1

  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'GridSet_mod' !< Module Name

contains
  subroutine GridSet_Init(configNmlFileName)

    ! モジュール引用; Use statement
    !
    use Constants_mod, only: RPlanet

    use PolyMesh_mod, only: &
         & getCellListSize, getFaceListSize, getPointListSize

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! local variable
    !
    character(STRING) :: gridFilePath


    ! 実行文; Executable statement
    !

    ! Read the path to grid data file from namelist.
    call read_nmlData(configNmlFileName, gridFilePath)

    ! Load grid data of from netcdf file. 
    ! If loading grid data succeed, a polyMesh object whose name is plMesh have read grid data.    
    call constract_grid(gridFilePath)
    call HexTriIcMesh_Init(htiMesh, plMesh, RPlanet)

    ! Initialize some modules for finite volume method
    !
    call MessageNotify( 'M', module_name, "Initialize some modules for finite volume method..")
    call fvMeshInfo_Init(fvmInfo, plMesh, dualMeshFlag=.true.)
    call HexTriIcMesh_configfvMeshInfo(htiMesh, fvmInfo)

    !
    nCell = getCellListSize(plMesh)
    nEdge = getFaceListSize(plMesh)
    nVertex = getPointListSize(plMesh)
    nVrLyr = nVzLyr + 1

  end subroutine GridSet_Init

  subroutine GridSet_Final()

    call fvMeshInfo_Final(fvmInfo)
    call HexTriIcMesh_Final(htiMesh)
    call PolyMesh_Final(plMesh)

  end subroutine GridSet_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_nmlData( configNmlFileName, &
       & gridFilePath )

    ! モジュール引用; Use statement
    !

    ! ファイル入出力補助
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! 種別型パラメタ
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName
    character(STRING), intent(out) :: gridFilePath

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
    ! IOSTAT of NAMELIST read

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /grid_nml/ &
         & gridFilePath, nVzLyr

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    gridFilePath = ""
    nVzLyr    = 1

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = grid_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if


    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  gridFilePath    = %a', ca=(/ gridFilePath /))
    call MessageNotify( 'M', module_name, '  nVzLyr          = %d', i=(/ nVzLyr /))

  end subroutine read_nmlData

  subroutine constract_Grid(gridFilePath)

    ! モジュール引用; Use statement
    !
    use netcdfDataReader_mod, only: &
         & netcdfDataReader, &
         & netcdfDataReader_Init, netcdfDataReader_Final

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: gridFilePath

    ! 局所変数
    ! Local variables
    !
    type(netcdfDataReader) :: ncReader
    type(PolyMesh) :: readPlMesh

    ! 実行文; Executable statement
    !

    call MessageNotify( 'M', module_name, "Load grid data from '%a' ..", ca=(/ gridFilePath /) ) 

    call netcdfDataReader_Init(ncReader, gridFilePath, readPlMesh)
    call netcdfDataReader_Final(ncReader)

    call PolyMesh_Init(plMesh, &
         & readPlMesh%pointPosList, readPlMesh%cellPosList, readPlMesh%faceList, readPlMesh%cellList, nVzLyr)

    call PolyMesh_Final(readPlMesh)

  end subroutine constract_Grid

end module GridSet_mod
