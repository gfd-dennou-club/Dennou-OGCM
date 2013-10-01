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
  public :: GridSet_getLocalMesh, GridSet_getLocalMeshInfo
  public :: GridSet_getLocalFVMInfo

  ! 公開変数
  ! Public variables
  !
  type(PolyMesh), public, save :: globalMesh
  type(HexTriIcMesh), public, save :: htiMesh
  type(fvMeshInfo), public, save, allocatable, target :: fvmInfos(:)

  integer, public, save :: nLocalPMesh
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
    integer :: meshID
    type(HexTriIcLocalMesh), pointer :: localMesh

    ! 実行文; Executable statement
    !

    ! Read the path to grid data file from namelist.
    call read_nmlData(configNmlFileName, gridFilePath)

    ! Load grid data of from netcdf file. 
    ! If loading grid data succeed, a polyMesh object whose name is globalMesh have read grid data.    
    call constract_grid(gridFilePath)

    nVrLyr = nVzLyr + 1
    nLocalPMesh = 1
    call HexTriIcMesh_Init(htiMesh, globalMesh, RPlanet, nLocalPMesh)

    do meshID=1, nLocalPMesh
       localMesh => HexTriIcMesh_getLocalMesh(htiMesh, meshID)
       localMesh%localMeshId = meshId 
       call PolyMesh_Init(localMesh%mesh, &
            & globalMesh%pointPosList, globalMesh%cellPosList, globalMesh%faceList, globalMesh%cellList, nVzLyr)
    end do

    ! Initialize some modules for finite volume method
    !
    call MessageNotify( 'M', module_name, "Initialize some modules for finite volume method..")

    allocate( fvmInfos(nLocalPMesh) )
    do meshID=1, nLocalPMesh
       localMesh => HexTriIcMesh_getLocalMesh(htiMesh, meshID)
       call fvMeshInfo_Init(fvmInfos(meshId), localMesh%mesh, dualMeshFlag=.true.)
       call HexTriIcMesh_configfvMeshInfo(htiMesh, fvmInfos(meshId), meshId)
    end do

    !

  end subroutine GridSet_Init

  subroutine GridSet_Final()

    integer :: meshID
    type(HexTriIcLocalMesh), pointer :: localMesh

    
    do meshID=1, nLocalPMesh
       localMesh => HexTriIcMesh_getLocalMesh(htiMesh, meshID)
       call PolyMesh_Final(localMesh%mesh)
       call fvMeshInfo_Final(fvmInfos(meshID))
    end do
    deallocate(fvmInfos)

    call HexTriIcMesh_Final(htiMesh)
    call PolyMesh_Final(globalMesh)
    
  end subroutine GridSet_Final

  function GridSet_getLocalMesh( localMeshId ) result(plmesh)
    integer, intent(in) :: localMeshId
    type(HexTriIcLocalMesh), pointer :: lcMesh
    type(PolyMesh), pointer :: plMesh

    lcMesh => HexTriIcMesh_getLocalMesh(htiMesh, localMeshId)
    plMesh => lcMesh%mesh

  end function GridSet_getLocalMesh

  function GridSet_getLocalFVMInfo( localMeshId ) result(fvm)
    integer, intent(in) :: localMeshId

    type(fvMeshInfo), pointer :: fvm

    fvm => fvmInfos(localMeshId)

  end function GridSet_getLocalFVMInfo

  pure subroutine GridSet_getLocalMeshInfo( &
       & mesh, nCell, nEdge, nVertex)

    type(PolyMesh), intent(in) :: mesh
    integer, intent(out), optional :: nCell, nEdge, nVertex

    if(present(nCell))   nCell = getCellListSize(mesh)
    if(present(nEdge))   nEdge = getFaceListSize(mesh)
    if(present(nVertex)) nVertex = getPointListSize(mesh)

  end subroutine GridSet_getLocalMeshInfo

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

    call PolyMesh_Init(globalMesh, &
         & readPlMesh%pointPosList, readPlMesh%cellPosList, readPlMesh%faceList, readPlMesh%cellList, nVzLyr)

    call PolyMesh_Final(readPlMesh)

  end subroutine constract_Grid

end module GridSet_mod
