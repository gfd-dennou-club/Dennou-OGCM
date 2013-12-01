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
       & STRING, DP

  use dc_message, only: &
       & MessageNotify
       
  use wa_module, only: & 
       & x_Lon, y_Lat

  !  use SimParameters_mod, only: gridFilePath, Radius

  ! 宣言文 ; Declaration statements  
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedures
  !
  public :: GridSet_Init, GridSet_Final
  public :: GridSet_construct

  ! Cascade
  public :: x_Lon, y_Lat

  ! 公開変数
  ! Public variables
  !
  integer, public, save :: kMax
  integer, public, save :: tMax
  integer, parameter, public :: VHaloSize = 1

  integer, public, save :: nMax, lMax
  integer, public, save :: iMax, jMax
  integer, public, save :: jMaxGlobe
  integer, public, save :: nLon

  real(DP), public, save, allocatable :: xyz_Lon(:,:,:)
  real(DP), public, save, allocatable :: xyz_Lat(:,:,:)

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

    ! 実行文; Executable statement
    !

    ! Read the path to grid data file from namelist.
    call read_nmlData(configNmlFileName)


  end subroutine GridSet_Init

  subroutine GridSet_Final()

    ! 実行文; Executable statement
    !

    if(allocated(xyz_Lat)) &
         & deallocate(xyz_Lat, xyz_Lon)

  end subroutine GridSet_Final


  !> @brief 
  !!
  !!
  subroutine GridSet_construct()

    !
    !
    use SpmlUtil_mod, only: &
        & isSpmlUtilInitialzed=>isInitialzed, &
        & xy_Lon, xy_Lat
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    if( .not. isSpmlUtilInitialzed() ) &
         & call MessageNotify('E', module_name, &
         &  "GridSet_construct is called before SpmlUtil_mod is initialized.")

    allocate(xyz_Lon(0:iMax-1,1:jMax,0:kMax))
    allocate(xyz_Lat(0:iMax-1,1:jMax,0:kMax))

    xyz_Lon = spread(xy_Lon,3,kMax+1)
    xyz_Lat = spread(xy_Lat,3,kMax+1)

  end subroutine GridSet_construct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine read_nmlData( configNmlFileName )

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

    ! 局所変数
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
    ! IOSTAT of NAMELIST read

    integer :: nLat, nLon, nZ

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /grid_nml/ &
         & nLat, nLon, nZ

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

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

    iMax = nLon
    jMaxGlobe = nLat
    jMax = jMaxGlobe
    kMax = nZ

    nMax = ( iMax - 1 )/ 3
    lMax = ( nMax + 1 )**2
    tMax = kMax

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '    nmax        = %d', i = (/   nmax  /) )
    call MessageNotify( 'M', module_name, '    iMax        = %d', i = (/   iMax  /) )
    call MessageNotify( 'M', module_name, '    jMaxGlobe = %d', i = (/   jMaxGlobe /) )
    call MessageNotify( 'M', module_name, '    jMax        = %d', i = (/   jMax  /) )
    call MessageNotify( 'M', module_name, '    kMax        = %d', i = (/   kmax  /) )
    call MessageNotify( 'M', module_name, '    lMax        = %d', i = (/   lmax  /) )
    call MessageNotify( 'M', module_name, '    tMax        = %d', i = (/   tMax  /) )


  end subroutine read_nmlData


end module GridSet_mod
