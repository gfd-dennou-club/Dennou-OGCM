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
       & STRING, TOKEN, DP

  use dc_message, only: &
       & MessageNotify

       
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

  real(DP), public, allocatable :: xyz_Lon(:,:,:)
  real(DP), public, allocatable :: xyz_Lat(:,:,:)
  real(DP), public, allocatable :: z_LyrThickSig(:)

  
  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'GridSet_mod' !< Module Name

contains
  subroutine GridSet_Init(configNmlFileName)

    ! モジュール引用; Use statement
    !
    use Constants_mod, only: RPlanet

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
         & deallocate(xyz_Lat, xyz_Lon, z_LyrThickSig)

  end subroutine GridSet_Final


  !> @brief 
  !!
  !!
  subroutine GridSet_construct()

    ! モジュール引用; Use statements
    !
    use SpmlUtil_mod, only: &
        & isSpmlUtilInitialzed=>isInitialzed, &
        & get_HorizontalGrid, g_Sig
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    real(DP), dimension(0:iMax-1,jMax) :: xy_Lon, xy_Lat

    real(DP) :: A(kMax,0:kMax), B(kMax+1), Work(2*(kMax+1))
    integer :: k, info
    
    ! 実行文; Executable statement
    !

    
    if( .not. isSpmlUtilInitialzed() ) &
         & call MessageNotify('E', module_name, &
         &  "GridSet_construct is called before SpmlUtil_mod is initialized.")
    

    ! Allocation
    !
    
    allocate(xyz_Lon(0:iMax-1,1:jMax,0:kMax))
    allocate(xyz_Lat(0:iMax-1,1:jMax,0:kMax))
    allocate(z_LyrThickSig(0:kMax))

    ! Get the coordinates of horizontal grid from SpmlUtil_mod.
    !
    
    call get_HorizontalGrid(xy_Lon_=xy_Lon, xy_Lat_=xy_Lat)

    xyz_Lon(:,:,:) = spread(xy_Lon,3,kMax+1)
    xyz_Lat(:,:,:) = spread(xy_Lat,3,kMax+1)

    ! Calculate layer thickness.
    !
    
    A = 0d0; B = 0d0
    do k=1,kMax
       A(k,k-1:k) = 0.5d0
       B(k) = g_Sig(k-1) - g_Sig(k)
    end do
    A(1,0) = 1; A(kMax,kMax) = 1

    call DGELS('N', kMax, kMax+1, 1, A, kMax, B, kMax+1, Work, 2*(kMax+1), info)
    z_LyrThickSig(:) = b

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

#ifdef DSOGCM_MODE_AXISYM
    if (iMax /= 1) then
       call MessageNotify("E", module_name, &
            & "The number of grid points in longitude must be 1 in DSOGCM_MODE_AXISYM")
    end if
#endif

    ! Set truncated wave number
#ifdef DSOGCM_MODE_AXISYM
    nMax = ( 2*jMaxGlobe - 1 )/ 3
    lMax = nMax + 1 
#else
    nMax = ( iMax - 1 )/ 3
    lMax = ( nMax + 1 )**2
#endif

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
