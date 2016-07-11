!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Admin_Grid_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & STRING, TOKEN, DP

  use dc_message, only: &
       & MessageNotify



  ! 宣言文 ; Declaration statements  
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedures
  !
  public :: DOGCM_Admin_Grid_Init, DOGCM_Admin_Grid_Final

  public :: DOGCM_Admin_Grid_construct
  public :: DOGCM_Admin_Grid_UpdateVCoord
  
  ! Cascade


  ! 公開変数
  ! Public variables
  !

  integer, public, save :: IA, IS, IE, IM, IHALO
  integer, public, save :: JA, JS, JE, JM, JHALO
  integer, public, save :: KA, KS, KE, KM, KHALO
  
  integer, public, save :: kMax
  integer, public, save :: tMax
  integer, parameter, public :: VHaloSize = 1

  integer, public, save :: nMax, lMax
  integer, public, save :: iMax, jMax
  integer, public, save :: jMaxGlobe

  
  real(DP), public, allocatable :: xyz_Lon(:,:,:)
  real(DP), public, allocatable :: xyz_Lat(:,:,:)
  real(DP), public, allocatable :: xyz_Z(:,:,:)
  real(DP), public, allocatable :: xy_Topo(:,:)
  
!  real(DP), public, allocatable :: z_Sig(:)  
!  real(DP), public, allocatable :: z_LyrThickSig(:)

  real(DP), public, allocatable :: x_Lon_Weight(:)
  real(DP), public, allocatable :: y_Lat_Weight(:)
  real(DP), public, allocatable :: z_Sig_Weight(:)

  type, public :: AXIS_INFO
     character(TOKEN) :: name
     character(STRING) :: long_name
     character(TOKEN) :: units
     character(TOKEN) :: weight_units
  end type AXIS_INFO

  type(AXIS_INFO), public, save :: IAXIS_info
  type(AXIS_INFO), public, save :: JAXIS_info
  type(AXIS_INFO), public, save :: KAXIS_info
  type(AXIS_INFO), public, save :: TAXIS_info

  real(DP), public, allocatable :: x_CI(:) ! Position of cell center
  real(DP), public, allocatable :: y_CJ(:)
  real(DP), public, allocatable :: z_CK(:)  
  real(DP), public, allocatable :: x_FI(:) ! Position of face
  real(DP), public, allocatable :: y_FJ(:)
  real(DP), public, allocatable :: z_FK(:)  
  real(DP), public, allocatable :: x_IAXIS_Weight(:)
  real(DP), public, allocatable :: y_JAXIS_Weight(:)
  real(DP), public, allocatable :: z_KAXIS_Weight(:)
  
  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'DOGCM_Admin_Grid_mod' !< Module Name

  logical :: isInitialzed = .false.
  
contains
  subroutine DOGCM_Admin_Grid_Init(configNmlFileName)

    ! モジュール引用; Use statement
    !
    use DOGCM_Admin_Constants_mod, only: RPlanet

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

    isInitialzed = .true.
    
  end subroutine DOGCM_Admin_Grid_Init

  subroutine DOGCM_Admin_Grid_Final()

    ! 実行文; Executable statement
    !

    if(isInitialzed) then
       deallocate( xyz_Lon, xyz_Lat, xyz_Z )
       deallocate( x_Lon_Weight, y_Lat_Weight, z_Sig_Weight )
    end if
    
  end subroutine DOGCM_Admin_Grid_Final


  !> @brief 
  !!
  !!
  subroutine DOGCM_Admin_Grid_construct()

    ! モジュール引用; Use statements
    !
    
    ! 宣言文; Declaration statement
    !
    
        
!!$    real(DP) :: A(kMax,0:kMax), B(kMax+1), Work(2*(kMax+1))
!!$    integer :: k, info
    
    ! 実行文; Executable statement
    !

    ! Caclculate grid index
    !


!!$    call DOGCM_Admin_Grid_construct_Cartesian3DGrid()
    
    !
    !
    call DOGCM_Admin_Grid_construct_SPMGaussGrid()
    
  end subroutine DOGCM_Admin_Grid_construct

  subroutine DOGCM_Admin_Grid_construct_SPMGaussGrid()

    ! モジュール引用; Use statement
    !
    
    use SpmlUtil_mod, only: &
        & isSpmlUtilInitialzed=>isInitialzed, &
        & get_SpmlGridInfo

    ! 局所変数
    ! Local variables
    !
    
    real(DP), allocatable :: xy_LonSpml(:,:)
    real(DP), allocatable :: xy_LatSpml(:,:)
    real(DP), allocatable :: z_SigSpml(:)

    real(DP), parameter :: PI = acos(-1d0)
    
    ! 実行文; Executable statement
    !


    if( .not. isSpmlUtilInitialzed() ) &
         & call MessageNotify('E', module_name, &
         &  "DOGCM_Admin_Grid_construct is called before SpmlUtil_mod is initialized.")
    
    IHALO = 1; JHALO = 1; KHALO = 1
    
    IA = IM + 2*IHALO; IS = 1 + IHALO; IE = IA - IHALO
    JA = JM + 2*JHALO; JS = 1 + JHALO; JE = JA - JHALO
    KA = KM + 2*KHALO; KS = 1 + KHALO; KE = KA - KHALO

    call MessageNotify('M', module_name, "(IS,JS,KS)=(%d,%d,%d)", (/ IS, JS, KS /) )
    call MessageNotify('M', module_name, "(IE,JE,KE)=(%d,%d,%d)", (/ IE, JE, KE /) )
    call MessageNotify('M', module_name, "(IHALO,JHALO,KHALO)=(%d,%d,%d)", (/ IHALO, JHALO, KHALO /) )
    
    ! Allocation arrays for coordinates.
    !

    allocate( xyz_Lon(IA,JA,KA), xyz_Lat(IA,JA,KA), xyz_Z(IA,JA,KA) )
    allocate( xy_Topo(IA,JA) )
!!$    allocate(xyz_Lon(0:iMax-1,1:jMax,0:kMax), xyz_Lat(0:iMax-1,1:jMax,0:kMax), z_Sig(0:kMax))

    allocate( x_CI(IA), y_CJ(JA), z_CK(KA) )
    allocate( x_FI(IA), y_FJ(JA), z_FK(KA) )
    allocate( x_IAXIS_Weight(IA), y_JAXIS_Weight(JA), z_KAXIS_Weight(KA) )
    
    ! Get the coordinates of horizontal grid from SpmlUtil_mod.
    !

    allocate( xy_LonSpml(0:iMax-1,jMax), xy_LatSpml(0:iMax-1,jMax), z_SigSpml(0:kMax) )
    allocate( x_Lon_Weight(0:iMax-1), y_Lat_Weight(jMax), z_Sig_Weight(0:kMax) )
    call get_SpmlGridInfo( &
         & xy_Lon_=xy_LonSpml, xy_Lat_=xy_LatSpml, z_Sig_=z_SigSpml,  &
         & x_Lon_Weight_=x_Lon_Weight, y_Lat_Weight_=y_Lat_Weight, z_Sig_Weight_=z_Sig_Weight &
         & )
    
    xyz_Lon(IS:IE,JS:JE,KS:KE) = spread(xy_LonSpml,3,kMax+1)
    xyz_Lat(IS:IE,JS:JE,KS:KE) = spread(xy_LatSpml,3,kMax+1)

    IAXIS_info%name = 'lon'
    IAXIS_info%long_name = 'longitude'
    IAXIS_info%units = 'degree_east'
    IAXIS_info%weight_units = 'radian'

    JAXIS_info%name = 'lat'
    JAXIS_info%long_name = 'latitude'
    JAXIS_info%units = 'degree_north'
    JAXIS_info%weight_units = 'radian'

    KAXIS_info%name = 'sig'
    KAXIS_info%long_name = 'general vertical coordinate'
    KAXIS_info%units = '1'
    KAXIS_info%weight_units = '1'

    TAXIS_info%name = 'time'
    TAXIS_info%long_name = 'time'
    TAXIS_info%units = 'sec'
    TAXIS_info%weight_units = '1'
    
    x_CI(IS:IE) = xy_LonSpml(0:iMax-1,1) * 180d0 / PI
    y_CJ(JS:JE) = xy_LatSpml(0,1:jMax) * 180d0 / PI
    z_CK(KS:KE) = z_SigSpml(0:kMax)
    x_IAXIS_Weight(IS:IE) = x_Lon_Weight
    y_JAXIS_Weight(JS:JE) = y_Lat_Weight
    z_KAXIS_Weight(KS:KE) = z_Sig_Weight
    
    ! Calculate layer thickness.
    !
!    allocate( z_LyrThickSig(0:kMax) )    
    
!!$    A = 0d0; B = 0d0
!!$    do k=1,kMax
!!$       A(k,k-1:k) = 0.5d0
!!$       B(k) = z_Sig(k-1) - z_Sig(k)
!!$    end do
!!$    A(1,0) = 1; A(kMax,kMax) = 1
!!$
!!$    call DGELS('N', kMax, kMax+1, 1, A, kMax, B, kMax+1, Work, 2*(kMax+1), info)
!!$    z_LyrThickSig(:) = b

!    z_LyrThickSig = z_Sig_Weight
    
  end subroutine DOGCM_Admin_Grid_construct_SPMGaussGrid

  subroutine DOGCM_Admin_Grid_UpdateVCoord( xyz_H, &     ! (out)
       & xy_SSH                               ) ! (in)

    ! モジュール引用; Use statement
    !
    
    use SpmlUtil_mod, only: &
         & xyz_DSig_xyz
    
    ! 宣言文; Declaration statement
    !

    real(DP), intent(out) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH(IA,JA)

    ! 局所変数
    ! Local variables
    !
    integer :: k
    
    ! 実行文; Executable statement
    !

    do k = 1, KA
       xyz_Z(:,:,k) = xy_Topo(:,:) * z_CK(k)
    end do
    xyz_H(IS:IE,JS:JE,KS:KE) = xyz_DSig_xyz( xyz_Z(IS:IE,JS:JE,KS:KE) )
    
  end subroutine DOGCM_Admin_Grid_UpdateVCoord
  
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
         & nLat, nLon, nZ, &
         & IM, JM, KM

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

    IM = iMax
    JM = jMax
    KM = kMax + 1
    
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

    call MessageNotify( 'M', module_name, ' (IM,JM,KM) = (%d,%d,%d) ', i = (/ IM,JM,KM /) )
    
  end subroutine read_nmlData


end module DOGCM_Admin_Grid_mod
