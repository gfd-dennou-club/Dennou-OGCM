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

  !* Dennou-OGCM

  use DOGCM_Admin_GaussSpmGrid_mod, only: &
       & DOGCM_Admin_GaussSpmGrid_Init,               &
       & DOGCM_Admin_GaussSpmGrid_Final,              &
       & DOGCM_Admin_GaussSpmGrid_ConstructAxisInfo,  &
       & DOGCM_Admin_GaussSpmGrid_ConstructGrid,      &
       & iMax, jMax, jMaxGlobe, kMax,                 &
       & nMax, lMax, tMax

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

  integer, public :: IA, IS, IE, IM, IHALO
  integer, public :: JA, JS, JE, JM, JHALO
  integer, public :: KA, KS, KE, KM, KHALO
  
  integer, parameter, public :: VHaloSize = 1

  real(DP), public, allocatable :: xy_Lon(:,:)
  real(DP), public, allocatable :: xy_Lat(:,:)
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
  real(DP), public, allocatable :: x_CDI(:) 
  real(DP), public, allocatable :: y_CDJ(:)
  real(DP), public, allocatable :: z_CDK(:)  
  real(DP), public, allocatable :: x_FDI(:) ! 
  real(DP), public, allocatable :: y_FDJ(:)
  real(DP), public, allocatable :: z_FDK(:)  
  real(DP), public, allocatable :: x_IAXIS_Weight(:)
  real(DP), public, allocatable :: y_JAXIS_Weight(:)
  real(DP), public, allocatable :: z_KAXIS_Weight(:)

  real(DP), public, allocatable :: SCALEF_E1(:,:,:)
  real(DP), public, allocatable :: SCALEF_E2(:,:,:)
  real(DP), public, allocatable :: SCALEF_E3(:,:,:)
  
  ! Cascade

  public :: nMax, lMax, tMax
  public :: iMax, jMax, jMaxGlobe, kMax
  
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

    !---------------------------------------
    
    call DOGCM_Admin_GaussSpmGrid_ConstructAxisInfo( &
       & IAXIS_info%name, IAXIS_info%long_name, IAXIS_info%units, IAXIS_info%weight_units,    & ! (out)
       & IS, IE, IA, IHALO,                                                                   & ! (out)
       & JAXIS_info%name, JAXIS_info%long_name, JAXIS_info%units, JAXIS_info%weight_units,    & ! (out)
       & JS, JE, JA, JHALO,                                                                   & ! (out)
       & KAXIS_info%name, KAXIS_info%long_name, KAXIS_info%units, KAXIS_info%weight_units,    & ! (out)
       & KS, KE, KA, KHALO,                                                                   & ! (out)
       & IM, JM, KM                                                                           & ! (in)
       & )

    call MessageNotify('M', module_name, "(IS,JS,KS)=(%d,%d,%d)", (/ IS, JS, KS /) )
    call MessageNotify('M', module_name, "(IE,JE,KE)=(%d,%d,%d)", (/ IE, JE, KE /) )
    call MessageNotify('M', module_name, "(IHALO,JHALO,KHALO)=(%d,%d,%d)", (/ IHALO, JHALO, KHALO /) )

    TAXIS_info%name      = "time"
    TAXIS_info%long_name = "time"
    TAXIS_info%units     = "sec"
    
  end subroutine DOGCM_Admin_Grid_Init

  subroutine DOGCM_Admin_Grid_Final()

    ! 実行文; Executable statement
    !

    if(isInitialzed) then
       deallocate( x_CI, x_CDI, x_FI, x_FDI )
       deallocate( y_CJ, y_CDJ, y_FJ, y_FDJ )
       deallocate( z_CK, z_CDK, z_FK, z_FDK )
       deallocate( x_IAXIS_Weight, y_JAXIS_Weight, z_KAXIS_Weight )
       
       deallocate( xy_Lon, xy_Lat )
       deallocate( xyz_Lon, xyz_Lat, xyz_Z )
       deallocate( xy_Topo )
       deallocate( SCALEF_E1, SCALEF_E2 )
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
        
    ! 局所変数
    ! Local variables
    !
    integer :: k

    ! 実行文; Executable statement
    !

    !---------------------------------
    
    
    allocate( x_CI(IA), x_CDI(IA), x_FI(IA), x_FDI(IA), x_IAXIS_Weight(IA) )
    allocate( y_CJ(JA), y_CDJ(JA), y_FJ(JA), y_FDJ(JA), y_JAXIS_Weight(JA) )
    allocate( z_CK(KA), z_CDK(KA), z_FK(KA), z_FDK(KA), z_KAXIS_Weight(KA) )

    allocate( xy_Lon(IA,JA), xy_Lat(IA,JA) )
    allocate( xyz_Lon(IA,JA,KA), xyz_Lat(IA,JA,KA), xyz_Z(IA,JA,KA) )
    allocate( xy_Topo(IA,JA) )
    allocate( SCALEF_E1(IA,JA,4), SCALEF_E2(IA,JA,4) )

    call DOGCM_Admin_GaussSpmGrid_ConstructGrid( &
         & x_CI, x_CDI, x_FI, x_FDI, x_IAXIS_Weight,      & ! (out)
         & y_CJ, y_CDJ, y_FJ, y_FDJ, y_JAXIS_Weight,      & ! (out)
         & z_CK, z_CDK, z_FK, z_FDK, z_KAXIS_Weight,      & ! (out)
         & xy_Lon, xy_Lat,                                & ! (out)
         & SCALEF_E1, SCALEF_E2,                          & ! (out)
         & IS, IE, IA, IM, IHALO,                         & ! (in)
         & JS, JE, JA, JM, JHALO,                         & ! (in)
         & KS, KE, KA, KM, KHALO                          & ! (in)
         & ) 

    do k = 1, KA
       xyz_Lon(:,:,k) = xy_Lon(:,:)
       xyz_Lat(:,:,k) = xy_Lat(:,:)
    end do
    
    !---------------------------------------------

    isInitialzed = .true.
    
  end subroutine DOGCM_Admin_Grid_construct

  subroutine DOGCM_Admin_Grid_UpdateVCoord( xyz_H, &     ! (out)
       & xy_SSH                               )          ! (in)

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
  
!--------------------------------------------------------
    
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
         & IM, JM, KM

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    IM = 1
    JM = 64
    KM = 61
   
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
    call MessageNotify( 'M', module_name, ' (IM,JM,KM) = (%d,%d,%d) ', i = (/ IM,JM,KM /) )
    
  end subroutine read_nmlData


end module DOGCM_Admin_Grid_mod
