!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_Admin_Grid_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  
  use DSIce_Admin_GaussSpmGrid_mod, only: &
       & DSIce_Admin_GaussSpmGrid_Init, DSIce_Admin_GaussSpmGrid_Final, &
       & DSIce_Admin_GaussSpmGrid_ConstructAxisInfo,              &
       & DSIce_Admin_GaussSpmGrid_ConstructGrid
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DSIce_Admin_Grid_Init, DSIce_Admin_Grid_Final
  public :: DSIce_Admin_Grid_Construct
  
  ! 公開変数
  ! Public variables
  !
  
  integer, public :: IM
  integer, public :: IA
  integer, public :: IS
  integer, public :: IE
  integer, public :: IHALO
  
  integer, public :: JM  
  integer, public :: JA  
  integer, public :: JS
  integer, public :: JE
  integer, public :: JHALO

  integer, public :: KM  
  integer, public :: KA  
  integer, public :: KS
  integer, public :: KE
  integer, public :: KHALO
  

  real(DP), allocatable :: CX(:)
  real(DP), allocatable :: CDX(:)
  real(DP), allocatable :: FX(:)
  real(DP), allocatable :: FDX(:)

  real(DP), allocatable :: CY(:)
  real(DP), allocatable :: CDY(:)
  real(DP), allocatable :: FY(:)
  real(DP), allocatable :: FDY(:)

  real(DP), allocatable :: CK(:)
  real(DP), allocatable :: CDK(:)
  real(DP), allocatable :: FK(:)
  real(DP), allocatable :: FDK(:)
  
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

  real(DP), public, allocatable :: xy_Lon(:,:)
  real(DP), public, allocatable :: xy_Lat(:,:)

  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_Admin_Grid_mod' !< Module Name

  logical :: initedFlag = .false.
  
contains

  !>
  !!
  !!
  subroutine DSIce_Admin_Grid_Init( configNmlName )

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName
    
    ! 実行文; Executable statements
    !

    call read_nmlData( configNmlName )
        
    
  end subroutine DSIce_Admin_Grid_Init

  !>
  !!
  !!
  subroutine DSIce_Admin_Grid_Final()

    ! 実行文; Executable statements
    !

    if ( initedFlag  ) then
       deallocate( x_CI, x_CDI, x_FI, x_FDI )
       deallocate( y_CJ, y_CDJ, y_FJ, y_FDJ )
       deallocate( z_CK, z_CDK, z_FK, z_FDK )

       deallocate( xy_Lon, xy_Lat )
    end if
    
  end subroutine DSIce_Admin_Grid_Final

  !------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine DSIce_Admin_Grid_Construct()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !


    call DSIce_Admin_GaussSpmGrid_ConstructAxisInfo( &
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
    
    allocate( x_CI(IA), x_CDI(IA), x_FI(IA), x_FDI(IA), x_IAXIS_Weight(IA) )
    allocate( y_CJ(JA), y_CDJ(JA), y_FJ(JA), y_FDJ(JA), y_JAXIS_Weight(JA) )
    allocate( z_CK(KA), z_CDK(KA), z_FK(KA), z_FDK(KA), z_KAXIS_Weight(KA) )

    allocate( xy_Lon(IA,JA), xy_Lat(IA,JA) )

    
    call DSIce_Admin_GaussSpmGrid_ConstructGrid( &
       & x_CI, x_FI, x_IAXIS_Weight,                    & ! (out)
       & y_CJ, y_FJ, y_JAXIS_Weight,                    & ! (out)
       & z_CK, z_FK, z_KAXIS_Weight,                    & ! (out)
       & xy_Lon, xy_Lat,                                & ! (out)
       & IS, IE, IA, IM, JS, JE, JA, JM, KS, KE, KA, KM & ! (in)
       & ) 


    initedFlag = .true.
    
  end subroutine DSIce_Admin_Grid_Construct

  !-------------------------------------------------
  
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
    namelist /SeaIce_grid_nml/ &
         & IM, JM, KM

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    IM = 1
    JM = 32
    KM = 2
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                  ! (in)
            & nml = seaice_grid_nml,   &  ! (out)
            & iostat = iostat_nml )       ! (out)
       close( unit_nml )
    end if

  end subroutine read_nmlData
  
end module DSIce_Admin_Grid_mod

