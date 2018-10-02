!-------------------------------------------------------------
! Copyright (c) 2016-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module GridIndex_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & STRING, TOKEN, DP

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  
  ! 宣言文 ; Declaration statements  
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedures
  !
  public :: GridIndex_Init, GridIndex_Final
  
  ! Cascade


  ! 公開変数
  ! Public variables
  !

  integer, public :: IA, IS, IE, IM, IHALO, IBLOCK
  integer, public :: JA, JS, JE, JM, JHALO, JBLOCK
  integer, public :: KA, KS, KE, KM, KHALO, KBLOCK
  
  integer, parameter, public :: VHaloSize = 1

  integer, public :: iMax
  integer, public :: jMax
  integer, public :: jMaxGlobal
  integer, public :: kMax
  integer, public :: lMax
  integer, public :: nMax
  integer, public :: tMax

  integer, public, parameter :: XDIR = 1
  integer, public, parameter :: YDIR = 2
  integer, public, parameter :: ZDIR = 3

  integer, public, parameter :: I_XY = 1
  integer, public, parameter :: I_UY = 2
  integer, public, parameter :: I_XV = 3
  integer, public, parameter :: I_UV = 4
  
  ! Cascade

  
  ! 非公開変数
  ! Private variables
  !

  character(*), parameter:: module_name = 'GridIndex_mod' !< Module Name

  logical :: isInitialzed = .false.
  
contains
  subroutine GridIndex_Init(configNmlFileName)

    ! モジュール引用; Use statement
    !
    
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

  end subroutine GridIndex_Init

  subroutine GridIndex_Final()

    ! 実行文; Executable statement
    !
    
  end subroutine GridIndex_Final


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
         & IM, JM, KM,             &
         & IBLOCK, JBLOCK, KBLOCK

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    IM = 1
    JM = 64
    KM = 61
    
    IBLOCK = -1
    JBLOCK = -1
    KBLOCK = -1

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

    if(IBLOCK == -1) IBLOCK = IM
    if(JBLOCK == -1) JBLOCK = JM
    if(KBLOCK == -1) KBLOCK = KM
    
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, ' (IM,JM,KM) = (%d,%d,%d) ', i = (/ IM,JM,KM /) )
    call MessageNotify( 'M', module_name, ' (IBLCOK,JBLOCK,KBLOCK) = (%d,%d,%d) ', &
         &                                                 i = (/ IBLOCK,JBLOCK,KBLOCK /) )
    
  end subroutine read_nmlData


end module GridIndex_mod
