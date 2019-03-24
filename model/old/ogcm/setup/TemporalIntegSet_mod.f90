!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module TemporalIntegSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: TemporalIntegSet_Init, TemporalIntegSet_Final
  public :: TemporalIntegSet_AdvanceLongTStep
  public :: TemporalIntegSet_AdvanceShortTStep

  ! 公開変数
  ! Public variable
  !
  integer, parameter, public :: nLongTimeLevel = 3
  integer, save, public :: Bl, Nl, Al

  integer, parameter, public :: nShortTimeLevel = 3
  integer, save, public :: Bs, Ns, As

  real(DP), save, public :: DelTime
  integer, save, public :: SubCycleNum
  real(DP), save, public :: CurrentTime
  real(DP), save, public :: StartTime
  real(DP), save, public :: TotalIntegTime
  integer, save, public :: CurrentTimeStep
  integer, save, public :: CurrentShortTimeStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'TemporalIntegSet_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine TemporalIntegSet_Init(configNmlFileName)

    ! モジュール引用; Use statement
    !


    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 実行文; Executable statements
    !

    ! Set DelTime, SubCycleNum by reading namelist. 
    call read_nmlData( configNmlFileName  )

    currentTime = StartTime
    CurrentTimeStep = 1
    CurrentShortTimeStep = 1

    Bl = 1; Nl = 2; Al = 3
    Bs = 1; Ns = 2; As = 3

  end subroutine TemporalIntegSet_Init

  !>
  !!
  !!
  subroutine TemporalIntegSet_Final()

    ! 実行文; Executable statements
    !

  end subroutine TemporalIntegSet_Final

  !> @brief 
  !!
  !!
  subroutine TemporalIntegSet_AdvanceLongTStep()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    integer :: oldBl
    
    ! 実行文; Executable statement
    !
 
    currentTime = currentTime + delTime
    CurrentTimeStep = CurrentTimeStep + 1
    CurrentShortTimeStep = 1

    oldBl = Bl
    Bl = Nl
    Nl = Al
    Al = oldBl
    
    call MessageNotify( 'M', module_name, "Advance time step.. current time=%d [sec]", i=(/ int(CurrentTime) /)) 
  end subroutine TemporalIntegSet_AdvanceLongTStep

  subroutine TemporalIntegSet_AdvanceShortTStep()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    integer :: oldBs    
    
    ! 実行文; Executable statement
    !

    oldBs = Bs
    Bs = Ns
    Ns = As
    As = oldBs

    CurrentShortTimeStep = CurrentShortTimeStep + 1

  end subroutine TemporalIntegSet_AdvanceShortTStep


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /temporalInteg_nml/ &
         & DelTIme, SubCycleNum, StartTime, TotalIntegTime

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    DelTIme     = 100d0
    SubCycleNum = 2
    StartTime  = 0d0
    TotalIntegTime  = DelTime

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &          ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                  ! (in)
            & nml = temporalInteg_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if


    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '  DelTime              = %f [sec]', d=(/ DelTime /))
    call MessageNotify( 'M', module_name, '  SubCycleNum          = %d [time]', i=(/ SubCycleNum /))
    call MessageNotify( 'M', module_name, '  StartTime            = %f [sec]', d=(/ StartTime /))
    call MessageNotify( 'M', module_name, '  TotalIntegTime       = %f [sec]', d=(/ TotalIntegTime /))

  end subroutine read_nmlData

end module TemporalIntegSet_mod

