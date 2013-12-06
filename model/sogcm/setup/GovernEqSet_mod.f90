!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module GovernEqSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & TOKEN, STRING, DP

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: GovernEqSet_Init, GovernEqSet_Final

  ! 公開変数
  ! Public variable
  !

  integer, public, save :: EOSType
  integer, public, parameter :: EOSTYPE_LINEAR = 1
  integer, public, parameter :: EOSTYPE_JM95   = 2
  character(*), public, parameter :: EOSTYPENAME_LINEAR = 'EOS_LINEAR'
  character(*), public, parameter :: EOSTYPENAME_JM95   = 'EOS_JM95'

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'GovernEqSet_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine GovernEqSet_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    call read_nmlData(configNmlName)

  end subroutine GovernEqSet_Init

  !>
  !!
  !!
  subroutine GovernEqSet_Final()

    ! 実行文; Executable statements
    !

  end subroutine GovernEqSet_Final


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

    character(TOKEN) :: EOSTypeName

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /governEq_nml/ &
         & EOSTypeName

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    EOSTypeName = EOSTYPENAME_LINEAR

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                                  ! (in)
            & nml = governEq_nml, iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    select case(EOSTypeName)
       case (EOSTYPENAME_LINEAR)
          EOSType = EOSTYPE_LINEAR
       case (EOSTYPENAME_JM95)
          EOSType = EOSTYPE_JM95
       case default
          call MessageNotify('E', module_name, &
               & 'The specified EOSType ''%c'' is invalid.', c1=EOSTypeName )
    end select

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '    EOSType      = %c', c1 = EOSTypeName )

  end subroutine read_nmlData

end module GovernEqSet_mod

