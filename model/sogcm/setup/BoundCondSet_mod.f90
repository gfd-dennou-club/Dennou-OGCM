!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module BoundCondSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & TOKEN

  use dc_message,only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: BoundCondSet_Init, BoundCondSet_Final

  ! 非公開手続き
  ! Private procedure
  !
  integer, public, save :: KinBC_Surface
  integer, public, save :: DynBC_Surface
  integer, public, save :: DynBC_Bottom
  integer, public, parameter :: DynBCTYPE_NoSlip = 1
  integer, public, parameter :: DynBCTYPE_Slip = 2
  integer, public, parameter :: KinBCTYPE_FreeSurf = 3
  integer, public, parameter :: KinBCTYPE_RigidLid = 4

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'BoundCondSet_mod' !< Module Name


contains

  !>
  !!
  !!
  subroutine BoundCondSet_Init(configNmlFileName)

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


  end subroutine BoundCondSet_Init

  !>
  !!
  !!
  subroutine BoundCondSet_Final()

    ! 実行文; Executable statements
    !

  end subroutine BoundCondSet_Final


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

    character(TOKEN) :: KinBCSurface, DynBCSurface, DynBCBottom

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /boundaryCondition_nml/ &
         & KinBCSurface, DynBCSurface, DynBCBottom

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
            & nml = boundaryCondition_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    select case(KinBCSurface)
       case("Free")
          KinBC_Surface = KinBCTYPE_FreeSurf
       case("Rigid")
          KinBC_Surface = KinBCTYPE_RigidLid
       case default
          call MessageNotify("E", module_name, &
               & "The Kinetic boundary condition '%c' imposed in the surface is invalid.", c1=KinBCSurface) 
    end select

    select case(DynBCSurface)
       case("Slip")
          DynBC_Surface = DynBCTYPE_Slip
       case("NoSlip")
          DynBC_Surface = DynBCTYPE_NoSlip
       case default
          call MessageNotify("E", module_name, &
               & "The dynamical boundary condition '%c' imposed in the surface is invalid.", c1=DynBCSurface) 
    end select

    select case(DynBCBottom)
       case("Slip")
          DynBC_Bottom = DynBCTYPE_Slip
       case("NoSlip")
          DynBC_Bottom = DynBCTYPE_NoSlip
       case default
          call MessageNotify("E", module_name, &
               & "The dynamical boundary condition '%c' imposed in the bottom is invalid.", c1=DynBCBottom) 
    end select

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, 'KinBC_Surface        = %c', c1 = KinBCSurface  )
    call MessageNotify( 'M', module_name, 'DynBC_Surface        = %c', c1 = DynBCSurface  )
    call MessageNotify( 'M', module_name, 'DynBC_Bottom         = %c', c1 = DynBCBottom  )

  end subroutine read_nmlData

end module BoundCondSet_mod
