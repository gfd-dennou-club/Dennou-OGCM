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
       & DP, TOKEN

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
  public :: inquire_VBCSpecType

  ! 非公開手続き
  ! Private procedure
  !

  ! Variables to save the boundary condition at the surface
  integer, public, save :: KinBC_Surface
  integer, public, save :: DynBC_Surface
  integer, public, save :: ThermBC_Surface
  integer, public, save :: SaltBC_Surface

  ! Variables to save the boundary condition at the bottom
  integer, public, save :: KinBC_Bottom
  integer, public, save :: DynBC_Bottom
  integer, public, save :: ThermBC_Bottom
  integer, public, save :: SaltBC_Bottom

  ! IDs to identify the type of boundary condition
  !
  integer, public, parameter :: DynBCTYPE_NoSlip = 101
  integer, public, parameter :: DynBCTYPE_Slip = 102
  integer, public, parameter :: DynBCTYPE_SpecStress = 103

  integer, public, parameter :: KinBCTYPE_FreeSurf = 201
  integer, public, parameter :: KinBCTYPE_RigidLid = 202

  integer, public, parameter :: ThermBCTYPE_Adiabat          = 301
  integer, public, parameter :: ThermBCTYPE_PrescFixedFlux   = 302
  integer, public, parameter :: ThermBCTYPE_PrescFlux        = 303
  integer, public, parameter :: ThermBCTYPE_PrescTemp        = 304
  integer, public, parameter :: ThermBCTYPE_TempRelaxed      = 305
  integer, public, parameter :: ThermBCTYPE_PresFlux_Han1984Method  = 306

  integer, public, parameter :: SaltBCTYPE_Adiabat          = 401
  integer, public, parameter :: SaltBCTYPE_PrescFixedFlux   = 402
  integer, public, parameter :: SaltBCTYPE_PrescFlux        = 403
  integer, public, parameter :: SaltBCTYPE_PrescSalt        = 404
  integer, public, parameter :: SaltBCTYPE_SaltRelaxed      = 405
  integer, public, parameter :: SaltBCTYPE_PresFlux_Han1984Method  = 406

  ! Labels to identify the type of boundary condition
  !
  character(*), public, parameter :: DynBCTYPELBL_NoSlip = 'NoSlip'
  character(*), public, parameter :: DynBCTYPELBL_Slip = 'Slip'
  character(*), public, parameter :: DynBCTYPELBL_SpecStress = 'SpecStress'

  character(*), public, parameter :: KinBCTYPELBL_FreeSurf = 'Free'
  character(*), public, parameter :: KinBCTYPELBL_RigidLid = 'Rigid'

  character(*), public, parameter :: ThermBCTYPELBL_Adiabat        = 'Adiabat'
  character(*), public, parameter :: ThermBCTYPELBL_PrescFlux      = 'PrescFlux'
  character(*), public, parameter :: ThermBCTYPELBL_PrescFixedFlux = 'PrescFixedFlux'
  character(*), public, parameter :: ThermBCTYPELBL_PrescTemp      = 'PrescTemp'
  character(*), public, parameter :: ThermBCTYPELBL_TempRelaxed    = 'TempRelaxed'
  character(*), public, parameter :: ThermBCTYPELBL_PrescFlux_Han1984Method = 'PrescFlux_Han1984'

  character(*), public, parameter :: SaltBCTYPELBL_Adiabat         = 'Adiabat'
  character(*), public, parameter :: SaltBCTYPELBL_PrescFlux       = 'PrescFlux'
  character(*), public, parameter :: SaltBCTYPELBL_PrescFixedFlux  = 'PrescFixedFlux'  
  character(*), public, parameter :: SaltBCTYPELBL_PrescSalt       = 'PrescSalt'
  character(*), public, parameter :: SaltBCTYPELBL_SaltRelaxed     = 'SaltRelaxed'
  character(*), public, parameter :: SaltBCTYPELBL_PrescFlux_Han1984Method = 'PrescFlux_Han1984'

  !
  !
  real(DP), public :: SurfTempRelaxedTime
  real(DP), public :: SurfSaltRelaxedTime

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


  !> @brief Inquire the type of boundary condition. 
  !! Inquire whether the type of vertical boundary condition is first-type or second-type. 
  !!
  !! @return 'D' or 'N'
  !!
  function inquire_VBCSpecType(VBCTypeID) result(VBCSpecType)
    
    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: VBCTypeID
    character :: VBCSpecType
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    select case(VBCTypeID)
       case(KinBCTYPE_FreeSurf)      !-- kinetic boundary condition ------------ 
       case(KinBCTYPE_RigidLid)
       case(DynBCTYPE_Slip)
          VBCSpecType = 'N'
       case(DynBCTYPE_NoSlip)        !-- dynamical boundary condition ------------ 
          VBCSpecType = 'D'
       case(DynBCTYPE_SpecStress)
          VBCSpecType = 'N'
       case(ThermBCTYPE_PrescTemp, ThermBCTYPE_TempRelaxed) !-- thermal boundary condition ------------ 
          VBCSpecType = 'D'
       case(ThermBCTYPE_PrescFixedFlux, ThermBCTYPE_PrescFlux, ThermBCTYPE_Adiabat, &
          & ThermBCTYPE_PresFlux_Han1984Method)
          VBCSpecType = 'N'
       case(SaltBCTYPE_PrescSalt, SaltBCTYPE_SaltRelaxed)     !-- Salinity boundary condition ------------ 
          VBCSpecType = 'D'
       case(SaltBCTYPE_PrescFixedFlux, SaltBCTYPE_PrescFlux, SaltBCTYPE_Adiabat, &
            & SaltBCTYPE_PresFlux_Han1984Method)
          VBCSpecType = 'N'
       case Default
          call MessageNotify("E", module_name, &
               & "The ID of inquired boundary condition '%d' is not registered.", i=(/VBCTypeID/) ) 
    end select

  end function inquire_VBCSpecType


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

    character(TOKEN) :: &
         & KinBCSurface, KinBCBottom, &
         & DynBCSurface, DynBCBottom, &
         & ThermBCSurface, ThermBCBottom, &
         & SaltBCSurface, SaltBCBottom

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /boundaryCondition_nml/ &
         & KinBCSurface, KinBCBottom, &
         & DynBCSurface, DynBCBottom, &
         & ThermBCSurface, ThermBCBottom, &
         & SaltBCSurface, SaltBCBottom, &
         & SurfTempRelaxedTime, SurfSaltRelaxedTime


    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !


    ! Set default boundary conditions
    !

    KinBCSurface = KinBCTYPELBL_RigidLid
    DynBCSurface = DynBCTYPELBL_NoSlip
    ThermBCSurface = ThermBCTYPELBL_Adiabat
    SaltBCSurface = SaltBCTYPELBL_Adiabat

    KinBCBottom = KinBCTYPELBL_RigidLid
    DynBCBottom = DynBCTYPELBL_NoSlip
    ThermBCBottom = ThermBCTYPELBL_Adiabat
    SaltBCBottom = SaltBCTYPELBL_Adiabat

    SurfTempRelaxedTime = -1d0
    SurfSaltRelaxedTime = -1d0
    
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

    ! Set the IDs of boundary condition on upper surface.
    !
    call label2ID_verticalBCType( &
         & KinBCSurface, DynBCSurface, ThermBCSurface, SaltBCSurface,     &  ! (in)
         & KinBC_Surface, DynBC_Surface, ThermBC_Surface, SaltBC_Surface, &  ! (out)
         & .true. )

    ! Set  the IDs of boundary condition on bottom surface.
    !
    call label2ID_verticalBCType( &
         & KinBCBottom, DynBCBottom, ThermBCBottom, SaltBCBottom,     & ! (in)
         & KinBC_Bottom, DynBC_Bottom, ThermBC_Bottom, SaltBC_Bottom, & ! (out)
         & .true. )


    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, 'KinBC_Surface        = %c', c1 = KinBCSurface  )
    call MessageNotify( 'M', module_name, 'DynBC_Surface        = %c', c1 = DynBCSurface  )
    call MessageNotify( 'M', module_name, 'ThermBC_Surface      = %c', c1 = ThermBCSurface  )
    call MessageNotify( 'M', module_name, 'SaltBC_Surface      = %c', c1 = SaltBCSurface  )
    call MessageNotify( 'M', module_name, 'KinBC_Bottom         = %c', c1 = KinBCBottom  )
    call MessageNotify( 'M', module_name, 'DynBC_Bottom         = %c', c1 = DynBCBottom  )
    call MessageNotify( 'M', module_name, 'ThermBC_Bottom       = %c', c1 = ThermBCBottom  )
    call MessageNotify( 'M', module_name, 'SaltBC_Bottom       = %c', c1 = SaltBCBottom  )

    if( ThermBC_Surface==ThermBCTYPE_TempRelaxed ) then
       if ( SurfTempRelaxedTime >= 0d0 ) then
          call MessageNotify('M', module_name, 'SurfTempRelaxedTime  = %f [sec]', d=(/ SurfTempRelaxedTime /))
       else
          call MessageNotify('E', module_name, ' &
               & `ThermBC_Surface=ThermBCTYPE_TempRelaxed`, but `SurfTempRelaxedTime` is not specified.' &
               & // 'Set the value of this parameter.' )
       end if
    end if

    if( SaltBC_Surface==SaltBCTYPE_SaltRelaxed ) then
       if ( SurfSaltRelaxedTime >= 0d0 ) then
          call MessageNotify('M', module_name, 'SurfSaltRelaxedTime  = %f [sec]', d=(/ SurfSaltRelaxedTime /))
       else
          call MessageNotify('E', module_name, ' &
               & `SaltBC_Surface=SaltCTYPE_SaltRelaxed`, but `SurfSaltRelaxedTime` is not specified.' &
               & // 'Set the value of this parameter.' )
       end if
    end if

  end subroutine read_nmlData

  !> @brief 
  !!
  !!
  subroutine label2ID_verticalBCType( &
       & KinBCLBL, DynBCLBL, ThermBCLBL, SaltBCLBL, &   ! (in)
       & KinBCID, DynBCID, ThermBCID, SaltBCID,     &   ! (out)
       & isExceptionCatch )                             ! (in)
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: KinBCLBL, DynBCLBL, ThermBCLBL, SaltBCLBL
    integer, intent(out) :: KinBCID, DynBCID, ThermBCID, SaltBCID
    logical, intent(in) :: isExceptionCatch

    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    select case(KinBCLBL)
       case(KinBCTYPELBL_FreeSurf)
          KinBCID = KinBCTYPE_FreeSurf
       case(KinBCTYPELBL_RigidLid)
          KinBCID = KinBCTYPE_RigidLid
       case default
          call MessageNotify("E", module_name, &
               & "The Kinetic boundary condition '%c' is not available.", c1=trim(KinBCLBL) ) 
    end select

    select case(DynBCLBL)
       case(DynBCTYPELBL_Slip)
          DynBCID = DynBCTYPE_Slip
       case(DynBCTYPELBL_NoSlip)
          DynBCID = DynBCTYPE_NoSlip
       case(DynBCTYPELBL_SpecStress)
          DynBCID = DynBCTYPE_SpecStress
       case default
          call MessageNotify("E", module_name, &
               & "The dynamical boundary condition '%c' is not available.", c1=trim(DynBCLBL) )
    end select

    select case (ThermBCLBL)
       case(ThermBCTYPELBL_Adiabat)
          ThermBCID = ThermBCTYPE_Adiabat
       case(ThermBCTYPELBL_PrescTemp)
          ThermBCID = ThermBCTYPE_PrescTemp
       case(ThermBCTYPELBL_PrescFixedFlux)
          ThermBCID = ThermBCTYPE_PrescFixedFlux
       case(ThermBCTYPELBL_PrescFlux)
          ThermBCID = ThermBCTYPE_PrescFlux
       case(ThermBCTYPELBL_PrescFlux_Han1984Method)
          ThermBCID = ThermBCTYPE_PresFlux_Han1984Method
       case(ThermBCTYPELBL_TempRelaxed)
          ThermBCID = ThermBCTYPE_TempRelaxed
       case default
          call MessageNotify("E", module_name, &
               & "The thermal boundary condition '%c' is not available.", c1=trim(ThermBCLBL) )
    end select

    select case (SaltBCLBL)
       case(SaltBCTYPELBL_Adiabat)
          SaltBCID = SaltBCTYPE_Adiabat
       case(SaltBCTYPELBL_PrescSalt)
          SaltBCID = SaltBCTYPE_PrescSalt
       case(SaltBCTYPELBL_PrescFixedFlux)
          SaltBCID = SaltBCTYPE_PrescFixedFlux
       case(SaltBCTYPELBL_PrescFlux)
          SaltBCID = SaltBCTYPE_PrescFlux
       case(SaltBCTYPELBL_PrescFlux_Han1984Method)
          SaltBCID = SaltBCTYPE_PresFlux_Han1984Method
       case(SaltBCTYPELBL_SaltRelaxed)
          SaltBCID = SaltBCTYPE_SaltRelaxed
       case default
          call MessageNotify("E", module_name, &
               & "The salinity boundary condition '%c' is not available.", c1=trim(SaltBCLBL) )
    end select

  end subroutine label2ID_verticalBCType

end module BoundCondSet_mod
