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

  !
    use EOS_Linear_mod, only: &
         & EOSTYPE_LINEAR, EOSTYPENAME_LINEAR

    use EOS_SimpleNonLinear_mod, only: &
         & EOSTYPE_SIMPLENONLINEAR, EOSTYPENAME_SIMPLENONLINEAR

    use EOS_JM95_mod, only: &
         & EOSTYPE_JM95, EOSTYPENAME_JM95
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: GovernEqSet_Init, GovernEqSet_Final
  public :: isPhysicsCompActived

  ! 公開変数
  ! Public variable
  !

  !* The definition of IDs associated with the type of each equation. 
  !

  ! For dynamical equation
  character(*), public, parameter :: GOVERNEQSET_DYN_HYDROBOUSSINESQ_NAME = "HydroBoussinesq"
  integer, public, parameter :: GOVERNEQSET_DYN_HYDROBOUSSINESQ = 1

  ! For equation of state
  integer, public, parameter :: GOVERNEQSET_EOS_LINEAR = EOSTYPE_LINEAR
  integer, public, parameter :: GOVERNEQSET_EOS_SIMPLENONLINEAR = EOSTYPE_SIMPLENONLINEAR
  integer, public, parameter :: GOVERNEQSET_EOS_JM95 = EOSTYPE_JM95

  ! For parametarization of sub-grid scale eddy mixing
  character(*), public, parameter :: GOVERNEQSET_SGSEDDYMIX_Redi_NAME = "Redi"
  integer, public, parameter :: GOVERNEQSET_SGSEDDYMIX_Redi = 1
  character(*), public, parameter :: GOVERNEQSET_SGSEDDYMIX_GM90_NAME = "GM90"
  integer, public, parameter :: GOVERNEQSET_SGSEDDYMIX_GM90 = 2

  ! For convective adjustment
  character(*), public, parameter :: GOVERNEQSET_SGSCONVADJUST_INSTANT_NAME = "Instantaneous"
  integer, public, parameter :: GOVERNEQSET_SGSCONVADJUST_INSTANT = 1
  character(*), public, parameter :: GOVERNEQSET_SGSCONVADJUST_SLOW_NAME = "Slow"
  integer, public, parameter :: GOVERNEQSET_SGSCONVADJUST_SLOW = 2

  !
  character(*), public, parameter :: GOVERNEQSET_TYPENOSPEC_NAME = "UnActivated"
  integer, public, parameter :: GOVERNEQSET_TYPENOSPEC = -1


  !
  integer, public, save :: DynEqType       !< The type of dynamical equations
  integer, public, save :: EOSType         !< The type of equation of state
  integer, public, save :: SGSConvAdjustType  !< The type of convective adjustment
  integer, public, save :: SGSEddyMixType  !< The type of parameterization for sub-grid scale eddy mixing 

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

  logical function isPhysicsCompActived(physicsCompType)
    integer, intent(in) :: physicsCompType

    if(physicsCompType==GOVERNEQSET_TYPENOSPEC) then
       isPhysicsCompActived = .false.
    else
       isPhysicsCompActived = .true.
    end if

  end function isPhysicsCompActived

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

    character(TOKEN) :: DynEqTypeName
    character(TOKEN) :: EOSTypeName
    character(TOKEN) :: SGSConvAdjustTypeName
    character(TOKEN) :: SGSEddyMixTypeName

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /governEq_nml/ &
         & DynEqTypeName, &
         & EOSTypeName, &
         & SGSConvAdjustTypeName, &
         & SGSEddyMixTypeName

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    DynEqTypeName = GOVERNEQSET_DYN_HYDROBOUSSINESQ_NAME
    EOSTypeName = EOSTYPENAME_LINEAR
    SGSConvAdjustTypeName = GOVERNEQSET_SGSCONVADJUST_INSTANT_NAME
    SGSEddyMixTypeName = GOVERNEQSET_TYPENOSPEC_NAME

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

    ! Convert the type name into the corresponding ID. 
    !

    DynEqType = typeName2ID( DynEqTypeName, &
         & (/ GOVERNEQSET_DYN_HYDROBOUSSINESQ_NAME /), &
         & (/ GOVERNEQSET_DYN_HYDROBOUSSINESQ /), &
         & "DynEqType", .false. )

    EOSType = typeName2ID( EOSTypeName, &
         & (/ EOSTYPENAME_LINEAR, EOSTYPENAME_SIMPLENONLINEAR, EOSTYPENAME_JM95 /), &
         & (/ GOVERNEQSET_EOS_LINEAR, GOVERNEQSET_EOS_SIMPLENONLINEAR, GOVERNEQSET_EOS_JM95 /), &
         & "EOSType", .false. )
    
    SGSConvAdjustType = typeName2ID( SGSConvAdjustTypeName, &
         & (/ GOVERNEQSET_SGSCONVADJUST_INSTANT_NAME, GOVERNEQSET_SGSCONVADJUST_SLOW_NAME /), &
         & (/ GOVERNEQSET_SGSCONVADJUST_INSTANT, GOVERNEQSET_SGSCONVADJUST_SLOW /), &
         & "SGSConvAdjustType", .true. )

    SGSEddyMixType = typeName2ID( SGSEddyMixTypeName, &
         & (/ GOVERNEQSET_SGSEDDYMIX_Redi_NAME, GOVERNEQSET_SGSEDDYMIX_GM90_NAME /), &
         & (/ GOVERNEQSET_SGSEDDYMIX_Redi, GOVERNEQSET_SGSEDDYMIX_GM90 /), &
         & "SGSEddyMixType", .true. )


    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '    DynEqType         = %c', c1 = DynEqTypeName ) 
    call MessageNotify( 'M', module_name, '    EOSType           = %c', c1 = EOSTypeName )
    call MessageNotify( 'M', module_name, '    SGSConvAdjustType = %c', c1 = SGSConvAdjustTypeName )
    call MessageNotify( 'M', module_name, '    SGSEddyMixType    = %c', c1 = SGSEddyMixTypeName )

  end subroutine read_nmlData

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! * Subroutines to associate the type name with the coressponding ID. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function typeName2ID(typeName, typeNameList, typeIDList, specTypeInfo, isNoSpecOK) result(typeID)
    character(*), intent(in) :: typeName
    character(*), intent(in) :: typeNameList(:)
    integer, intent(in) :: typeIDList(size(typeNameList))
    character(*), intent(in) :: specTypeInfo
    logical, intent(in) :: isNoSpecOK
    integer :: typeID

    integer :: m

    if(isNoSpecOK) then
       if(trim(typeName)==GOVERNEQSET_TYPENOSPEC_NAME) then
          typeID = GOVERNEQSET_TYPENOSPEC; return
       end if
    end if

    do m=1, size(typeNameList)
       if(trim(typeName)==trim(typeNameList(m))) then
          typeID = typeIDList(m); exit
       end if
       ! The specified type name has not been found in the list. 
       if(m == size(typeNameList)) then
          call MessageNotify('E', module_name, &
               & 'The specified %a ''%a'' is invalid.', ca=(/ trim(specTypeInfo), trim(typeName) /) )
       end if
    end do

  end function typeName2ID

end module GovernEqSet_mod

