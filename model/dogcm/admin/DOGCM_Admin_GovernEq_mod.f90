!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_Admin_GovernEq_mod 

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
  public :: DOGCM_Admin_GovernEq_Init, DOGCM_Admin_GovernEq_Final
  public :: isPhysicsCompActivated

  ! 公開変数
  ! Public variable
  !

  !* The definition of IDs associated with the type of each equation. 
  !

  ! For dynamical equation

  character(*), public, parameter :: OCNGOVERNEQ_DYN_HYDROBOUSSINESQ_NAME = "HydroBoussinesq"
  integer, public, parameter :: OCNGOVERNEQ_DYN_HYDROBOUSSINESQ = 1

  character(*), public, parameter :: OCNGOVERNEQ_DYN_NONDYN_MIXEDLYR_NAME = "NonDynMixedLyr"
  integer, public, parameter :: OCNGOVERNEQ_DYN_NONDYN_MIXEDLYR = 2

  character(*), public, parameter :: OCNGOVERNEQ_DYN_NONDYN_NAME = "NonDyn"
  integer, public, parameter :: OCNGOVERNEQ_DYN_NONDYN = 3
  
  ! For equation of state
  integer, public, parameter :: OCNGOVERNEQ_EOS_LINEAR = EOSTYPE_LINEAR
  integer, public, parameter :: OCNGOVERNEQ_EOS_SIMPLENONLINEAR = EOSTYPE_SIMPLENONLINEAR
  integer, public, parameter :: OCNGOVERNEQ_EOS_JM95 = EOSTYPE_JM95


  !
  character(*), public, parameter :: OCNGOVERNEQ_LPHYS_MIXMOM_NAME  = 'LMixMOM'
  character(*), public, parameter :: OCNGOVERNEQ_LPHYS_MIXTRC_NAME  = 'LMixTRC'
  character(*), public, parameter :: OCNGOVERNEQ_LPHYS_REDIGM_NAME  = 'RediGM'

  character(*), public, parameter :: OCNGOVERNEQ_VPHYS_MIXMOM_NAME  = 'VMixMOM'
  character(*), public, parameter :: OCNGOVERNEQ_VPHYS_MIXTRC_NAME  = 'VMixTRC'
  character(*), public, parameter :: OCNGOVERNEQ_VPHYS_CONVEC_NAME  = 'Convect'

  !
  character(*), public, parameter :: OCNGOVERNEQ_TYPENOSPEC_NAME = "UnActivated"
  integer, public, parameter :: OCNGOVERNEQ_TYPENOSPEC = -1


  !
  integer, public, save :: DynEqType       !< The type of dynamical equations
  integer, public, save :: EOSType         !< The type of equation of state
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_Admin_GovernEq_mod' !< Module Name

  logical :: LPHYS_LMixMOM_Flag
  logical :: LPHYS_LMixTRC_Flag
  logical :: LPHYS_RediGM_Flag

  logical :: VPHYS_VMixMOM_Flag
  logical :: VPHYS_VMixTRC_Flag
  logical :: VPHYS_Convect_Flag
  
contains

  !>
  !!
  !!
  subroutine DOGCM_Admin_GovernEq_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    call read_nmlData(configNmlName)

  end subroutine DOGCM_Admin_GovernEq_Init

  !>
  !!
  !!
  subroutine DOGCM_Admin_GovernEq_Final()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_Admin_GovernEq_Final

  logical function isPhysicsCompActivated(componentName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: componentName

    ! 実行文; Executable statements
    !
    
    select case(componentName)
       !-------- LPHYS --------------------       
    case (OCNGOVERNEQ_LPHYS_MIXMOM_NAME)
       isPhysicsCompActivated = LPHYS_LMixMOM_Flag
    case (OCNGOVERNEQ_LPHYS_MIXTRC_NAME)
       isPhysicsCompActivated = LPHYS_LMixTRC_Flag
    case (OCNGOVERNEQ_LPHYS_REDIGM_NAME)
       isPhysicsCompActivated = LPHYS_RediGM_Flag
       !-------- VPHYS -------------------
    case (OCNGOVERNEQ_VPHYS_MIXMOM_NAME) 
       isPhysicsCompActivated = VPHYS_VMixMOM_Flag
    case (OCNGOVERNEQ_VPHYS_MIXTRC_NAME)
       isPhysicsCompActivated = VPHYS_VMixTRC_Flag
    case (OCNGOVERNEQ_VPHYS_CONVEC_NAME)
       isPhysicsCompActivated = VPHYS_Convect_Flag
    end select
    
  end function isPhysicsCompActivated
  
!-----------------------------------------------------------

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

    !
    use dc_string, only: Split, Replace, StrInclude

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
    character(STRING) :: LPhysNames
    character(STRING) :: VPhysNames

    character(TOKEN), pointer :: LPhysNameList(:)
    character(TOKEN), pointer :: VPhysNameList(:)
    
    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /GovernEq_nml/ &
         & DynEqTypeName,   &
         & EOSTypeName,     &
         & LPhysNames,      &
         & VPhysNames


    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    DynEqTypeName = OCNGOVERNEQ_DYN_HYDROBOUSSINESQ_NAME
    EOSTypeName = EOSTYPENAME_LINEAR

    LPhysNames = 'LMixMom, LMixTRC'
    VPhysNames = 'VMixMOM, VMixTRC'
    
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

    ! - Convert the type name into the corresponding ID ---------
    !

    ! Specify the governing equations used in dynamical core
    
    DynEqType = typeName2ID( DynEqTypeName, &
         & (/ OCNGOVERNEQ_DYN_HYDROBOUSSINESQ_NAME,     &
         &    OCNGOVERNEQ_DYN_NONDYN_MIXEDLYR_NAME,     &
         &    OCNGOVERNEQ_DYN_NONDYN_NAME           /), &
         & (/ OCNGOVERNEQ_DYN_HYDROBOUSSINESQ,          &
         &    OCNGOVERNEQ_DYN_NONDYN_MIXEDLYR,          &
         &    OCNGOVERNEQ_DYN_NONDYN                /), &
         & "DynEqType", .false. )

    EOSType = typeName2ID( EOSTypeName, &
         & (/ EOSTYPENAME_LINEAR, EOSTYPENAME_SIMPLENONLINEAR, EOSTYPENAME_JM95 /), &
         & (/ OCNGOVERNEQ_EOS_LINEAR, OCNGOVERNEQ_EOS_SIMPLENONLINEAR, OCNGOVERNEQ_EOS_JM95 /), &
         & "EOSType", .false. )
    
    ! Set the flag to represent whether each package is activated.
    !

    ! Lateral ocean physics
    
    LPhysNames = Replace( LPhysNames, " ", "")
    LPhysNameList => null()
    call Split(trim(LPhysNames), LPhysNameList, ",")

    if ( StrInclude( LPhysNameList, OCNGOVERNEQ_LPHYS_MIXMOM_NAME) ) then
       LPHYS_LMixMOM_Flag = .true.
    end if
    if ( StrInclude( LPhysNameList, OCNGOVERNEQ_LPHYS_MIXTRC_NAME) ) then
       LPHYS_LMixTRC_Flag = .true.
    end if
    if ( StrInclude( LPhysNameList, OCNGOVERNEQ_LPHYS_REDIGM_NAME) ) then
       LPHYS_RediGM_Flag = .true.
    end if

    deallocate( LPhysNameList )
    
    ! Vertical ocean physics
    
    VPhysNames = Replace( VPhysNames, " ", "")
    VPhysNameList => null()
    call Split(trim(VPhysNames), VPhysNameList, ",")

    if ( StrInclude( VPhysNameList, OCNGOVERNEQ_VPHYS_MIXMOM_NAME) ) then
       VPHYS_VMixMOM_Flag = .true.
    end if
    if ( StrInclude( VPhysNameList, OCNGOVERNEQ_VPHYS_MIXTRC_NAME) ) then
       VPHYS_VMixTRC_Flag = .true.
    end if
    if ( StrInclude( VPhysNameList, OCNGOVERNEQ_VPHYS_CONVEC_NAME) ) then
       VPHYS_Convect_Flag = .true.
    end if
    
    deallocate( VPhysNameList )

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '< Ocean Dynamics             >')
    call MessageNotify( 'M', module_name, '  - DynEqType            = %c', c1 = DynEqTypeName ) 
    call MessageNotify( 'M', module_name, '  - EOSType              = %c', c1 = EOSTypeName )
    call MessageNotify( 'M', module_name, '< Ocean Lateral Physics      >')
    call MessageNotify( 'M', module_name, '  - LMixMOM_Flag         = %b', l = (/ LPHYS_LMixMOM_Flag /))
    call MessageNotify( 'M', module_name, '  - LMixTRC_Flag         = %b', l = (/ LPHYS_LMixTRC_Flag /))
    call MessageNotify( 'M', module_name, '  - RediGM_Flag         = %b', l = (/ LPHYS_RediGM_Flag /))
    call MessageNotify( 'M', module_name, '< Ocean Vertical Physics     >')
    call MessageNotify( 'M', module_name, '  - VMixMOM_Flag         = %b', l = (/ VPHYS_VMixMOM_Flag /))
    call MessageNotify( 'M', module_name, '  - VMixTRC_Flag         = %b', l = (/ VPHYS_VMixTRC_Flag /))
    call MessageNotify( 'M', module_name, '  - Convect_Flag         = %b', l = (/ VPHYS_Convect_Flag /))


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
       if(trim(typeName)==OCNGOVERNEQ_TYPENOSPEC_NAME) then
          typeID = OCNGOVERNEQ_TYPENOSPEC; return
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

end module DOGCM_Admin_GovernEq_mod

