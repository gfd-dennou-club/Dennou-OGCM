!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_Admin_GovernEq_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & TOKEN, STRING, DP

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM/Sea-ice
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_Admin_GovernEq_Init, DSIce_Admin_GovernEq_Final
  public :: isPhysicsCompActivated

  ! 公開変数
  ! Public variable
  !

  !* The definition of IDs associated with the type of each equation. 
  !

  ! For dynamical equation

  character(*), public, parameter :: SICEGOVERNEQ_THMDYN_W00_NAME = "Winton2000"
  integer, public, parameter :: SICEGOVERNEQ_THMDYN_W00 = 1
    
  !
  character(*), public, parameter :: SICEGOVERNEQ_TYPENOSPEC_NAME = "UnActivated"
  integer, public, parameter :: SICEGOVERNEQ_TYPENOSPEC = -1


  !
  integer, public, save :: THMDynType       !< The type of thermodynamics sea-ice model
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_Admin_GovernEq_mod' !< Module Name

  
contains

  !>
  !!
  !!
  subroutine DSIce_Admin_GovernEq_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    call read_nmlData(configNmlName)

  end subroutine DSIce_Admin_GovernEq_Init

  !>
  !!
  !!
  subroutine DSIce_Admin_GovernEq_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_Admin_GovernEq_Final

  logical function isPhysicsCompActivated(componentName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: componentName

    ! 実行文; Executable statements
    !
    
    select case(componentName)
    end select

    isPhysicsCompActivated = .false.
    call MessageNotify( 'E', module_name, &
         & 'Unexcepted componentName(=%a) is specified. Check!', ca=(/ componentName /) )
    
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

    character(TOKEN) :: ThermoDynName
    
    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /seaice_governEq_nml/ &
         & ThermoDynName


    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !

    ThermoDynName = SICEGOVERNEQ_THMDYN_W00_NAME

    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                                         ! (in)
            & nml = seaice_governEq_nml, iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    ! - Convert the type name into the corresponding ID ---------
    !

    ! Specify the governing equations used in thermodynamics model

    THMDynType = typeName2ID( ThermoDynName, &
         & (/ SICEGOVERNEQ_THMDYN_W00_NAME /), &
         & (/ SICEGOVERNEQ_THMDYN_W00 /),      &
         & "THMDynType", .false. )
    
    

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, '< Thermodynamics             >')
    call MessageNotify( 'M', module_name, '  - ThermoDynName          = %c', c1 = ThermoDynName ) 

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
       if(trim(typeName)==SICEGOVERNEQ_TYPENOSPEC_NAME) then
          typeID = SICEGOVERNEQ_TYPENOSPEC; return
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

end module DSIce_Admin_GovernEq_mod

