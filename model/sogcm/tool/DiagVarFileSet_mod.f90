!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DiagVarFileSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_historyauto

  use VariableSet_mod
  use DiagVarSet_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DiagVarFileSet_Init, DiagVarFileSet_Final
  public :: DiagVarFileSet_OutputVar


  type, public :: gtool_historyauto_info
     real(DP) :: IntValue
     character(TOKEN) :: IntUnit
     character(STRING) :: FilePrefix
     character(STRING) :: Name
  end type gtool_historyauto_info


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DiagVarFileSet_mod' !< Module Name
  

contains

  !>
  !!
  !!
  subroutine DiagVarFileSet_Init(configNmlFileName, diagVar_gthsInfo)

    ! モジュール引用; Use statement
    !

    use Constants_mod, only: RPlanet, PI

    use GridSet_mod, only: &
         & iMax, jMax, kMax, &
         & x_Lon, y_Lat

    use SpmlUtil_mod, only: g_Sig

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo

    ! 局所変数
    ! Local variable
    !
    character(12) :: intUnit

    ! 実行文; Executable statements
    !


    !
    !
    intUnit = diagVar_gthsInfo%intUnit

    call HistoryAutoCreate( &                            ! ヒストリー作成
         & title='OGCM Diagnose Variables Output',             &
         & source='OGCM Output',    &
         & institution='GFD_Dennou Club OGCM project',    &
         & dims=(/'lon','lat','sig','t  '/), dimsizes=(/iMax,jMax,kMax+1,0/),       &
         & longnames=(/'longitude','latitude ','sigma    ', 'time     '/),      &
         & units=(/'degree_east ','degree_north','(1)         ', intUnit /), &
         & namelist_filename=configNmlFileName )    

    call HistoryAutoPutAxis('lon', x_Lon*180/PI)
    call HistoryAutoAddAttr('lon', 'topology', 'circular')
    call HistoryAutoAddAttr('lon', 'modulo', 360.0)
    call HistoryAutoPutAxis('lat', y_Lat*180/PI)
    call HistoryAutoPutAxis('sig', g_Sig)

    call HistoryAutoAddVariable( &
         varname='Chi', dims=(/'lon','lat','sig','t  '/), & 
         longname='velocity potential ', units='m2/s')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_DIV, dims=(/'lon','lat','sig','t  '/), & 
         longname='divergence ', units='s-1')

    call HistoryAutoAddVariable( &
         varname='Psi', dims=(/'lon','lat','sig','t  '/), & 
         longname='stream function', units='m2/s')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_VOR, dims=(/'lon','lat','sig','t  '/), & 
         longname='vorcity', units='s-1')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_TOTPRESS, dims=(/'lon','lat','sig','t  '/), & 
         longname='pressure deviation from RefDens*Grav*z', units='Pa')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_MASSSTREAMFUNC, dims=(/'lat','sig','t  '/), & 
         longname='mass stream function on meriodinal plane', units='kg.m2.s-1')

    call HistoryAutoAllVarFix()

  end subroutine DiagVarFileSet_Init

  !>
  !!
  !!
  subroutine DiagVarFileSet_Final()


    ! 宣言文; Declaration statement
    !

    ! 実行文; Executable statements
    !

    call HistoryAutoClose()

  end subroutine DiagVarFileSet_Final

  !> @brief 
  !!
  !!
  subroutine DiagVarFileSet_OutputVar(CurrentTime, varName, varScalar, var1D, var2D, var3D)

    ! モジュール引用; Use statement
    !
    use GridSet_mod, only: &
         & iMax, jMax, kMax, lMax, &
         & xyz_Lat

    use SpmlUtil_mod

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: CurrentTime
    character(*), intent(in) :: varName
    real(DP), intent(in), optional :: varScalar
    real(DP), intent(in), optional :: var1D(:)
    real(DP), intent(in), optional :: var2D(:,:)
    real(DP), intent(in), optional :: var3D(:,:,:)

    ! 局所変数
    ! Local variables
    !
    

    ! 実行文; Executable statement
    !


    call MessageNotify("M", module_name, "Output data of field '%c' at %d [sec] ..", c1=varName, i=(/ int(CurrentTime) /))

    if(present(varScalar)) &
         call HistoryAutoPut(CurrentTime, varName, varScalar)
    
    if(present(var1D)) &
         call HistoryAutoPut(CurrentTime, varName, var1D)

    if(present(var2D)) &
         call HistoryAutoPut(CurrentTime, varName, var2D)

    if(present(var3D)) &
         call HistoryAutoPut(CurrentTime, varName, var3D)

  end subroutine DiagVarFileSet_OutputVar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

!!$  subroutine read_nmlData( configNmlFileName, &
!!$       & outputFileName )
!!$
!!$    ! モジュール引用; Use statement
!!$    !
!!$
!!$    ! ファイル入出力補助
!!$    ! File I/O support
!!$    !
!!$    use dc_iounit, only: FileOpen
!!$
!!$    ! 種別型パラメタ
!!$    ! Kind type parameter
!!$    !
!!$    use dc_types, only: STDOUT ! 標準出力の装置番号. Unit number of standard output
!!$
!!$    ! 宣言文; Declaration statement
!!$    !
!!$    character(*), intent(in) :: configNmlFileName
!!$    character(STRING), intent(out) :: outputFileName
!!$
!!$    ! 局所変数
!!$    ! Local variables
!!$    !
!!$    integer:: unit_nml        ! NAMELIST ファイルオープン用装置番号. 
!!$                              ! Unit number for NAMELIST file open
!!$
!!$    integer:: iostat_nml      ! NAMELIST 読み込み時の IOSTAT. 
!!$                              ! IOSTAT of NAMELIST read
!!$
!!$    real(DP) :: outputIntrval
!!$
!!$    ! NAMELIST 変数群
!!$    ! NAMELIST group name
!!$    !
!!$    namelist /datafile_nml/ &
!!$         & outputFileName, outputIntrval
!!$
!!$    ! 実行文; Executable statements
!!$
!!$    ! デフォルト値の設定
!!$    ! Default values settings
!!$    !
!!$    outputFileName = "data.nc"
!!$    outputIntrval  = -1d0
!!$
!!$    ! NAMELIST からの入力
!!$    ! Input from NAMELIST
!!$    !
!!$    if ( trim(configNmlFileName) /= '' ) then
!!$       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
!!$       call FileOpen( unit_nml, &             ! (out)
!!$            & configNmlFileName, mode = 'r' ) ! (in)
!!$
!!$       rewind( unit_nml )
!!$       read( unit_nml, &           ! (in)
!!$            & nml = datafile_nml, &  ! (out)
!!$            & iostat = iostat_nml )   ! (out)
!!$       close( unit_nml )
!!$    end if
!!$
!!$
!!$
!!$    ! 印字 ; Print
!!$    !
!!$    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
!!$    call MessageNotify( 'M', module_name, '  outputFileName        = %a', ca=(/ outputFileName /))
!!$    call MessageNotify( 'M', module_name, '  outputIntrval        = %d [sec]', i=(/ int(this%outputIntrval) /))
!!$
!!$  end subroutine read_nmlData

end module DiagVarFileSet_mod

