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
     real(DP) :: origin
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
         & xyz_Lon, xyz_Lat

    use SpmlUtil_mod, only: g_Sig

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName
    type(gtool_historyauto_info), intent(in) :: diagVar_gthsInfo

    ! 局所変数
    ! Local variable
    !
    character(12) :: intUnit
    character(TOKEN) :: dims_YZT(3), dims_XYT(3), dims_XYZT(4)
    
    ! 実行文; Executable statements
    !


    !
    !
    intUnit = diagVar_gthsInfo%intUnit
    dims_YZT(1) = 'lat'; dims_YZT(2) = 'sig'; dims_YZT(3) = 'time'
    dims_XYT(1) = 'lon'; dims_XYT(2) = 'lat'; dims_XYT(3) = 'time'
    dims_XYZT(1) = 'lon'; dims_XYZT(2) = 'lat'; dims_XYZT(3) = 'sig'; dims_XYZT(4) = 'time'
    
    
    call HistoryAutoCreate( &                            ! ヒストリー作成
         & title='OGCM Diagnose Variables Output',             &
         & source='OGCM Output',    &
         & institution='GFD_Dennou Club OGCM project',    &
         & dims=(/'lon ','lat ','sig ','time'/), dimsizes=(/iMax,jMax,kMax+1,0/),       &
         & longnames=(/'longitude','latitude ','sigma    ', 'time     '/),      &
         & units=(/'degree_east ','degree_north','(1)         ', intUnit /), &
         & origin=real(diagVar_gthsInfo%origin), &
         & namelist_filename=configNmlFileName )    

    call HistoryAutoPutAxis('lon', xyz_Lon(:,1,0)*180/PI)
    call HistoryAutoAddAttr('lon', 'topology', 'circular')
    call HistoryAutoAddAttr('lon', 'modulo', 360.0)
    call HistoryAutoPutAxis('lat', xyz_Lat(0,:,0)*180/PI)
    call HistoryAutoPutAxis('sig', g_Sig)

    call HistoryAutoAddVariable( &
         varname='Chi', dims=dims_XYZT, & 
         longname='velocity potential ', units='m2/s')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_DIV, dims=dims_XYZT, &
         longname='divergence ', units='s-1')

    call HistoryAutoAddVariable( &
         varname='Psi', dims=dims_XYZT, &
         longname='stream function', units='m2/s')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_VOR, dims=dims_XYZT, &
         longname='vorcity', units='s-1')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_DENSEDD, dims=dims_XYZT, &
         longname='density deviation from refrence density(RefDens)', units='kg.m-3')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_DENSPOT, dims=dims_XYZT, &
         longname='potential density deviation from refrence density(RefDens). Refrence pressure sets zero.', &
         units='kg.m-3')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_PRESSEDD, dims=dims_XYZT, &
         longname='pressure deviation from static pressure(RefDens*Grav*z)', units='Pa')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_MASSSTREAMFUNC, dims=dims_YZT, &
         longname='mass stream function on meriodinal plane', units='Sv')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_STATICSTABILITY, dims=dims_XYZT, &
         longname='static stability(N^2)', units='s-2')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_PTEMP, dims=dims_XYZT, &
         longname='potential temperature', units='K')

    call HistoryAutoAddVariable( &
         varname=DVARKEY_TEMP, dims=dims_XYZT, & 
         longname='temperature', units='K')



    !
    !
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

    call MessageNotify("M", module_name, "Output data of field '%c' at %f [sec] ..", c1=varName, d=(/ CurrentTime /))

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

