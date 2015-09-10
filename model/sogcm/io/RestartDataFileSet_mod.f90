!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module RestartDataFileSet_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_history, only: &
       & gt_history

  use VariableSet_mod, only: &
       & xyz_UN, xyz_VN, xyz_PTempEddN, xyz_SaltN, xy_SurfHeightN, &
       & VARSET_KEY_U, VARSET_KEY_V, VARSET_KEY_PTEMPEDD, VARSET_KEY_SALT, VARSET_KEY_SURFHEIGHT, &
       & xyz_UB, xyz_VB, xyz_PTempEddB, xyz_SaltB, xy_SurfHeightB, &
       & VARSET_KEY_UB, VARSET_KEY_VB, VARSET_KEY_PTEMPEDDB, VARSET_KEY_SALTB, VARSET_KEY_SURFHEIGHTB, &
       & z_PTempBasic, VARSET_KEY_PTEMPBASIC

  use VarSetSeaice_mod, only: &
       & xy_SIceConN, xy_IceThickN, xy_SnowThickN, xyz_SIceTempN, xy_SIceSurfTempN, &
       & VARSET_KEY_SICECON, VARSET_KEY_ICETHICK, VARSET_KEY_SNOWTHICK, VARSET_KEY_SICETEMP, VARSET_KEY_SICESURFTEMP, &
       & xy_SIceConB, xy_IceThickB, xy_SnowThickB, xyz_SIceTempB, xy_SIceSurfTempB, &
       & VARSET_KEY_SICECONB, VARSET_KEY_ICETHICKB, VARSET_KEY_SNOWTHICKB, VARSET_KEY_SICETEMPB, VARSET_KEY_SICESURFTEMPB
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: RestartDataFileSet_Init, RestartDataFileSet_Final
  public :: RestartDataFileSet_Input, RestartDataFileSet_Output

  character(STRING), save, public :: InputFileName
  character(STRING), save, public :: OutputFileName
  real(DP), save, public :: IntValue
  character(TOKEN), save, public :: IntUnit
  real(DP), save, public :: IntTime
  logical, save, public :: RestartFlag
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'RestartDataFileSet_mod' !< Module Name
  type(gt_history), save :: gthst_rst
  
contains

  !>
  !!
  !!
  subroutine RestartDataFileSet_Init(configNmlFileName)

    use gtool_history

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 実行文; Executable statements
    !

    !
    RestartFlag = .false.
    
    !
    call read_nmlData(configNmlFileName)
 
    !
    call createRestartFile()

  end subroutine RestartDataFileSet_Init

  !>
  !!
  !!
  subroutine RestartDataFileSet_Final()

    use gtool_history, only: HistoryClose

    ! 実行文; Executable statements
    !

    call HistoryClose(history=gthst_rst)

  end subroutine RestartDataFileSet_Final

  !> @brief 
  !!
  !!
  subroutine RestartDataFileSet_Output()
    
    ! モジュール引用; Use statement
    !
    use gtool_history, only: &
         & HistorySetTime, HistoryPut

    use dc_calendar, only: &
         & DCCalDateInquire

    use TemporalIntegSet_mod, only: &
         & CurrentTime, InitDate

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: dateStr
    
    ! 実行文; Executable statement
    !
    if( mod(CurrentTime, IntTime) /= 0 ) return 

    call MessageNotify('M', module_name, 'Output data to restart..')


    ! Output data of some variables for restarting a simulation.
    
    call HistorySetTime(timed=CurrentTime, history=gthst_rst)
    call DCCalDateInquire(dateStr, &
         & elapse_sec=CurrentTime, date=InitDate)
    
    call HistoryPut('datetime', dateStr, history=gthst_rst)

    ! Output prognostic variables in OGCM
    !
    call HistoryPut(VARSET_KEY_U, xyz_UN, history=gthst_rst)
    call HistoryPut(VARSET_KEY_V, xyz_VN, history=gthst_rst)
    call HistoryPut(VARSET_KEY_PTEMPEDD, xyz_PTempEddN, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SALT, xyz_SaltN, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SURFHEIGHT, xy_SurfHeightN, history=gthst_rst)
    
    call HistoryPut(VARSET_KEY_UB, xyz_UB, history=gthst_rst)
    call HistoryPut(VARSET_KEY_VB, xyz_VB, history=gthst_rst)
    call HistoryPut(VARSET_KEY_PTEMPEDDB, xyz_PTempEddB, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SALTB, xyz_SaltB, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SURFHEIGHTB, xy_SurfHeightB, history=gthst_rst)

    call HistoryPut(VARSET_KEY_PTEMPBASIC, z_PTempBasic, history=gthst_rst)
    
    ! Output prognostic variables in seaice model
    !
    call HistoryPut(VARSET_KEY_SICECON, xy_SIceConN, history=gthst_rst)
    call HistoryPut(VARSET_KEY_ICETHICK, xy_IceThickN, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SNOWTHICK, xy_SnowThickN, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SICETEMP, xyz_SIceTempN, history=gthst_rst)

    call HistoryPut(VARSET_KEY_SICECONB, xy_SIceConB, history=gthst_rst)
    call HistoryPut(VARSET_KEY_ICETHICKB, xy_IceThickB, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SNOWTHICKB, xy_SnowThickB, history=gthst_rst)
    call HistoryPut(VARSET_KEY_SICETEMPB, xyz_SIceTempB, history=gthst_rst)
    
    
  end subroutine RestartDataFileSet_Output

  !> @brief 
  !!
  !!
  subroutine RestartDataFileSet_Input(isRestartDataInput)
    
    
    ! モジュール引用; Use statement
    !
    use dc_string, only: &
         & toChar

    use gtool_history, only: &
         & HistoryGetAttr, HistoryGet
    
    use TemporalIntegSet_mod, only: &
         & RestartTime

    use GridSet_mod, only: iMax, jMax, kMax
    
    ! 宣言文; Declaration statement
    !
    logical, intent(out), optional :: isRestartDataInput
    
    ! 局所変数
    ! Local variables
    !
    
    character(STRING) :: timeRange
    character(TOKEN) :: dummyStr
    logical :: getError
    
    ! 実行文; Executable statement
    !
    
    if( InputFileName == '' ) then
       RestartFlag = .false.
       return
    end if

    ! Input data of prognostic variables output at the previous simulation.
    
    RestartFlag = .true.

    call MessageNotify( 'M', module_name, &
         & "Initial/Restart data is input from '%c'..", c1=InputFileName )

    call HistoryGetAttr(InputFileName, 'lon', 'units', &
         & dummyStr, err=getError)
    if( getError ) then
       call MessageNotify('E', module_name, &
            & "Initial/Restart data file '%c' is not found.", c1=trim(InputFileName) ) 
    end if
    
    timeRange = 'time=' // toChar(RestartTime)

    ! Input prognostic variables in OGCM.
    !
    call input_XYZVar(VARSET_KEY_U, xyz_UN)
    call input_XYZVar(VARSET_KEY_V, xyz_VN)
    call input_XYZVar(VARSET_KEY_PTEMPEDD, xyz_PTempEddN)
    call input_XYZVar(VARSET_KEY_SALT, xyz_SaltN)
    
    call input_XYVar(VARSET_KEY_SURFHEIGHT, xy_SurfHeightN)
    
    call input_XYZVar(VARSET_KEY_UB, xyz_UB)
    call input_XYZVar(VARSET_KEY_VB, xyz_VB)
    call input_XYZVar(VARSET_KEY_PTEMPEDDB, xyz_PTempEddB)
    call input_XYZVar(VARSET_KEY_SALTB, xyz_SaltB)
    call input_XYVar(VARSET_KEY_SURFHEIGHTB, xy_SurfHeightB)

    call input_ZVar(VARSET_KEY_PTEMPBASIC, z_PTempBasic)
    
    ! Input prognostic variables in seaice model.
    !
    call input_XYVar(VARSET_KEY_SICECON, xy_SIceConN)
    call input_XYVar(VARSET_KEY_ICETHICK, xy_IceThickN)
    call input_XYVar(VARSET_KEY_SNOWTHICK, xy_SnowThickN)
    call HistoryGet( InputFileName, VARSET_KEY_SICETEMP, range=timeRange, &
         & array=xyz_SIceTempN )

    call input_XYVar(VARSET_KEY_SICECONB, xy_SIceConB)
    call input_XYVar(VARSET_KEY_ICETHICKB, xy_IceThickB)
    call input_XYVar(VARSET_KEY_SNOWTHICKB, xy_SnowThickB)
    call HistoryGet( InputFileName, VARSET_KEY_SICETEMPB, range=timeRange, &
         & array=xyz_SIceTempB )
    
  contains
    subroutine input_XYZVar( varName, xyz )
      
      character(*), intent(in) :: varName
      real(DP), intent(out) :: xyz(0:iMax-1,jMax,0:kMax)
      
      call HistoryGet( InputFileName, varName, range=timeRange, &
           & array=xyz )
      
    end subroutine input_XYZVar

    subroutine input_XYVar( varName, xy )
      
      character(*), intent(in) :: varName
      real(DP), intent(out) :: xy(0:iMax-1,jMax)
      
      call HistoryGet( InputFileName, varName, range=timeRange, &
           & array=xy )
      
    end subroutine input_XYVar

    subroutine input_ZVar( varName, z )
      
      character(*), intent(in) :: varName
      real(DP), intent(out) :: z(0:kMax)
      
      call HistoryGet( InputFileName, varName, range=timeRange, &
           & array=z )
      
    end subroutine input_ZVar

  end subroutine RestartDataFileSet_Input


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

    use dc_calendar, only: DCCalConvertByUnit

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
    namelist /restartFile_nml/ &
         & InputFileName, OutputFileName, IntValue, IntUnit

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
            & nml = restartFile_nml, &  ! (out)
            & iostat = iostat_nml )   ! (out)
       close( unit_nml )
    end if

    IntTime = DCCalConvertByUnit(IntValue, IntUnit, 'sec')

    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, ' OutputFileName = %c', c1=OutputFileName )
    call MessageNotify( 'M', module_name, ' InputFileName  = %c', c1=InputFileName )
    call MessageNotify( 'M', module_name, ' IntValue       = %f', d=(/ IntValue /) )
    call MessageNotify( 'M', module_name, ' IntUnit        = %c', c1=IntUnit  )
    call MessageNotify( 'M', module_name, ' IntTime        = %f', d=(/ IntTime  /) )

  end subroutine read_nmlData

  !> @brief 
  !!
  !!
  subroutine createRestartFile()
    
    !
    !
    use gtool_history, only: &
         & HistoryCreate, &
         & HistoryAddAttr, HistoryAddVariable, &
         & HistoryPut

    use TemporalIntegSet_mod, only: &
         & RestartTime

    use Constants_mod, only: &
         & PI

    use GridSet_mod, only: &
         & iMax, jMax, kMax, &
         & xyz_Lon, xyz_Lat

    use SpmlUtil_mod, only: &
         & g_Sig

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: Z(1), XYT(3), XYZT(4), XYZ2T(4)
    
    ! 実行文; Executable statement
    !
    
    !

    call HistoryCreate( &
         & file=OutputFileName, title=' ', source=' ',  institution=' ',      &
         & dims=(/ 'lon    ', 'lat    ', 'sig    ', 'sig2   ', 'timestr', 'time   ' /),  &
         & dimsizes=(/ iMax, jMax, kMax+1, 2, TOKEN, 0 /),                   &
         & longnames =(/ 'longitude                            ',   &
         &               'latitude                             ',   &
         &               'sigma at collocation points          ',   &
         &               'sigma at collocation points(sea-ice) ',   &         
         &               'number of characters for datetime    ',   &
         &               'time                                 '/), &
         & units=(/'degree_east ', 'degree_north', '1           ',          &
         &         '1           ', '1           ', 'sec         ' /),       &
         & xtypes=(/'double','double','double','int   ', 'double'/),      &
         & origind=RestartTime, intervald=IntTime, history=gthst_rst )
 
    ! 年月日時分秒形式の日時情報変数の設定
    ! Set a date information variable with year-month-day hour:minute:second format
    !   
    call HistoryAddVariable( &
         & varname='datetime', dims=(/'timestr','time   ' /), &
         & longname='time represented as strings', units='1', xtype='char', &
         & history=gthst_rst )

    ! 座標データの設定
    ! Axes data settings
    !
    call HistoryAddAttr( &
      & 'lon', attrname = 'standard_name', &     ! (in)
      & value = 'longitude', &                   ! (in)
      & history = gthst_rst )                    ! (inout)
    call HistoryAddAttr( &
      & 'lat', attrname = 'standard_name', &     ! (in)
      & value = 'latitude', &                    ! (in)
      & history = gthst_rst )                    ! (inout)
    call HistoryAddAttr( &
      & 'sig', attrname = 'standard_name', &     ! (in)
      & value = 'ocean_sigma_coordinate', & ! (in)
      & history = gthst_rst )                    ! (inout)
    call HistoryAddAttr( &
      & 'time', attrname = 'standard_name', &    ! (in)
      & value = 'time', &                        ! (in)
      & history = gthst_rst )                    ! (inout)
    call HistoryAddAttr( &
      & 'sig', attrname = 'positive', &          ! (in)
      & value = 'down', &                        ! (in)
      & history = gthst_rst )                    ! (inout)

    call HistoryPut( &
         & 'lon', xyz_Lon(:,1,1)/PI*180d0, & ! (in)
         & history = gthst_rst )             ! (inout)
    call HistoryPut( &
         & 'lat', xyz_Lat(0,:,1)/PI*180d0, & ! (in)
         & history = gthst_rst )             ! (inout)
    call HistoryPut( &
         & 'sig', g_Sig, &                 ! (in)
         & history = gthst_rst )           ! (inout)
    call HistoryPut( &
         & 'sig2', (/ -0.25d0, -0.75d0 /), & ! (in)
         & history = gthst_rst )             ! (inout)
    ! 
    !
    Z(1) = 'sig'
    XYT(1)  = 'lon'; XYT(2)  = 'lat'; XYT(3) = 'time' 
    XYZT(1) = 'lon'; XYZT(2) = 'lat'; XYZT(3) = 'sig'; XYZT(4) = 'time'
    XYZ2T(1) = 'lon'; XYZ2T(2) = 'lat'; XYZ2T(3) = 'sig2'; XYZ2T(4) = 'time'

    ! Regist variables in  OGCM
    call registVar(VARSET_KEY_U, XYZT, 'zonal velociy', 'm s-1')
    call registVar(VARSET_KEY_V, XYZT, 'meridional velociy', 'm s-1')
    call registVar(VARSET_KEY_PTEMPBASIC, Z, 'basic state of potential temperature', 'K')    
    call registVar(VARSET_KEY_PTEMPEDD, XYZT, 'eddy component of potential temerature', 'K')
    call registVar(VARSET_KEY_SALT, XYZT, 'salinity', 'psu')
    call registVar(VARSET_KEY_SURFHEIGHT, XYT, 'surface height', 'm')
    
    ! Regist variables in sea-ice model
    call registVar(VARSET_KEY_SICECON, XYT, 'Concentration of sea-ice', '(1)')
    call registVar(VARSET_KEY_SNOWTHICK, XYT, 'Snow thickness', 'm')
    call registVar(VARSET_KEY_ICETHICK, XYT, 'Ice thickness', 'm')
    call registVar(VARSET_KEY_SICESURFTEMP, XYT, 'Surface temperature of snow or ice layer', 'degC')
    call registVar(VARSET_KEY_SICETEMP, XYZ2T, 'Temperature of ice layer', 'degC')

  contains

    subroutine registVar(varName, dim, long_name, units)
      character(*), intent(in) :: varName
      character(TOKEN), intent(in) :: dim(:)
      character(*), intent(in) :: long_name, units


      call HistoryAddVariable(varName, dim, &
           & trim(long_name)//' (at t)', units, xtype='double', history=gthst_rst)
      call HistoryAddVariable(trim(varName)//'B', dim, &
           & trim(long_name)//' (at t-1)', units, xtype='double', history=gthst_rst)

    end subroutine registVar
    
  end subroutine createRestartFile

end module RestartDataFileSet_mod
