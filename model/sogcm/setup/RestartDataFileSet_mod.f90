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
       & xyz_UN, xyz_VN, xyz_PTempEddN, xyz_SaltN, &
       & xyz_UB, xyz_VB, xyz_PTempEddB, xyz_SaltB

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

    call HistoryPut('UN', xyz_UN, history=gthst_rst)
    call HistoryPut('VN', xyz_VN, history=gthst_rst)
    call HistoryPut('PTempEddN', xyz_PTempEddN, history=gthst_rst)
    call HistoryPut('SaltN', xyz_SaltN, history=gthst_rst)

    call HistoryPut('UB', xyz_UB, history=gthst_rst)
    call HistoryPut('VB', xyz_VB, history=gthst_rst)
    call HistoryPut('PTempEddB', xyz_PTempEddB, history=gthst_rst)
    call HistoryPut('SaltB', xyz_SaltB, history=gthst_rst)

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

    !
    call HistoryGet( InputFileName, 'UN', range=timeRange, &
         & array=xyz_UN )
    call HistoryGet( InputFileName, 'VN', range=timeRange, &
         & array=xyz_VN )
    call HistoryGet( InputFileName, 'PTempEddN', range=timeRange, &
         & array=xyz_PTempEddN )
    call HistoryGet( InputFileName, 'SaltN', range=timeRange, &
         & array=xyz_SaltN )

    call HistoryGet( InputFileName, 'UB', range=timeRange, &
         & array=xyz_UB )
    call HistoryGet( InputFileName, 'VB', range=timeRange, &
         & array=xyz_VB )
    call HistoryGet( InputFileName, 'PTempEddB', range=timeRange, &
         & array=xyz_PTempEddB )
    call HistoryGet( InputFileName, 'SaltB', range=timeRange, &
         & array=xyz_SaltB )


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
    
    
    ! 実行文; Executable statement
    !
    
    !

    call HistoryCreate( &
         & file=OutputFileName, title=' ', source=' ',  institution=' ',      &
         & dims=(/ 'lon    ', 'lat    ', 'sig    ', 'timestr', 'time   ' /),  &
         & dimsizes=(/ iMax, jMax, kMax+1, TOKEN, 0 /),                   &
         & longnames =(/ 'longitude                        ',   &
         &               'latitude                         ',   &
         &               'sigma at collocation points      ',   &
         &               'number of characters for datetime',   &
         &               'time                             '/), &
         & units=(/'degree_east ', 'degree_north', '1           ',          &
         &         '1           ', 'sec         ' /),                      &
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
      & history = gthst_rst )           ! (inout)
    call HistoryPut( &
      & 'lat', xyz_Lat(0,:,1)/PI*180d0, & ! (in)
      & history = gthst_rst )           ! (inout)
    call HistoryPut( &
      & 'sig', g_Sig, &               ! (in)
      & history = gthst_rst )           ! (inout)

    !
    !
    call regist_Variable('UN', 'zonal velociy (at t)', 'm s-1')
    call regist_Variable('VN', 'meridional velociy (at t)', 'm s-1')
    call regist_Variable('PTempEddN', 'eddy component of potential temerature (at t)', 'K')
    call regist_Variable('SaltN', 'salinity (at t)', 'psu')

    call regist_Variable('UB', 'zonal velociy (at t-1)', 'm s-1')
    call regist_Variable('VB', 'meridional velociy (at t-1)', 'm s-1')
    call regist_Variable('PTempEddB', 'eddy component of potential temerature (at t-1)', 'K')
    call regist_Variable('SaltB', 'salinity (at t-1)', 'psu')

  contains

    subroutine regist_Variable(varName, long_name, units)
      character(*), intent(in) :: varName, long_name, units

      call HistoryAddVariable(varName, (/ 'lon ', 'lat ', 'sig ', 'time' /), &
           & long_name, units, xtype='double', history=gthst_rst)
    end subroutine regist_Variable
    
  end subroutine createRestartFile

end module RestartDataFileSet_mod
