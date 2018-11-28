!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_IO_Restart_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_history, only: &
       & gt_history,         &
       & HistorySetTime,                       &
       & HistoryCreate, HistoryClose,          &
       & HistoryAddVariable, HistoryAddAttr,   &       
       & HistoryPut,                           &
       & HistoryGet, HistoryGetAttr


  !* Dennou-OGCM / SIce

  use DSIce_Admin_TInteg_mod, only: &
       & RestartTime, &
       & TIMELV_ID_N, TIMELV_ID_B, TIMELV_ID_A  

  use DSIce_Admin_Grid_mod, only: &
       & IS, IE, IM, &
       & JS, JE, JM, &
       & KS, KE, KM
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_IO_Restart_Init, DSIce_IO_Restart_Final

  public :: DSIce_IO_Restart_Create
  public :: DSIce_IO_Restart_RegistVar
  public :: DSIce_IO_Restart_Output
  interface DSIce_IO_Restart_HistPut
     module procedure DSIce_IO_Restart_HistPut1D
     module procedure DSIce_IO_Restart_HistPut2D
     module procedure DSIce_IO_Restart_HistPut3D
  end interface DSIce_IO_Restart_HistPut
  public :: DSIce_IO_Restart_HistPut  


  public :: DSIce_IO_Restart_Input  
  interface DSIce_IO_Restart_HistGet
     module procedure DSIce_IO_Restart_HistGet1D
     module procedure DSIce_IO_Restart_HistGet2D
     module procedure DSIce_IO_Restart_HistGet3D
  end interface DSIce_IO_Restart_HistGet
  public :: DSIce_IO_Restart_HistGet
  
  public :: DSIce_IO_Restart_IsOutputTiming
  
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
  character(*), parameter:: module_name = 'DSIce_IO_Restart_mod' !< Module Name
  type(gt_history), save :: gthst_rst

  integer :: HstDimsID_I = 1
  integer :: HstDimsID_J = 2
  integer :: HstDimsID_K = 3
  integer :: HstDimsID_S = 4  
  integer :: HstDimsID_T = 5
  character(TOKEN) :: HstDimsList(5)
  
  character(STRING) :: InputTimeRange
  
contains

  !>
  !!
  !!
  subroutine DSIce_IO_Restart_Init(configNmlFileName)

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
 

  end subroutine DSIce_IO_Restart_Init

  !>
  !!
  !!
  subroutine DSIce_IO_Restart_Final()

    use gtool_history, only: HistoryClose

    ! 実行文; Executable statements
    !

    call HistoryClose(history=gthst_rst)

  end subroutine DSIce_IO_Restart_Final

    !-----------------------------------------

  logical function DSIce_IO_Restart_isOutputTiming(CurrentTimeSec)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: CurrentTimeSec

    
    ! 実行文; Executable statement
    !
    
    DSIce_IO_Restart_isOutputTiming &
         & = ( mod(CurrentTimeSec, IntTime) == 0 )

  end function DSIce_IO_Restart_isOutputTiming

  !-----------------------------------------------
  !----------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DSIce_IO_Restart_HistPut1D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:)
    
    ! 実行文; Executable statement
    !    
    call HistoryPut(varName, var, history=gthst_rst)

  end subroutine DSIce_IO_Restart_HistPut1D

  !> @brief 
  !!
  !!
  subroutine DSIce_IO_Restart_HistPut2D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:,:)
    
    ! 実行文; Executable statement
    !    
    call HistoryPut(varName, var, history=gthst_rst)

  end subroutine DSIce_IO_Restart_HistPut2D

  !> @brief 
  !!
  !!
  subroutine DSIce_IO_Restart_HistPut3D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:,:,:)
    
    ! 実行文; Executable statement
    !    
    call HistoryPut(varName, var, history=gthst_rst)

  end subroutine DSIce_IO_Restart_HistPut3D

  !-----------------------------------------------------------------
  
  subroutine DSIce_IO_Restart_RegistVar( &
       & varName, varDimsName, varLongName, varUnits  & ! (in)
       & )

    
    character(*), intent(in) :: varName
    character(*), intent(in) :: varDimsName
    character(*), intent(in) :: varLongName
    character(*), intent(in) :: varUnits

    character(TOKEN) :: varDims(size(HstDimsList))
    integer :: varDimsLen
    
    select case (varDimsName)
    case ('K')
       varDimsLen = 1       
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_K) /)
    case ('IJ')
       varDimsLen = 2
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_I),  HstDimsList(HstDimsID_J) /)
    case ('IJK')
       varDimsLen = 3
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_I),  HstDimsList(HstDimsID_J), HstDimsList(HstDimsID_K) /)
    case ('IJT')
       varDimsLen = 3
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_I),  HstDimsList(HstDimsID_J), HstDimsList(HstDimsID_T) /)
    case ('IJKT')   
       varDimsLen = 4
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_I),  HstDimsList(HstDimsID_J), HstDimsList(HstDimsID_K), &
            &                     HstDimsList(HstDimsID_T) /)
    case default
       call MessageNotify( 'E', module_name, &
            & "Specified varDimsName(=%a) is not supported. Check!", ca=(/ trim(varDimsName) /) )
    end select

    call HistoryAddVariable( varname=trim(varName), &
         & dims=varDims(1:varDimsLen), longname=varLongName, units=trim(varUnits), xtype='double', history=gthst_rst &
         & )
    
  end subroutine DSIce_IO_Restart_RegistVar

  !-----------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DSIce_IO_Restart_Output()
    
    ! モジュール引用; Use statement
    !

    use dc_calendar, only: &
         & DCCalDateInquire

    use DSIce_Admin_TInteg_mod, only: &
         & CurrentTime, InitDate

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    character(TOKEN) :: dateStr
    
    ! 実行文; Executable statement
    !
    if( .not. DSIce_IO_Restart_isOutputTiming( CurrentTime ) ) return

    call MessageNotify('M', module_name, 'Output data to restart..')


    ! Output data of some variables for restarting a simulation.
    
    call HistorySetTime(timed=CurrentTime, history=gthst_rst)
    call DCCalDateInquire(dateStr, &
         & elapse_sec=CurrentTime, date=InitDate)
    
    call HistoryPut('datetime', dateStr, history=gthst_rst)
    
  end subroutine DSIce_IO_Restart_Output

  subroutine DSIce_IO_Restart_HistGet1D(       &
       & varName,                               & ! (in)
       & var                                    & ! (in)
       & )

    real(DP), intent(inout) :: var(:)
    character(*), intent(in) :: varName

    call HistoryGet(InputFileName, trim(varName), range=InputTimeRange, &
         & array=var )
    
  end subroutine DSIce_IO_Restart_HistGet1D

  subroutine DSIce_IO_Restart_HistGet2D(       &
       & varName,                               & ! (in)
       & var                                    & ! (in)
       & )

    real(DP), intent(inout) :: var(:,:)
    character(*), intent(in) :: varName

    call HistoryGet(InputFileName, trim(varName), range=InputTimeRange, &
         & array=var )
    
  end subroutine DSIce_IO_Restart_HistGet2D

  subroutine DSIce_IO_Restart_HistGet3D(       &
       & varName,                               & ! (in)
       & var                                    & ! (in)
       & )

    real(DP), intent(inout) :: var(:,:,:)
    character(*), intent(in) :: varName

    call HistoryGet(InputFileName, trim(varName), range=InputTimeRange, &
         & array=var )
    
  end subroutine DSIce_IO_Restart_HistGet3D
  
  !> @brief 
  !!
  !!
  subroutine DSIce_IO_Restart_Input()
    
    
    ! モジュール引用; Use statement
    !
    use dc_string, only: &
         & toChar    
    
    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
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
    
    InputTimeRange = 'time=' // toChar(RestartTime)

  end subroutine DSIce_IO_Restart_Input


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
    namelist /SeaIce_IO_Restart_nml/ &
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
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &                      ! (in)
            & nml = SeaIce_IO_Restart_nml, &  ! (out)
            & iostat = iostat_nml )           ! (out)
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
    call MessageNotify( 'M', module_name, ' IntTime        = %f [sec]', d=(/ IntTime  /) )

  end subroutine read_nmlData

  !> @brief 
  !!
  !!
  subroutine DSIce_IO_Restart_Create( configNmlFileName )
    
    !
    !
    use DSIce_Admin_Grid_mod, only: &
         & AXIS_INFO, &
         & IAXIS_info, JAXIS_info, KAXIS_info, TAXIS_info, &
         & x_CI, y_CJ, z_CK, x_FI, y_FJ, z_FK, &
         & x_IAXIS_Weight, y_JAXIS_Weight, z_KAXIS_Weight

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName
    
    ! 局所変数
    ! Local variables
    !

    character(TOKEN) :: dims_K(1)
    character(TOKEN) :: dims_IJ(2)
    character(TOKEN) :: dims_IJK(3)
    character(TOKEN) :: dims_IJT(3)
    character(TOKEN) :: dims_KT(2)
    character(TOKEN) :: dims_IJKT(4)
    character(TOKEN) :: dims_IJKAT(5)

    character(STRING) :: HstLongNameList(5)
    character(TOKEN) :: HstUnitsList(5)
    character(TOKEN) :: HstTypeList(5)
    
    ! 実行文; Executable statement
    !

    ! Set arrays storing the name of axises
    !

    dims_K = (/ KAXIS_info%name /)
    dims_KT = (/ KAXIS_info%name, TAXIS_info%name /)
    dims_IJ = (/ IAXIS_info%name, JAXIS_info%name /)
    dims_IJK = (/ IAXIS_info%name, JAXIS_info%name, KAXIS_info%name /)
    dims_IJT = (/ IAXIS_info%name, JAXIS_info%name, TAXIS_info%name /)
    dims_IJKT = (/ IAXIS_info%name, JAXIS_info%name, KAXIS_info%name, TAXIS_info%name /)
    dims_IJKT = (/ IAXIS_info%name, JAXIS_info%name, KAXIS_info%name, TAXIS_info%name /)
    
    !
    HstDimsList(1) = IAXIS_info%name
    HstDimsList(2) = JAXIS_info%name
    HstDimsList(3) = KAXIS_info%name
    HstDimsList(4) = 'timestr'
    HstDimsList(5) = TAXIS_info%name

    HstLongNameList(1) = IAXIS_info%long_name
    HstLongNameList(2) = JAXIS_info%long_name
    HstLongNameList(3) = KAXIS_info%long_name
    HstLongNameList(4) = 'number of characters for datetime    '
    HstLongNameList(5) = TAXIS_info%long_name

    HstUnitsList(1) = IAXIS_info%units
    HstUnitsList(2) = JAXIS_info%units
    HstUnitsList(3) = KAXIS_info%units
    HstUnitsList(4) = '1'
    HstUnitsList(5) = TAXIS_info%units

    HstTypeList(1) = 'double'
    HstTypeList(2) = 'double'
    HstTypeList(3) = 'double'
    HstTypeList(4) = 'int'
    HstTypeList(5) = 'double'
    
    call HistoryCreate( &
         & file=OutputFileName, title=' ', source=' ',  institution=' ',      &
         & dims=HstDimsList, dimsizes=(/ IM, JM, KM, TOKEN, 0 /),             &
         & longnames=HstLongNameList, units=HstUnitsList, xtypes=HstTypeList, &
         & origind=RestartTime, intervald=IntTime, history=gthst_rst )
 
    ! 座標データの設定
    ! Axes data settings

    call regist_axis(IAXIS_info, x_CI(IS:IE), x_IAXIS_Weight(IS:IE))
    call regist_axis(JAXIS_info, y_CJ(JS:JE), y_JAXIS_Weight(JS:JE))
    call regist_axis(KAXIS_info, z_CK(KS:KE), z_KAXIS_Weight(KS:KE))

    ! 年月日時分秒形式の日時情報変数の設定
    ! Set a date information variable with year-month-day hour:minute:second format
    !    
    call HistoryAddVariable( &
         & varname='datetime', dims=(/HstDimsList(HstDimsID_S), HstDimsList(HstDimsID_T)/), &
         & longname='time represented as strings', units='1', xtype='char', &
         & history=gthst_rst )
    

  contains
    subroutine regist_axis(axis, cell_pos, cell_weight)
      type(AXIS_INFO), intent(in) :: axis
      real(DP), intent(in) :: cell_pos(:)
      real(DP), intent(in) :: cell_weight(:)
      
      call HistoryAddAttr(axis%name, 'standard_name', axis%long_name, &
           & history = gthst_rst )
      call HistoryPut(axis%name, cell_pos, &
           & history = gthst_rst )      
    end subroutine regist_axis

  end subroutine DSIce_IO_Restart_Create

end module DSIce_IO_Restart_mod
