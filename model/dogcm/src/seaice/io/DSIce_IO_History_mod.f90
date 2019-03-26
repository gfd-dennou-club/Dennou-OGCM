!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_IO_History_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5

  use dc_types, only:               &
       & DP, TOKEN, STRING

  use dc_string, only:              &
       & StrInclude
  
  use dc_message, only:             &
       & MessageNotify

  use gtool_history, only: &
       & gt_history,         &
       & HistorySetTime,                       &
       & HistoryCreate, HistoryClose,          &
       & HistoryAddVariable, HistoryAddAttr,   &       
       & HistoryPut                            

  use dc_calendar, only: &
       & DCCalConvertByUnit
  
  !* Dennou-OGCM / SIce

  use DSIce_Admin_Grid_mod, only: &
       & IM, IS, IE, IA, &
       & JM, JS, JE, JA, &
       & KM, KS, KE, KA

  use DSIce_Admin_TInteg_mod, only: &
       & CurrentTime, RestartTime, &
       & IntegTime, EndTime
  
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_IO_History_Init, DSIce_IO_History_Final
  public :: DSIce_IO_History_Create
  public :: DSIce_IO_History_RegistVar
  public :: DSIce_IO_History_Output
  public :: DSIce_IO_History_IsOutputTiming

  interface DSIce_IO_History_HistPut
     module procedure DSIce_IO_History_HistPut1D
     module procedure DSIce_IO_History_HistPut2D
     module procedure DSIce_IO_History_HistPut3D
  end interface DSIce_IO_History_HistPut
  public :: DSIce_IO_History_HistPut
  
  
  ! 公開変数
  ! Public variable
  !  
  character(TOKEN), public, save :: FilePrefix
  real(DP), public, save :: OutputIntrvalSec
  character(TOKEN), public, save :: OutputIntUnit
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_IO_History_mod' !< Module Name

  integer :: HstDimsID_I  = 1
  integer :: HstDimsID_J  = 2
  integer :: HstDimsID_Jm = 3 
  integer :: HstDimsID_K  = 4
  integer :: HstDimsID_T  = 5  
  character(TOKEN) :: HstDimsList(5)

  integer, public :: HstVarsNum
  integer, parameter :: HstVarsNumMax = 50
  character(TOKEN) :: HstVarsList(HstVarsNumMax)
  logical :: FlagTimeAverageList(HstVarsNumMax)
  
  logical :: isHistoryCreated

  type(gt_history), save :: hst_seaice

  real(DP) :: CurrentTimeOutputUnit
  
contains

  !>
  !!
  !!
  subroutine DSIce_IO_History_Init(configNmlFileName)


    ! モジュール引用; Use statement
    !
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlFileName

    ! 局所変数
    ! Local variable
    !
    character(STRING) :: outputFileName

    ! 実行文; Executable statements
    !

    !
    HstVarsNum = 0
    
    call read_nmlData( configNmlFileName )

    isHistoryCreated = .false.
    
  end subroutine DSIce_IO_History_Init

  !>
  !!
  !!
  subroutine DSIce_IO_History_Final()

    ! 実行文; Executable statements
    !
    if (isHistoryCreated) then
      call HistoryClose( hst_seaice )
    end if
    
  end subroutine DSIce_IO_History_Final


  
  logical function DSIce_IO_History_isOutputTiming(CurrentTimeSec)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: CurrentTimeSec

    
    ! 実行文; Executable statement
    !

    DSIce_IO_History_isOutputTiming &
         & = ( mod(CurrentTimeSec, OutputIntrvalSec) == 0 )

  end function DSIce_IO_History_isOutputTiming

  !----------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DSIce_IO_History_HistPut1D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:)
    
    ! 実行文; Executable statement
    !

    logical :: isOutputTiming
    
    if(is_HistoryPut_called(varName, isOutputTiming)) then
       call HistoryPut( varName, var, history=hst_seaice, &
            & timed=CurrentTimeOutputUnit,               &
            & time_average_store=.not.isOutputTiming )
    end if

  end subroutine DSIce_IO_History_HistPut1D

  !> @brief 
  !!
  !!
  subroutine DSIce_IO_History_HistPut2D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:,:)
    
    ! 実行文; Executable statement
    !    

    logical :: isOutputTiming
    
    if(is_HistoryPut_called(varName, isOutputTiming)) then
       call HistoryPut( varName, var, history=hst_seaice, &
            & timed=CurrentTimeOutputUnit,               &
            & time_average_store=.not.isOutputTiming )
    end if

  end subroutine DSIce_IO_History_HistPut2D

  !> @brief 
  !!
  !!
  subroutine DSIce_IO_History_HistPut3D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:,:,:)

    ! 実行文; Executable statement
    !
    
    logical :: isOutputTiming
    
    if(is_HistoryPut_called(varName, isOutputTiming)) then
       call HistoryPut( varName, var, history=hst_seaice, &
            & timed=CurrentTimeOutputUnit,               &
            & time_average_store=.not.isOutputTiming )
    end if

  end subroutine DSIce_IO_History_HistPut3D

  function is_HistoryPut_called(varName, isOutputTiming) result(ret)
    character(*), intent(in) :: varName
    logical :: isOutputTiming
    logical :: ret
    
    integer :: varListID

    ret = .false.

    varListID = get_VarListID(varName)
    if (varListID < 0) return

    isOutputTiming = DSIce_IO_History_isOutputTiming(CurrentTime)
    if ( .not. FlagTimeAverageList(varListID) .and..not. isOutputTiming ) &
         return

    ret = .true.
    return
    
  end function is_HistoryPut_called
  
  !-----------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DSIce_IO_History_Output()

    ! モジュール引用; Use statement
    !

    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    ! 実行文; Executable statement
    !


    CurrentTimeOutputUnit =  DCCalConvertByUnit(CurrentTime, "sec", OutputIntUnit)

    if(  DSIce_IO_History_isOutputTiming(CurrentTime) ) then
       call HistorySetTime( &
            & timed=CurrentTimeOutputUnit,                                 & 
            & history=hst_seaice                                           &
            & )
    end if
    
  end subroutine DSIce_IO_History_Output

  !---------------------------------------
  
  subroutine DSIce_IO_History_RegistVar( &
       & varName, varDimsName, varLongName, varUnits  & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    character(*), intent(in) :: varName
    character(*), intent(in) :: varDimsName
    character(*), intent(in) :: varLongName
    character(*), intent(in) :: varUnits

    ! 局所変数
    ! Local variables
    !    
    character(TOKEN) :: varDims(size(HstDimsList))
    integer :: varDimsLen

    integer :: VarListID

    varListID = get_VarListID(varName)
    if (varListID < 0) return
    
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
    case ('IJmT')
       varDimsLen = 3
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_I),  HstDimsList(HstDimsID_Jm), HstDimsList(HstDimsID_T) /)
    case ('IJKT')   
       varDimsLen = 4
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_I),  HstDimsList(HstDimsID_J), HstDimsList(HstDimsID_K), &
            &                     HstDimsList(HstDimsID_T) /)
    case default
       call MessageNotify( 'E', module_name, &
            & "Specified varDimsName(=%a) is not supported. Check!", ca=(/ trim(varDimsName) /) )
    end select

    call HistoryAddVariable( varname=trim(varName), &
         & dims=varDims(1:varDimsLen), longname=varLongName, units=trim(varUnits), &
         & history=hst_seaice, time_average=FlagTimeAverageList(VarListID)         &
         & )

  end subroutine DSIce_IO_History_RegistVar
  
  !------------------------------------------------------
  
  subroutine DSIce_IO_History_Create( configNmlFileName )

    ! モジュール引用; Use statement
    !
    use DSIce_Admin_Grid_mod, only: &
         & AXIS_INFO, &
         & IAXIS_info, JAXIS_info, JAXIS_half_info, KAXIS_info, TAXIS_info, &
         & x_CI, y_CJ, z_CK, x_FI, y_FJ, z_FK,                              &
         & x_IAXIS_Weight, y_JAXIS_Weight, y_JAXIS_half_Weight, z_KAXIS_Weight

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

    character(STRING) :: HstLongNameList(5)
    character(TOKEN) :: HstUnitsList(5)
  
    ! 実行文; Executable statements
    !

    ! Set arrays storing the name of axises
    !
    dims_K = (/ KAXIS_info%name /)
    dims_KT = (/ KAXIS_info%name, TAXIS_info%name /)
    dims_IJ = (/ IAXIS_info%name, JAXIS_info%name /)
    dims_IJK = (/ IAXIS_info%name, JAXIS_info%name, KAXIS_info%name /)
    dims_IJT = (/ IAXIS_info%name, JAXIS_info%name, TAXIS_info%name /)
    dims_IJKT = (/ IAXIS_info%name, JAXIS_info%name, KAXIS_info%name, TAXIS_info%name /)

    ! Initialize gtool objects
    !
    HstDimsList(HstDimsID_I) = IAXIS_info%name
    HstDimsList(HstDimsID_J) = JAXIS_info%name
    HstDimsList(HstDimsID_Jm) = JAXIS_half_info%name
    HstDimsList(HstDimsID_K) = KAXIS_info%name
    HstDimsList(HstDimsID_T) = TAXIS_info%name

    HstLongNameList(HstDimsID_I) = IAXIS_info%long_name
    HstLongNameList(HstDimsID_J) = JAXIS_info%long_name
    HstLongNameList(HstDimsID_Jm) = JAXIS_half_info%long_name
    HstLongNameList(HstDimsID_K) = KAXIS_info%long_name
    HstLongNameList(HstDimsID_T) = TAXIS_info%long_name

    HstUnitsList(HstDimsID_I) = IAXIS_info%units
    HstUnitsList(HstDimsID_J) = JAXIS_info%units
    HstUnitsList(HstDimsID_Jm) = JAXIS_half_info%units
    HstUnitsList(HstDimsID_K) = KAXIS_info%units
    HstUnitsList(HstDimsID_T) = outputIntUnit    !TAXIS_info%units

    call HistoryCreate( &                            ! ヒストリー作成
         & file   = trim(FilePrefix) // "history_sice.nc",                  & 
         & title  = 'DSIce Output',                                         &
         & source = 'DSIce Output',                                         &
         & institution='GFD_Dennou Club OGCM and sea-ice model  project',   &
         & dims=HstDimsList, dimsizes=(/ IM, JM, JM+1, KM, 0 /),            &
         & longnames=HstLongNameList, units=HstUnitsList,                   &
         & origind=DCCalConvertByUnit(RestartTime, "sec", OutputIntUnit),   &
         & history=hst_seaice )    

    isHistoryCreated = .true.

    ! 座標データの設定
    ! Axes data settings

    call regist_axis(IAXIS_info, x_CI(IS:IE), x_IAXIS_Weight(IS:IE))
    call regist_axis(JAXIS_info, y_CJ(JS:JE), y_JAXIS_Weight(JS:JE))
    call regist_axis(JAXIS_half_info, y_FJ(JS-1:JE), y_JAXIS_half_Weight(JS-1:JE))
    call regist_axis(KAXIS_info, z_CK(KS:KE), z_KAXIS_Weight(KS:KE))

  contains
    subroutine regist_axis(axis, cell_pos, cell_weight)
      type(AXIS_INFO), intent(in) :: axis
      real(DP), intent(in) :: cell_pos(:)
      real(DP), intent(in) :: cell_weight(:)
      
      call HistoryAddAttr(axis%name, 'standard_name', axis%long_name, &
           & history = hst_seaice )
      call HistoryPut(axis%name, cell_pos, &
           & history = hst_seaice )      
    end subroutine regist_axis
    
  end subroutine DSIce_IO_History_Create

  !---------------------------------------

  function get_VarListID(varname) result(varListID)
    character(*), intent(in) :: varname
    integer :: varListID

    integer :: n

    varListID = -1
    do n = 1,  HstVarsNum
       if(trim(varName) == trim(HstVarsList(n))) then
          varListID = n; exit
       end if
    end do

  end function get_VarListID
  
  subroutine read_nmlData( configNmlFileName )

    ! モジュール引用; Use statement
    !

    use dc_string, only: &
         & Replace, Split
    
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

    character(TOKEN) :: pos_nml

    real(DP) :: IntValue
    character(TOKEN) :: IntUnit
    character(STRING) :: Name
    logical :: TimeAverage
    
    character(TOKEN), pointer :: HstVarsList_ptr(:)
    integer :: n
    
    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /seaice_io_history_nml/ &
         & IntValue, IntUnit, FilePrefix, &
         & Name, TimeAverage


    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    IntValue = 0d0
    IntUnit = "sec"
    FilePrefix = ''
    
    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNmlFileName) /= '' ) then
       call MessageNotify( 'M', module_name, "reading namelist '%a'", ca=(/ configNmlFileName /))
       call FileOpen( unit_nml, &             ! (out)
            & configNmlFileName, mode = 'r' ) ! (in)

       pos_nml = ''; iostat_nml = 0
       do while ( .true. )

          TimeAverage = .false.
          Name        = ""
          
          read( unit_nml, &           ! (in)
               & nml = seaice_io_history_nml, iostat = iostat_nml )   ! (out)
          inquire( unit_nml, &      !(in)
               & position=pos_nml ) !(out)
          if( trim(pos_nml) == 'APPEND' .or. iostat_nml /= 0 ) exit
          

          HstVarsList_ptr => null()
          if (len(trim(Name)) > 0) then
               Name = Replace( trim(Name), " ", "", .true.)
               call Split( trim(Name), HstVarsList_ptr, ",")
          end if

          if (associated(HstVarsList_ptr)) then
             do n = 1, size(HstVarsList_ptr)
                if ( .not. StrInclude( HstVarsList(1:HstVarsNum), trim(HstVarsList_ptr(n)) )) then
                   HstVarsNum = HstVarsNum + 1
                   if (HstVarsNum > HstVarsNumMax) then
                      call MessageNotify( 'E', module_name, "Exceed the number of variables allowed to regist. Check!")
                   end if
                   HstVarsList(HstVarsNum) = trim(HstVarsList_ptr(n))
                   FlagTimeAverageList(HstVarsNum) = TimeAverage
                end if
             end do
             deallocate( HstVarsList_ptr )                             
          end if
          
       end do
       close( unit_nml )
    end if

    outputIntUnit = IntUnit
    OutputIntrvalSec = DCCalConvertByUnit(IntValue, IntUnit, "sec")
 
    ! 印字 ; Print
    !
    call MessageNotify( 'M', module_name, '----- Initialization Messages -----' )
    call MessageNotify( 'M', module_name, 'Output Interval %f [%c]', d=(/IntValue/), c1=trim(IntUnit) )

  end subroutine read_nmlData
  
end module DSIce_IO_History_mod

