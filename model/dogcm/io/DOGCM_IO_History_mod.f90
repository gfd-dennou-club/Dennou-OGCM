!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_IO_History_mod 

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_historyauto, only:  &
       & HistoryAutoPut,        &
       & HistoryAutoAddVariable
  
  !* Dennou-OGCM
  
  use DOGCM_Admin_Grid_mod, only: &
       & IM, IS, IE, IA, &
       & JM, JS, JE, JA, &
       & KM, KS, KE, KA, &
       & xyz_Z,          &
       & AXIS_INFO,                                      &
       & IAXIS_info, JAXIS_info, KAXIS_info, TAXIS_info, &
       & x_CI, y_CJ, z_CK, x_FI, y_FJ, z_FK,             &
       & x_IAXIS_Weight, y_JAXIS_Weight, z_KAXIS_Weight
  
  use DOGCM_Admin_TInteg_mod, only: &
       & CurrentTime, EndTime,       &
       & RestartTime,                &
       & TIMELV_ID_N, TIMELV_ID_B
  
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_IO_History_Init, DOGCM_IO_History_Final
  public :: DOGCM_IO_History_Create
  public :: DOGCM_IO_History_RegistVar
  public :: DOGCM_IO_History_IsOutputTiming

  interface DOGCM_IO_History_HistPut
     module procedure DOGCM_IO_History_HistPut0D
     module procedure DOGCM_IO_History_HistPut1D
     module procedure DOGCM_IO_History_HistPut2D
     module procedure DOGCM_IO_History_HistPut3D
  end interface DOGCM_IO_History_HistPut
  public :: DOGCM_IO_History_HistPut
  
  
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
  character(*), parameter:: module_name = 'DOGCM_IO_History_mod' !< Module Name

  integer :: HstDimsID_I = 1
  integer :: HstDimsID_J = 2
  integer :: HstDimsID_K = 3
  integer :: HstDimsID_T = 4  
  character(TOKEN) :: HstDimsList(4)
  
contains

  !>
  !!
  !!
  subroutine DOGCM_IO_History_Init(configNmlFileName)


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
    call read_nmlData(configNmlFileName, outputFileName)


  end subroutine DOGCM_IO_History_Init

  !>
  !!
  !!
  subroutine DOGCM_IO_History_Final()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_IO_History_Final


  
  logical function DOGCM_IO_History_isOutputTiming(CurrentTimeSec)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: CurrentTimeSec

    
    ! 実行文; Executable statement
    !
    
    DOGCM_IO_History_isOutputTiming &
         & = ( mod(CurrentTimeSec, OutputIntrvalSec) == 0 )

  end function DOGCM_IO_History_isOutputTiming

  !----------------------------------------------------------------

  !> @brief 
  !!
  !!
  subroutine DOGCM_IO_History_HistPut0D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var
    
    ! 実行文; Executable statement
    !    
    call HistoryAutoPut(CurrentTime, varName, var)

  end subroutine DOGCM_IO_History_HistPut0D
  
  !> @brief 
  !!
  !!
  subroutine DOGCM_IO_History_HistPut1D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:)
    
    ! 実行文; Executable statement
    !    
    call HistoryAutoPut(CurrentTime, varName, var)

  end subroutine DOGCM_IO_History_HistPut1D

  !> @brief 
  !!
  !!
  subroutine DOGCM_IO_History_HistPut2D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:,:)
    
    ! 実行文; Executable statement
    !    
    call HistoryAutoPut(CurrentTime, varName, var)

  end subroutine DOGCM_IO_History_HistPut2D

  !> @brief 
  !!
  !!
  subroutine DOGCM_IO_History_HistPut3D(varName, var)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: varName
    real(DP), intent(in) :: var(:,:,:)
    
    ! 実行文; Executable statement
    !    
    call HistoryAutoPut(CurrentTime, varName, var)
    
  end subroutine DOGCM_IO_History_HistPut3D
  
  !-----------------------------------------------------------------
  
  !> @brief 
  !!
  !!
  subroutine DOGCM_IO_History_Create( configNmlFileName )

    ! モジュール引用; Use statement
    !

    use gtool_historyauto, only: &
         & HistoryAutoCreate,      &
         & HistoryAutoAddAttr,     &
         & HistoryAutoAddWeight,   &
         & HistoryAutoPutAxis         

    
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

    character(STRING) :: HstLongNameList( size(HstDimsList) )
    character(TOKEN) :: HstUnitsList( size(HstDimsList) )
    character(TOKEN) :: HstTypeList( size(HstDimsList) )
    
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

    ! Initialize gtool objects
    !
    HstDimsList(1) = IAXIS_info%name
    HstDimsList(2) = JAXIS_info%name
    HstDimsList(3) = KAXIS_info%name
    HstDimsList(4) = TAXIS_info%name

    
    HstLongNameList(1) = IAXIS_info%long_name
    HstLongNameList(2) = JAXIS_info%long_name
    HstLongNameList(3) = KAXIS_info%long_name
    HstLongNameList(4) = TAXIS_info%long_name

    HstUnitsList(1) = IAXIS_info%units
    HstUnitsList(2) = JAXIS_info%units
    HstUnitsList(3) = KAXIS_info%units
    HstUnitsList(4) = TAXIS_info%units

    call HistoryAutoCreate( &                            ! ヒストリー作成
         & title  = 'DOGCM Output',             &
         & source = 'DOGCM Output',                                         &
         & institution='GFD_Dennou Club OGCM project',                      &
         & dims=HstDimsList, dimsizes=(/ IM, JM, KM, 0 /),                  &
         & longnames=HstLongNameList, units=HstUnitsList,                   &
         & origin=RestartTime, interval=outputIntrvalSec, terminus=EndTime, &
         & namelist_filename=configNmlFileName )    
    
    ! Regist coordinates

    call MessageNotify('M', module_name, 'Regist axis information..')
    call regist_axis(IAXIS_info, x_CI(IS:IE), x_IAXIS_Weight(IS:IE) )
    call regist_axis(JAXIS_info, y_CJ(JS:JE), y_JAXIS_Weight(JS:JE) )
    call regist_axis(KAXIS_info, z_CK(KS:KE), z_KAXIS_Weight(KS:KE) )    

  contains
    subroutine regist_axis(axis, cell_pos, cell_weight)
      type(AXIS_INFO), intent(in) :: axis
      real(DP), intent(in) :: cell_pos(:)
      real(DP), intent(in) :: cell_weight(:)
      
      call HistoryAutoAddAttr(axis%name, 'standard_name', axis%long_name)
      call HistoryAutoPutAxis(axis%name, cell_pos)
      call HistoryAutoAddWeight(axis%name, cell_weight, axis%weight_units, xtype='double')
      
    end subroutine regist_axis
    
  end subroutine DOGCM_IO_History_create

  !---------------------------------------
  
  subroutine DOGCM_IO_History_RegistVar( &
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
    
    select case (varDimsName)
    case ('T')
       varDimsLen = 1
       varDims(1:varDimsLen) = (/ HstDimsList(HstDimsID_T) /)
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
    
    call HistoryAutoAddVariable( varname=trim(varName), &
         & dims=varDims(1:varDimsLen), longname=varLongName, units=trim(varUnits) &
         & )
    
  end subroutine DOGCM_IO_History_RegistVar

  !---------------------------------------
  
  subroutine read_nmlData( configNmlFileName, &
       & outputFileName )

    ! モジュール引用; Use statement
    !
    use dc_calendar, only: &
         & DCCalConvertByUnit

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
    character(STRING), intent(out) :: outputFileName

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
    real(DP) :: OriginValue
    character(TOKEN) :: OriginUnit
    real(DP) :: TerminusValue
    character(TOKEN) :: TerminusUnit

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /gtool_historyauto_nml/ &
         & IntValue, IntUnit, Name, FilePrefix,                  &
         & OriginValue, OriginUnit, TerminusValue, TerminusUnit


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
       do while ( trim(pos_nml) /= 'APPEND' .and. iostat_nml == 0 ) 
          read( unit_nml, &           ! (in)
               & nml = gtool_historyauto_nml, iostat = iostat_nml )   ! (out)
          inquire( unit_nml, &      !(in)
               & position=pos_nml ) !(out)
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
  
end module DOGCM_IO_History_mod

