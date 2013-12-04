program diagVar_main

  use dc_types
  use dc_message

  use Constants_mod
  use GridSet_mod
  use SpmlUtil_mod
  use VariableSet_mod

  use DiagVarSet_mod
  use DiagVarFileSet_mod
  use DiagVarEval_mod

  use gtool_history
   use dc_calendar, only: &
         & DCCalCreate, DCCalDateCreate, &
         & DCCalConvertByUnit, DCCalDateDifference, &
         & DCCalDateInquire, DCCalDateEval

  implicit none

  ! 宣言文; Declaration statement
  !

  character(*), parameter :: PROGRAM_NAME = 'diagVar'
  character(*), parameter :: DEFAULT_CONFIG_NML = 'defaultConfigDiagVar.nml' 
  character(TOKEN), pointer :: ogcmOutputVarsName(:) => null()
  character(TOKEN), pointer :: diagVarsName(:) => null()

  real(DP) :: CurrentTimeSec, EndTimeSec, TimeIntSec

  type(gtool_historyauto_info) :: diagVar_gthsInfo
  type(gtool_historyauto_info) :: ogcm_gthsInfo

  real(DP), pointer :: ogcm_outputTime(:) => null()

  
  ! 実行文; Executable statement
  !

  !
  call MessageNotify('M', PROGRAM_NAME, "Setup...")
  call setup()

  !
  write(*,*) ogcmOutputVarsName
  call HistoryGetPointer( &
       & trim(ogcm_gthsInfo%FilePrefix) // trim(ogcmOutputVarsName(1)) // '.nc', &
       & 't', ogcm_outputTime)
  write(*,*) 'outputTime', size(ogcm_outputTime), ogcm_outputTime
  EndTimeSec = DCCalConvertByUnit(ogcm_outputTime(size(ogcm_outputTime)), ogcm_gthsInfo%intUnit, "sec")
  TimeIntSec = DCCalConvertByUnit(diagVar_gthsInfo%intValue, diagVar_gthsInfo%intUnit, 'sec')

  CurrentTimeSec = 0
  do while(CurrentTimeSec <= EndTimeSec)
     call MessageNotify('M', PROGRAM_NAME, &
          & "%f [sec]..", d=(/CurrentTimeSec/) )

     call diagnose_outputVariables()

     CurrentTimeSec = CurrentTimeSec + TimeIntSec
  end do

  !
  call MessageNotify('M', PROGRAM_NAME, "Shutdown...")
  call shutdown()

contains

  !> @brief 
  !!
  !!
  subroutine setup()

    ! モジュール引用; Use statement
    !
    use OptionParser_mod, only: &
         & OptionParser_Init, OptionParser_Final, &
         & OptionParser_GetInfo

 
#ifdef _OPENMP
    use omp_lib
#endif

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    character(STRING) :: configNmlFile
    character(STRING) :: ogcmConfigNmlFile

    integer :: nThread
    integer :: varID

    ! 実行文; Executable statement
    !

    call OptionParser_Init()
    call OptionParser_GetInfo(configNmlFile,   & !(out)
         & defaultConfigNml=DEFAULT_CONFIG_NML ) 
    call OptionParser_Final()
    
    call MessageNotify('M', PROGRAM_NAME, &
         & "Read configure file ' %c ' ..", c1=trim(configNmlFile) )
    call readNml(configNmlFile, &  ! (in)
         & ogcmConfigNmlFile ) ! (out)
    call readOgcmNml(ogcmConfigNmlFile)

    call Constants_Init(ogcmConfigNmlFile)
    call GridSet_Init(ogcmConfigNmlFile)

#ifdef _OPENMP
    !$omp parallel
    !$omp single
    nThread = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    call MessageNotify('M', PROGRAM_NAME, "Execute as Thread Parallel Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, np=nThread)
#else
    call MessageNotify('M', PROGRAM_NAME, "Execute as Serial  Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
#endif


    call GridSet_construct()
    call VariableSet_Init()

    !
    call MessageNotify('M', PROGRAM_NAME, &
         & "Some modules is initialized.. ")
    call DiagVarSet_Init(diagVarsName)
    call DiagVarEval_Init()
    call DiagVarFileSet_Init(configNmlFile, diagVar_gthsInfo)

    !

  end subroutine setup

  !> @brief 
  !!
  !!
  subroutine shutdown()
    
    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    call DiagVarEval_Final()
    call DiagVarFileSet_Final()
    call DiagVarSet_Final()

    if( associated(diagVarsName) ) deallocate(diagVarsName)
    if( associated(ogcmOutputVarsName) ) deallocate(ogcmOutputVarsName)
    call VariableSet_Final()
    call SpmlUtil_Final()
    call GridSet_Final()
    call Constants_Final()

  end subroutine shutdown

  !> @brief 
  !!
  !!
  subroutine readNml(configNml, ogcmConfigNml)

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
    
    use dc_string, only: Split, Replace

    ! 宣言文; Declaration statement
    !    
    character(*), intent(in) :: configNml
    character(*), intent(out) :: ogcmConfigNml

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
    character(STRING) :: Name, FilePrefix

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /diagVar_nml/ &
         & ogcmConfigNml

    namelist /gtool_historyauto_nml/ &
         & IntValue, IntUnit, Name, FilePrefix

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    Name = ''
    FilePrefix = ''

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(configNml) /= '' ) then
       call MessageNotify( 'M',PROGRAM_NAME, "reading namelist '%a'", ca=(/ configNml /))
       call FileOpen( unit_nml, &          ! (out)
            & configNml, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = diagVar_nml, iostat = iostat_nml )   ! (out)

       pos_nml = ''
       do while ( trim(pos_nml) /= 'APPEND' .and. iostat_nml == 0 ) 
          read( unit_nml, &           ! (in)
               & nml = gtool_historyauto_nml, iostat = iostat_nml )   ! (out)
          inquire( unit_nml, & !(in)
               & position=pos_nml ) !(out)
       end do

       close( unit_nml )
    end if

    Name = Replace(Name, " ", "")
    call Split(trim(Name), diagVarsName, ",")

    diagVar_gthsInfo = gtool_historyauto_info( intValue=intValue, intUnit=intUnit, &
         &                                     FilePrefix=FilePrefix, Name=Name )

    ! 印字 ; Print
    !
    call MessageNotify( 'M', PROGRAM_NAME, '----- Initialization Messages -----' )
    call MessageNotify( 'M', PROGRAM_NAME, ' ogcm configure file   = %c', c1=trim(ogcmConfigNml) )
    call MessageNotify( 'M', PROGRAM_NAME, ' diagnostic variables  = %c', c1=Name )
    
  end subroutine readNml

  !> @brief 
  !!
  !!
  subroutine readOgcmNml(ogcmConfigNml)

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
    
    use dc_string, only: Split, Replace

    ! 宣言文; Declaration statement
    !    
    character(*), intent(in) :: ogcmConfigNml

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
    character(STRING) :: Name, FilePrefix

    ! NAMELIST 変数群
    ! NAMELIST group name
    !
    namelist /gtool_historyauto_nml/ &
         & IntValue, IntUnit, Name, FilePrefix

    ! 実行文; Executable statements

    ! デフォルト値の設定
    ! Default values settings
    !
    FilePrefix = ''
    Name = ''

    ! NAMELIST からの入力
    ! Input from NAMELIST
    !
    if ( trim(ogcmConfigNml) /= '' ) then
       call MessageNotify( 'M',PROGRAM_NAME, "reading ogcm namelist '%a'", ca=(/ ogcmConfigNml /))
       call FileOpen( unit_nml, &          ! (out)
            & ogcmConfigNml, mode = 'r' ) ! (in)

       pos_nml = ''
       do while ( trim(pos_nml) /= 'APPEND' .and. iostat_nml == 0 ) 
          read( unit_nml, &           ! (in)
               & nml = gtool_historyauto_nml, iostat = iostat_nml )   ! (out)
          inquire( unit_nml, & !(in)
               & position=pos_nml ) !(out)
       end do
       close( unit_nml )
    end if


    Name = Replace(Name, " ", "")
    call Split(trim(Name), ogcmOutputVarsName, ",")

    ogcm_gthsInfo = gtool_historyauto_info( intValue=intValue, intUnit=intUnit, &
         &                                     FilePrefix=FilePrefix, Name=Name )

    ! 印字 ; Print
    !
    call MessageNotify( 'M', PROGRAM_NAME, '----- Initialization Messages -----' )
    call MessageNotify( 'M', PROGRAM_NAME, ' ogcm output variables   = %c', c1=Name )

  end subroutine readOgcmNml


  !> @brief 
  !!
  !!
  subroutine diagnose_outputVariables()

    ! モジュール引用; Use statement
    !
    use dc_string, only: CPrintf

    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
    integer :: varID
    character(STRING) :: rangeStr
    real(DP) :: CurrentTime
    real(DP) :: xyz_BarocPress(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_Press(0:iMax-1,jMax,0:kMax)

    ! 実行文; Executable statement
    !
    
    CurrentTime = DCCalConvertByUnit(CurrentTimeSec, 'sec', ogcm_gthsInfo%intUnit)
    rangeStr = CPrintf('t=%f', d=(/ currentTime /))
    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // "u.nc", &
         & 'u', xyz_UN, range=rangeStr )
    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // "v.nc", &
         & 'v', xyz_VN, range=rangeStr )
    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // "SurfPress.nc", &
         & 'SurfPress', xy_SurfPress, range=rangeStr )
    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // "BarocPress.nc", &
         & 'BarocPress', xyz_BarocPress, range=rangeStr )

    do varID=1, size(diagVarsName)
       select case(diagVarsName(varID))
          case(DVARKEY_DIV)
             xyz_Div = eval_Div(xyz_UN, xyz_VN)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_Div, var3D=xyz_Div)
          case(DVARKEY_VOR)
             xyz_Vor = eval_Vor(xyz_UN, xyz_VN)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_Vor, var3D=xyz_Vor)
          case (DVARKEY_PRESS)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_PRESS, &
                  & var3D= eval_totPress(xy_SurfPress, xyz_BarocPress) )
       end select
    end do

    
  end subroutine diagnose_outputVariables

end program diagVar_main
