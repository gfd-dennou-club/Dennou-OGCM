program diagVar_main

  ! �⥸�塼�����; Use statements
  !

  !* gtool
  !
  use dc_types
  use dc_message
  use gtool_history
  use dc_calendar, only: &
         & DCCalCreate, DCCalDateCreate, &
         & DCCalConvertByUnit, DCCalDateDifference, &
         & DCCalDateInquire, DCCalDateEval

  !* Dennou-OGCM
  !
  use Constants_mod
  use TemporalIntegSet_mod, only: &
       & TemporalIntegSet_Init, TemporalIntegSet_Final, &
       & RestartTime

  use GridSet_mod
  use SpmlUtil_mod
  use BoundCondSet_mod

  use GovernEqSet_mod, only: &
       & GovernEqSet_Init, GovernEqSet_Final, &
       & EOSType

  use EOSDriver_mod

  use Exp_WindDrivenCirculation_mod, only: &
       & Exp_Init => Exp_WindDrivenCirculation_Init, &
       & Exp_Final => Exp_WindDrivenCirculation_Final, &
       & Exp_SetInitCond => SetInitCondition

  !* DiagVar
  !

  use DiagVarSet_mod
  use DiagVarFileSet_mod
  use DiagVarEval_mod
  use BudgetAnalysis_mod
  use SpectralAnalysis_mod



  implicit none

  ! ���ʸ; Declaration statement
  !

  character(*), parameter :: PROGRAM_NAME = 'diagVar'
  character(*), parameter :: DEFAULT_CONFIG_NML = 'defaultConfigDiagVar.nml' 
  character(TOKEN), pointer :: ogcmOutputVarsName(:) => null()
  character(TOKEN), pointer :: diagVarsName(:) => null()
  character(TOKEN), pointer :: BudgetTypesName(:) => null()
  character(TOKEN), pointer :: SpectralTypesName(:) => null()

  real(DP) :: CurrentTimeSec, EndTimeSec, TimeIntSec

  type(gtool_historyauto_info) :: diagVar_gthsInfo
  type(gtool_historyauto_info) :: ogcm_gthsInfo


  
  ! �¹�ʸ; Executable statement
  !

  !
  call MessageNotify('M', PROGRAM_NAME, "Setup...")
  call setup()

  !

  do while(CurrentTimeSec <= EndTimeSec)
     call MessageNotify('M', PROGRAM_NAME, "%f [%c]..", &
          & d=(/ DCCalConvertByUnit(CurrentTimeSec, 'sec', ogcm_gthsInfo%intUnit) /), c1=trim(ogcm_gthsInfo%intUnit) )

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

    ! �⥸�塼�����; Use statement
    !
    use OptionParser_mod, only: &
         & OptionParser_Init, OptionParser_Final, &
         & OptionParser_GetInfo

    use VariableSet_mod, only: VariableSet_Init
 
#ifdef _OPENMP
    use omp_lib
#endif


    ! ���ʸ; Declaration statement
    !
    
    
    ! �ɽ��ѿ�
    ! Local variables
    !
    character(STRING) :: configNmlFile
    character(STRING) :: ogcmConfigNmlFile

    integer :: nThread
    integer :: varID

    ! �¹�ʸ; Executable statement
    !

    EndTimeSec = -1d0

    call OptionParser_Init()
    call OptionParser_GetInfo(configNmlFile,   & !(out)
         & defaultConfigNml=DEFAULT_CONFIG_NML ) 
    call OptionParser_Final()
    
    call MessageNotify('M', PROGRAM_NAME, &
         & "Read configure file ' %c ' ..", c1=trim(configNmlFile) )

    call readNml(configNmlFile, &  ! (in)
         & ogcmConfigNmlFile )     ! (out)
    call readOgcmNml(ogcmConfigNmlFile)

    call Constants_Init(ogcmConfigNmlFile)
    call TemporalIntegSet_Init(ogcmConfigNmlFile)
    call BoundCondSet_Init(ogcmConfigNmlFile)
    call GridSet_Init(ogcmConfigNmlFile)
    call GovernEqSet_Init(ogcmConfigNmlFile)
    call EOSDriver_Init(EOSType)

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
    call Exp_Init(ogcmConfigNmlFile)
    call Exp_SetInitCond()

    !
    call MessageNotify('M', PROGRAM_NAME, &
         & "Some modules is initialized..")

    call DiagVarSet_Init(diagVarsName)
    call DiagVarEval_Init()
    call DiagVarFileSet_Init(configNmlFile, diagVar_gthsInfo)
    call BudgetAnalysis_Init(diagVar_gthsInfo, BudgetTypesName)
    call SpectralAnalysis_Init(diagVar_gthsInfo, SpectralTypesName)

     !
     call MessageNotify('M', PROGRAM_NAME, &
          & "Period %f - % f [%c]..", &
          & d=(/ DCCalConvertByUnit(CurrentTimeSec, 'sec', ogcm_gthsInfo%intUnit),    &
          &      DCCalConvertByUnit(EndTimeSec, 'sec', ogcm_gthsInfo%intUnit) /),     &
          & c1=trim(ogcm_gthsInfo%intUnit) )

  end subroutine setup

  !> @brief 
  !!
  !!
  subroutine shutdown()

    ! �⥸�塼�����
    ! Use statements
    use VariableSet_mod, only: VariableSet_Final

    ! ���ʸ; Declaration statement
    !
    
    
    ! �ɽ��ѿ�
    ! Local variables
    !
    
    
    ! �¹�ʸ; Executable statement
    !

    call DiagVarEval_Final()
    call DiagVarFileSet_Final()
    call BudgetAnalysis_Final()
    call SpectralAnalysis_Final()

    call DiagVarSet_Final()

    if( associated(diagVarsName) ) deallocate(diagVarsName)
    if( associated(ogcmOutputVarsName) ) deallocate(ogcmOutputVarsName)
    if( associated(BudgetTypesName) ) deallocate(BudgetTypesName)

    call EOSDriver_Final()
    call GovernEqSet_Final()
    call VariableSet_Final()
    call SpmlUtil_Final()
    call GridSet_Final()
    call BoundCondSet_Final()
    call TemporalIntegSet_Final()
    call Constants_Final()

  end subroutine shutdown

  !> @brief 
  !!
  !!
  subroutine readNml(configNml, ogcmConfigNml)

    ! �⥸�塼�����; Use statement
    !

    ! �ե��������������
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! ���̷��ѥ�᥿
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! ɸ����Ϥ������ֹ�. Unit number of standard output
    
    use dc_string, only: Split, Replace

    ! ���ʸ; Declaration statement
    !    
    character(*), intent(in) :: configNml
    character(STRING), intent(out) :: ogcmConfigNml

    ! �ɽ��ѿ�
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST �ե����륪���ץ��������ֹ�. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST �ɤ߹��߻��� IOSTAT. 
    ! IOSTAT of NAMELIST read

    character(TOKEN) :: pos_nml

    real(DP) :: StartTime, EndTime, IntValue
    character(TOKEN) :: TimeUnit, IntUnit
    character(STRING) :: Name, FilePrefix

    character(STRING) :: BudgetTypes
    character(STRING) :: SpectralTypes

    ! NAMELIST �ѿ���
    ! NAMELIST group name
    !
    namelist /diagVar_nml/ &
         & ogcmConfigNml, StartTime, EndTime, TimeUnit

    namelist /gtool_historyauto_nml/ &
         & IntValue, IntUnit, Name, FilePrefix
    
    namelist /budgetAnalysis_nml/ &
         & BudgetTypes

    namelist /spectralAnalysis_nml/ &
         & SpectralTypes

    ! �¹�ʸ; Executable statements

    ! �ǥե�����ͤ�����
    ! Default values settings
    !
    
    ogcmConfigNml = 'config.nml'
    EndTime = -1d0
    StartTime = 0d0
    TimeUnit = 'day'

    Name = ''
    FilePrefix = ''
    BudgetTypes = ''
    SpectralTypes = ''
    

    ! NAMELIST ���������
    ! Input from NAMELIST
    !
    if ( trim(configNml) /= '' ) then
       call MessageNotify( 'M',PROGRAM_NAME, "reading namelist '%a'", ca=(/ configNml /))
       call FileOpen( unit_nml, &          ! (out)
            & configNml, mode = 'r' ) ! (in)

       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = diagVar_nml, iostat = iostat_nml )   ! (out)

       !
       pos_nml = ''
       do while ( trim(pos_nml) /= 'APPEND' .and. iostat_nml == 0 ) 
          read( unit_nml, &           ! (in)
               & nml = gtool_historyauto_nml, iostat = iostat_nml )   ! (out)
          inquire( unit_nml, & !(in)
               & position=pos_nml ) !(out)
       end do
       
       !
       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = budgetAnalysis_nml, iostat = iostat_nml )   ! (out)

       !
       rewind( unit_nml )
       read( unit_nml, &           ! (in)
            & nml = spectralAnalysis_nml, iostat = iostat_nml )   ! (out)

       close( unit_nml )
    end if


    Name = Replace(Name, " ", "")
    call Split(trim(Name), diagVarsName, ",")

    diagVar_gthsInfo = gtool_historyauto_info( intValue=intValue, intUnit=intUnit, &
         &                                     FilePrefix=FilePrefix, Name=Name, origin=0d0 )

    BudgetTypes = Replace(BudgetTypes, " ", "")
    call Split(trim(BudgetTypes), BudgetTypesName, ",")

    SpectralTypes = Replace(SpectralTypes, " ", "")
    call Split(trim(SpectralTypes), SpectralTypesName, ",")


    CurrentTimeSec = DCCalConvertByUnit(StartTime, TimeUnit, "sec")
    TimeIntSec = DCCalConvertByUnit(diagVar_gthsInfo%intValue, diagVar_gthsInfo%intUnit, 'sec')
    if ( EndTime >= 0d0 ) EndTimeSec = DCCalConvertByUnit(EndTime, TimeUnit, "sec")

    ! ���� ; Print
    !
    call MessageNotify( 'M', PROGRAM_NAME, '----- Initialization Messages -----' )
    call MessageNotify( 'M', PROGRAM_NAME, ' ogcm configure file   = %c', c1=trim(ogcmConfigNml) )
    call MessageNotify( 'M', PROGRAM_NAME, ' diagnostic variables  = %c', c1=Name )
    call MessageNotify( 'M', PROGRAM_NAME, ' budget analysis types  = %c', c1=BudgetTypes )
    call MessageNotify( 'M', PROGRAM_NAME, ' spectral analysis types  = %c', c1=SpectralTypes )
    
  end subroutine readNml

  !> @brief 
  !!
  !!
  subroutine readOgcmNml(ogcmConfigNml)

    ! �⥸�塼�����; Use statement
    !

    ! �ե��������������
    ! File I/O support
    !
    use dc_iounit, only: FileOpen

    ! ���̷��ѥ�᥿
    ! Kind type parameter
    !
    use dc_types, only: STDOUT ! ɸ����Ϥ������ֹ�. Unit number of standard output
    
    use dc_string, only: Split, Replace

    ! ���ʸ; Declaration statement
    !    
    character(*), intent(in) :: ogcmConfigNml

    ! �ɽ��ѿ�
    ! Local variables
    !
    integer:: unit_nml        ! NAMELIST �ե����륪���ץ��������ֹ�. 
    ! Unit number for NAMELIST file open

    integer:: iostat_nml      ! NAMELIST �ɤ߹��߻��� IOSTAT. 
    ! IOSTAT of NAMELIST read

    character(TOKEN) :: pos_nml

    real(DP) :: IntValue
    character(TOKEN) :: IntUnit
    character(STRING) :: Name, FilePrefix

    real(DP), pointer :: ogcm_outputTime(:) => null()

    ! NAMELIST �ѿ���
    ! NAMELIST group name
    !
    namelist /gtool_historyauto_nml/ &
         & IntValue, IntUnit, Name, FilePrefix

    ! �¹�ʸ; Executable statements

    ! �ǥե�����ͤ�����
    ! Default values settings
    !
    FilePrefix = ''
    Name = ''

    ! NAMELIST ���������
    ! Input from NAMELIST
    !
    if ( trim(ogcmConfigNml) /= '' ) then
       call MessageNotify( 'M',PROGRAM_NAME, "reading ogcm namelist '%a'", ca=(/ ogcmConfigNml /))
       call FileOpen( unit_nml, &          ! (out)
            & ogcmConfigNml, mode = 'r' ) ! (in)

       pos_nml = ''; iostat_nml = 0
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
         &                                     FilePrefix=FilePrefix, Name=Name, origin=0d0 )


    if(EndTimeSec < 0d0 ) then
       call HistoryGetPointer( &
            & trim(ogcm_gthsInfo%FilePrefix) // trim(ogcmOutputVarsName(1)) // '.nc', &
            & 't', ogcm_outputTime)
       write(*,*) 'outputTime', size(ogcm_outputTime), ogcm_outputTime

       ogcm_gthsInfo%origin = ogcm_outputTime(1)
       diagVar_gthsInfo%origin = ogcm_outputTime(1)
       CurrentTimeSec = DCCalConvertByUnit(ogcm_outputTime(1), ogcm_gthsInfo%intUnit, "sec")
       EndTimeSec = DCCalConvertByUnit(ogcm_outputTime(size(ogcm_outputTime)), ogcm_gthsInfo%intUnit, "sec")
    end if

    ! ���� ; Print
    !
    call MessageNotify( 'M', PROGRAM_NAME, '----- Initialization Messages -----' )
    call MessageNotify( 'M', PROGRAM_NAME, ' ogcm output variables   = %c', c1=Name )

  end subroutine readOgcmNml


  !> @brief 
  !!
  !!
  subroutine diagnose_outputVariables()

    ! �⥸�塼�����; Use statement
    !
    use dc_string, only: CPrintf
    
    use VariableSet_mod, only: &
         & VARSET_KEY_U, VARSET_KEY_V, VARSET_KEY_SURFPRESS, VARSET_KEY_BAROCPRESS, VARSET_KEY_SURFHEIGHT, VARSET_KEY_PTEMPEDD, &
         & VARSET_KEY_UB, VARSET_KEY_VB, VARSET_KEY_PTEMPEDDB, &
         & VARSET_KEY_SIGDOT, VARSET_KEY_PTEMPBASIC, VARSET_KEY_TOTDEPTHBASIC, &
         & xyz_UN, xyz_VN, xy_SurfPressN, xy_SurfHeightN, xyz_PTempEddN, xyz_SaltN, &
         & xyz_UB, xyz_VB, xyz_PTempEddB, &
         & xyz_SigDot, z_PTempBasic, xy_totDepthBasic

    use DiagVarSet_mod, only: &
         & xyz_Div, xyz_Vor, xyz_BarocPress, xyz_PressEdd, &
         & xyz_DensEdd!, xy_totDepth

    ! ���ʸ; Declaration statement
    !
    
    ! �ɽ��ѿ�
    ! Local variables
    !
    integer :: varID
    character(STRING) :: rangeStr
    real(DP) :: CurrentTime
    real(DP) :: xy_totDepth(0:iMax-1,jMax), xyz_PTemp(0:iMax-1,jMax,0:kMax)
    character(*), parameter :: NCEXT =".nc"
    integer :: k

    ! �¹�ʸ; Executable statement
    !
    
    CurrentTime = DCCalConvertByUnit(CurrentTimeSec, 'sec', ogcm_gthsInfo%intUnit)
    rangeStr = CPrintf('t=%f', d=(/ currentTime /))

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_U) // NCEXT, &
         & VARSET_KEY_U, xyz_UN, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_V) // NCEXT, &
         & VARSET_KEY_V, xyz_VN, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_PTEMPEDD) // NCEXT, &
         & VARSET_KEY_PTEMPEDD, xyz_PTempEddN, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_UB) // NCEXT, &
         & VARSET_KEY_UB, xyz_UB, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_VB) // NCEXT, &
         & VARSET_KEY_VB, xyz_VB, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_PTEMPEDDB) // NCEXT, &
         & VARSET_KEY_PTEMPEDDB, xyz_PTempEddB, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_SIGDOT) // NCEXT, &
         & VARSET_KEY_SIGDOT, xyz_SigDot, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_SURFPRESS) // NCEXT, &
         & VARSET_KEY_SURFPRESS, xy_SurfPressN, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_BAROCPRESS) // NCEXT, &
         & VARSET_KEY_BAROCPRESS, xyz_BarocPress, range=rangeStr )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_PTEMPBASIC) // NCEXT, &
         & VARSET_KEY_PTEMPBASIC, z_PTempBasic )

    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // trim(VARSET_KEY_TOTDEPTHBASIC) // NCEXT, &
         & VARSET_KEY_TOTDEPTHBASIC, xy_totDepthBasic )

!!$    call HistoryGet( trim(ogcm_gthsInfo%FilePrefix) // "SurfHeight.nc", &
!!$         & 'SurfHeight', xy_SurfHeightN, range=rangeStr )

    xyz_SaltN = 0d0
    xy_totDepth = xy_totDepthBasic  !+ xy_SurfHeightN

    forAll(k=0:kMax) xyz_PTemp(:,:,k) = z_PTempBasic(k) + xyz_PTempEddN(:,:,k)

    xyz_PressEdd = eval_PressEdd(xy_SurfPressN, xyz_BarocPress)
    xyz_DensEdd = eval_DensEdd(xyz_PTemp, xyz_SaltN, xy_totDepth )

    do varID=1, size(diagVarsName)
       select case(diagVarsName(varID))
          case(DVARKEY_DIV)
             xyz_Div = eval_Div(xyz_UN, xyz_VN)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_Div, var3D=xyz_Div)
          case(DVARKEY_VOR)
             xyz_Vor = eval_Vor(xyz_UN, xyz_VN)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_Vor, var3D=xyz_Vor)
          case (DVARKEY_PRESSEDD)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_PRESSEDD, var3D=xyz_PressEdd )
          case (DVARKEY_MASSSTREAMFUNC)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_MASSSTREAMFUNC, &
                  & var2D=eval_MassStreamFunc(xyz_VN, xy_totDepth) )
          case (DVARKEY_STATICSTABILITY)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_STATICSTABILITY, &
                  & var3D=eval_StaticStability(xyz_PTemp, xy_totDepth) )
          case (DVARKEY_PTEMP)
             call DiagVarFileSet_OutputVar(CurrentTimeSec, DVARKEY_PTEMP, &
                  & var3D=xyz_PTemp )                
       end select
    end do

    ! Call some subroutines associated with budget and spectral analaysis.

    call SpectralAnalysis_Perform()
    call BudgetAnalysis_perform()
    
  end subroutine diagnose_outputVariables

end program diagVar_main
