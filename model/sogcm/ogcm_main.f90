program ogcm_main

  ! �⥸�塼�����; Use statement
  !
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use Constants_mod
  use TemporalIntegSet_mod
  use BoundCondSet_mod
  use GridSet_mod
  use VariableSet_mod
  use DataFileSet_mod
  use RestartDataFileSet_mod
  use GovernEqSet_mod

  use GovernEqSolverDriver_mod

!!$  use Exp_W94_Case2_mod, only: &
!!$  use Exp_BarotRossbyWave_mod, only: &
!!$       & Exp_Init => Exp_BarotRossbyWave_Init, &
!!$       & Exp_Final => Exp_BarotRossbyWave_Final, &
!!$       & Exp_SetInitCond => SetInitCondition 
!!$  use Exp_InternalGravWave_mod, only: &
!!$       & Exp_Init => Exp_InternalGravWave_Init, &
!!$       & Exp_Final => Exp_InternalGravWave_Final, &
!!$       & Exp_SetInitCond => SetInitCondition
  
  use Exp_WindDrivenCirculation_mod, only: &
       & Exp_Init => Exp_WindDrivenCirculation_Init, &
       & Exp_Final => Exp_WindDrivenCirculation_Final, &
       & Exp_SetInitCond => SetInitCondition


  use InitCond_mod

  ! ���ʸ; Declaration statement
  !
  implicit none


  ! �ɽ��ѿ�
  ! Local variables
  !

  character(*), parameter :: PROGRAM_NAME = "ogcm_main"
  type(DataFileSet) :: datFile

  ! �¹�ʸ; Executable statement
  !

  call MessageNotify("M", PROGRAM_NAME, "Start..")

  !***********************************
  ! Set up
  !***********************************

  call ogcm_setup()

  !***********************************
  ! Set initial condition
  !***********************************

  call InitCond_Init()

  call InitCond_Set(Exp_SetInitCond)
  call RestartDataFileSet_Input()

  call DataFileSet_OutputBasicData()
  call DataFileSet_OutputData(datFile)

  call InitCond_Final()

  !***********************************
  ! The loop for temporal integration
  !************************************


  call MessageNotify("M", PROGRAM_NAME, "[==== Start temporal integration ====]")

  do while( .not. EndTemporalInteg() )

     !

     call GovernEqSolverDriver_AdvanceTStep()



     !
     call TemporalIntegSet_AdvanceLongTStep()
     call VariableSet_AdvanceTStep()

     !
     !

     call DataFileSet_OutputData(datFile)
     call RestartDataFileSet_Output()

  end do

  call MessageNotify("M", PROGRAM_NAME, "[==== Finish temporal integration ====]")

  !*************************************
  ! Finalize
  !*************************************

  call ogcm_finalize()

  call MessageNotify("M", PROGRAM_NAME, "..End")


contains

  !> @brief
  !!
  !!
  subroutine ogcm_setup()

    ! �⥸�塼�����; Use statements
    !
    use OptionParser_mod

    use GridSet_mod, only: &
         & iMax, jMax, kMax, &
         & nMax, tMax, &
         & GridSet_construct

    use SpmlUtil_mod, only: &
         SpmlUtil_Init

#ifdef _OPENMP
    use omp_lib
#endif

    ! ���ʸ; Declaration statement
    !


    ! �ɽ��ѿ�
    ! Local variables
    !
    character(STRING) :: configNmlFile
    integer :: nThread

    ! �¹�ʸ; Executable statement
    !

    call OptionParser_Init()
    call OptionParser_GetInfo(configNmlFile)
    call OptionParser_Final()

    call Constants_Init(configNmlFile)

    call TemporalIntegSet_Init(configNmlFile)
    call BoundCondSet_Init(configNmlFile)
    call GridSet_Init(configNmlFile)

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
    call DataFileSet_Init(datFile, configNmlFile)
    call RestartDataFileSet_Init(configNmlFile)
    call GovernEqSet_Init(configNmlFile)

    call GovernEqSolverDriver_Init()

    ! 
    call Exp_Init(configNmlFile)

  end subroutine ogcm_setup

  !> @brief
  !!
  !!
  subroutine ogcm_finalize()

    use SpmlUtil_mod, only: &
         SpmlUtil_Final

    ! ���ʸ; Declaration statement
    !


    ! �ɽ��ѿ�
    ! Local variables
    !


    ! �¹�ʸ; Executable statement
    !

    call Exp_Final()
    call GovernEqSolverDriver_Final()
    call DataFileSet_Final(datFile)
    call RestartDataFileSet_Final()
    call VariableSet_Final()
    call SpmlUtil_Final()
    call GridSet_Final()
    call BoundCondSet_Final()
    call TemporalIntegSet_Final()
    call Constants_Final()

  end subroutine ogcm_finalize

end program ogcm_main
