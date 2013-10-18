program ogcm_main

  ! �⥸�塼�����; Use statement
  !
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use Constants_mod
  use TemporalIntegSet_mod
  use GridSet_mod
  use VariableSet_mod
  use DataFileSet_mod
  use GovernEqSolverDriver_mod

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

  call InitCond_Set()
  call DataFileSet_OutputData(datFile)
  
  call InitCond_Final()
  
  !***********************************
  ! The loop for temporal integration
  !************************************
  

  call MessageNotify("M", PROGRAM_NAME, "[==== Start temporal integration ====]")

  do while(CurrentTime < TotalIntegTime)
     
     !
     
     call GovernEqSolverDriver_AdvanceTStep()

     

     !
     call TemporalIntegSet_AdvanceLongTStep()
     call VariableSet_AdvanceTStep()
     
     !
     !
     
     call DataFileSet_OutputData(datFile)
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
         & nMax, tMax

    use SpmlUtil_mod, only: &
         SpmlUtil_Init

    use GridSet_mod, only: &
         & GridSet_construct

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
    call GridSet_Init(configNmlFile)

#ifdef _OPENMP
    !$omp parallel 
    !$omp single
    nThread = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, np=nThread)
#else
 
#endif

    call GridSet_construct()

    call VariableSet_Init()
    call DataFileSet_Init(datFile, configNmlFile)

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

    call DataFileSet_Final(datFile)
    call VariableSet_Final()
    call SpmlUtil_Final()
    call GridSet_Final() 
    call TemporalIntegSet_Final()
    call Constants_Final()

  end subroutine ogcm_finalize

end program ogcm_main
