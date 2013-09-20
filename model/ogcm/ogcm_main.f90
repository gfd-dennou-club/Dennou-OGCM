program ogcm_main

  ! モジュール引用; Use statement
  !
  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  use TemporalIntegSet_mod, only: &
       & CurrentTime, TotalIntegTime, &
       & TemporalIntegSet_AdvanceTime

  use DataFileSet_mod, only: &
       & DataFileSet_OutputData

  use GovernEqSolverDriver_mod, only: &
       & GovernEqSolverDriver_AdvanceTime

  ! 宣言文; Declaration statement
  !
  implicit none


  ! 局所変数
  ! Local variables
  !
  
  character(*), parameter :: PROGRAM_NAME = "ogcm_main"

  ! 実行文; Executable statement
  ! 

  call MessageNotify("M", PROGRAM_NAME, "Start..")

  !***********************************
  ! Set up
  !***********************************

  call ogcm_setup()

  !***********************************
  ! Set initial condition
  !***********************************
  
  call DataFileSet_OutputData()

  !***********************************
  ! The loop for temporal integration
  !************************************
  

  call MessageNotify("M", PROGRAM_NAME, "[==== Start temporal integration ====]")

  do while(CurrentTime < TotalIntegTime)
     
     !
     call GovernEqSolverDriver_AdvanceTime()
     
     !
     call TemporalIntegSet_AdvanceTime()

     !
     call DataFileSet_OutputData()
  end do

  call MessageNotify("M", PROGRAM_NAME, "[==== Finish temporal integration ====]")

  !*************************************
  ! Finalize 
  !*************************************

  call ogcm_finalize()

contains

  !> @brief 
  !!
  !!
  subroutine ogcm_setup()
    
    ! 引用文; Use statement
    !
    use OptionParser_mod, only: &
         & OptionParser_Init, OptionParser_Final, &
         & OptionParser_getInfo

    use Constants_mod, only: &
         & Constants_Init
    
    use GridSet_mod, only: &
         & GridSet_Init, fvmInfo

    use fvCalculus_mod, only: &
         & fvCalculus_Init

    use TemporalIntegSet_mod, only: &
         & TemporalIntegSet_Init

    use VariableSet_mod, only: &
         & VariableSet_Init

    use DataFileSet_mod, only: &
         & DataFileSet_Init

    use CGridFieldDataUtil_mod, only: &
         & CGridFieldDataUtil_Init

    use VGridFieldDataUtil_mod, only: &
         & VGridFieldDataUtil_Init

    use GovernEqSolverDriver_mod, only: &
         & GovernEqSolverDriver_Init

    ! 宣言文; Declaration statement
    !

    
    ! 局所変数
    ! Local variables
    !
    character(STRING) :: configNmlFile
    
    ! 実行文; Executable statement
    !

    call MessageNotify("M", PROGRAM_NAME, "[==== Set up some modules ====]")

    !
    call OptionParser_Init()
    call OptionParser_getInfo(configNmlFile)
    call OptionParser_Final()

    !
    call Constants_Init(configNmlFile)

    !
    call GridSet_Init(configNmlFile)

    !
    call fvCalculus_Init(fvmInfo)

    !
    call TemporalIntegSet_Init(configNmlFile)

    !
    call VariableSet_Init()

    !
    call DataFileSet_Init(configNmlFile)

    !
    call CGridFieldDataUtil_Init()

    !
    call VGridFieldDataUtil_Init()

    !
    call GovernEqSolverDriver_Init()

  end subroutine ogcm_setup

  !> @brief 
  !!
  !!
  subroutine ogcm_finalize()

    ! 引用文; Use statement
    !
    use Constants_mod, only: &
         & Constants_Final

    use GridSet_mod, only: &
         & GridSet_Final

    use fvCalculus_mod, only: &
         & fvCalculus_Final

    use TemporalIntegSet_mod, only: &
         & TemporalIntegSet_Final

    use VariableSet_mod, only: &
         & VariableSet_Final

    use DataFileSet_mod, only: &
         & DataFileSet_Final

    use CGridFieldDataUtil_mod, only: &
         & CGridFieldDataUtil_Final

    use VGridFieldDataUtil_mod, only: &
         & VGridFieldDataUtil_Final

    use GovernEqSolverDriver_mod, only: &
         & GovernEqSolverDriver_Final

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    call MessageNotify("M", PROGRAM_NAME, "[==== Finalize some modules ====]")

    !
    call GovernEqSolverDriver_Final()
    
    !
    call VGridFieldDataUtil_Final()

    !
    call CGridFieldDataUtil_Final()

    !
    call TemporalIntegSet_Final()

    !
    call DataFileSet_Final()

    !
    call VariableSet_Final()

    !
    call fvCalculus_Final()

    !
    call GridSet_Final()

    !
    call Constants_Final()


  end subroutine ogcm_finalize


end program ogcm_main
