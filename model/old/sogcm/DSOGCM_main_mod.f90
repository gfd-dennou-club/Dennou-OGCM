!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DSOGCM_main_mod 

  ! モジュール引用; Use statements
  !

  !+ gtool5

  use dc_types, only: &
       & DP, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  
  use Constants_mod
  use TemporalIntegSet_mod
  use BoundCondSet_mod

  use GridSet_mod
  use VariableSet_mod
  use BoundaryCondO_mod
  use GovernEqSet_mod
  use PhysicsDriver_mod

  use VarSetSeaice_mod  
  use SeaIceConstants_mod

  use GovernEqSolverDriver_mod
  
  use DataFileSet_mod
  use RestartDataFileSet_mod

!!$  use Exp_W94_Case2_mod, only: &
!------------------------------------------------
!!$  use Exp_BarotRossbyWave_mod, only: &
!!$       & Exp_Init => Exp_BarotRossbyWave_Init, &
!!$       & Exp_Final => Exp_BarotRossbyWave_Final, &
!!$       & Exp_SetInitCond => SetInitCondition 
!------------------------------------------------
!!$  use Exp_InternalGravWave_mod, only: &
!!$       & Exp_Init => Exp_InternalGravWave_Init, &
!!$       & Exp_Final => Exp_InternalGravWave_Final, &
!!$       & Exp_SetInitCond => SetInitCondition
!------------------------------------------------  
!!$  use Exp_WindDrivenCirculation_mod, only: &
!!$       & Exp_Init => Exp_WindDrivenCirculation_Init, &
!!$       & Exp_Final => Exp_WindDrivenCirculation_Final, &
!!$       & Exp_SetInitCond => SetInitCondition
!------------------------------------------------  
!!$  use Exp_APEOGCirc_mod, only: &
!!$       & Exp_Init => Exp_APEOGCirc_Init, &
!!$       & Exp_Final => Exp_APEOGCirc_Final, &
!!$       & Exp_SetInitCond => SetInitCondition
!------------------------------------------------  
!!$  use Exp_APEOGCircSeaice_mod, only: &
!!$       & Exp_Init => Exp_APEOGCircSeaice_Init, &
!!$       & Exp_Final => Exp_APEOGCircSeaice_Final, &
!!$       & Exp_SetInitCond => SetInitCondition

  use Exp_APECoupledAOGCMSeaice_mod, only: &
       & Exp_Init => Exp_APECoupledAOGCMSeaice_Init, &
       & Exp_Final => Exp_APECoupledAOGCMSeaice_Final, &
       & Exp_SetInitCond => SetInitCondition

  use InitCond_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSOGCM_main_Init, DSOGCM_main_Final
  public :: DSOGCM_main_setup, DSOGCM_main_shutdown
  public :: DSOGCM_main_advance_timestep
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !

  type(DataFileSet), save :: datFile
  logical :: CoupledRunFlag
  
  character(*), parameter:: module_name = 'DSOGCM_main_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DSOGCM_main_Init(isCoupledRun)

    logical, intent(in), optional :: isCoupledRun
    
    ! 実行文; Executable statements
    !

    CoupledRunFlag = .false.
    if(present(isCoupledRun)) CoupledRunFlag = .true.
    
  end subroutine DSOGCM_main_Init

  !>
  !!
  !!
  subroutine DSOGCM_main_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSOGCM_main_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> @brief 
  !!
  !!
  subroutine DSOGCM_main_advance_timestep( &
       & tstep, loop_end_flag  )
    
    ! 宣言文; Declaration statement
    !
    integer, intent(in) :: tstep
    logical, intent(inout) :: loop_end_flag
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    if( EndTemporalInteg() ) then
       loop_end_flag = .true.; return
    else
       loop_end_flag = .false.
    end if

    
    if(tstep == 0) then
       !***********************************
       ! Set initial condition
       !***********************************
       
       call InitCond_Init()

       call InitCond_Set(Exp_SetInitCond)
       call RestartDataFileSet_Input()

       call DataFileSet_OutputBasicData()
       call DataFileSet_OutputData(datFile)
       call PhysicsDriver_OutputData(datFile)

       call InitCond_Final()
       return
    end if

    
    ! Solve governing (partial diffrential) equation with line method. 
    ! Evaluate RHS of governing equation with some spatial discrete scheme, which  obtains the
    ! time tendency of prognostic variables. Then calculate values of prognostic variables at
    ! next time step using temporal scheme.
    !
     
    if(CurrentTimeStep == 1 .and. (.not. RestartFlag)) then
       call GovernEqSolverDriver_AdvanceTStep(isSelfStartSchemeUsed=.true.)
    else
       call GovernEqSolverDriver_AdvanceTStep()
    end if

    ! Update variables in order to advance next time step
    !
    call TemporalIntegSet_AdvanceLongTStep()
    call VariableSet_AdvanceTStep()
    call VarSetSeaice_AdvanceTStep()

    ! Output data of variables
    !
    call DataFileSet_OutputData(datFile)
    call PhysicsDriver_OutputData(datFile)
    call RestartDataFileSet_Output()
    
  end subroutine DSOGCM_main_advance_timestep

  !> @brief 
  !!
  !!
  subroutine DSOGCM_main_setup()


    ! モジュール引用; Use statements
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

    ! 宣言文; Declaration statement
    !


    ! 局所変数
    ! Local variables
    !
    character(STRING) :: configNmlFile
    integer :: nThread

    ! 実行文; Executable statement
    !

    !
    call OptionParser_Init()
    call OptionParser_GetInfo(configNmlFile)
    call OptionParser_Final()

    !
    call Constants_Init(configNmlFile)
    call SeaIceConstants_Init(configNmlFile)
    
    call TemporalIntegSet_Init(configNmlFile)

    !
    call GovernEqSet_Init(configNmlFile)    
    call BoundCondSet_Init(configNmlFile)
    call GridSet_Init(configNmlFile)

#ifdef _OPENMP
    !$omp parallel
    !$omp single
    nThread = omp_get_num_threads()
    !$omp end single
    !$omp end parallel

    call MessageNotify('M', module_name, "Execute as Thread Parallel Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet, np=nThread)
#else
    call MessageNotify('M', module_name, "Execute as Serial  Mode..")
    call SpmlUtil_Init(iMax, jMax, kMax, nMax, tMax, RPlanet)
#endif

    call GridSet_construct()
    call VariableSet_Init()
    call VarSetSeaice_Init(nIceLyr=2)    
    
    call DataFileSet_Init(datFile, configNmlFile)
    call RestartDataFileSet_Init(configNmlFile)

    call BoundaryCondO_Init()    
    call PhysicsDriver_Init(datFile, configNmlFile)
    call GovernEqSolverDriver_Init()
    
    call Exp_Init(configNmlFile)
    
  end subroutine DSOGCM_main_setup

  !> @brief 
  !!
  !!
  subroutine DSOGCM_main_shutdown()
    ! モジュール引用; Use statements
    !
    use SpmlUtil_mod, only: &
         SpmlUtil_Final

    ! 宣言文; Declaration statement
    !


    ! 局所変数
    ! Local variables
    !


    ! 実行文; Executable statement
    !

    call Exp_Final()
    call PhysicsDriver_Final()
    call GovernEqSolverDriver_Final()
    call BoundaryCondO_Final()
    call DataFileSet_Final(datFile)
    call RestartDataFileSet_Final()
    call VariableSet_Final()
    call VarSetSeaice_Final()
    call SpmlUtil_Final()
    call GridSet_Final()
    call BoundCondSet_Final()
    call TemporalIntegSet_Final()
    call Constants_Final()
    call SeaIceConstants_Final()
    
  end subroutine DSOGCM_main_shutdown

end module DSOGCM_main_mod

