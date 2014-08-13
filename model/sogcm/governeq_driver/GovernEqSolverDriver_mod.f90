!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module GovernEqSolverDriver_mod 

  ! モジュール引用; Use statements
  !
  use GovernEqSet_mod, only: &
       & EOSType

!!$  use HydroBouEqSolver_mod, only: &
!!$       & HydroBouEqSolver_Init, HydroBouEqSolver_Final, &
!!$       & HydroBouEqSolver_AdvanceTStep

  use HydroBoudEq_TimeInteg_mod, only: &
       & HydroBouEq_TimeInteg_Init, HydroBouEq_TimeInteg_Final, &
       & HydroBouEqSolver_AdvanceTStep

  use EOSDriver_mod, only: &
       & EOSDriver_Init, EOSDriver_Final

!!$  use HydroBouEqSolverSelfStart_mod, only: &
!!$       & HydroBouEqSolverSelfStart_Init, HydroBouEqSolverSelfStart_Final, &
!!$       & HydroBouEqSolverSelfStart_AdvanceTStep

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: GovernEqSolverDriver_Init, GovernEqSolverDriver_Final
  public :: GovernEqSolverDriver_AdvanceTStep

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'GovernEqSolverDriver_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine GovernEqSolverDriver_Init()

    ! 実行文; Executable statements
    !

    call EOSDriver_Init(EOSType)

!!$    call HydroBouEqSolver_Init()
    call HydroBouEq_TimeInteg_Init()

  end subroutine GovernEqSolverDriver_Init

  !> @brief 
  !!
  !!
  subroutine GovernEqSolverDriver_AdvanceTStep()
    
    ! モジュール引用; Use statement
    !
    use TemporalIntegUtil_mod, only: &
         & timeIntMode_RK4, timeIntMode_Euler

    use TemporalIntegSet_mod, only: &
         & CurrentTimeStep, DelTime, &
         & barocTimeIntMode, nStage_BarocTimeInt, isVarBUsed_BarocTimeInt, &
         & SemiImplicitFlag

    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
    integer :: timeIntMode
    integer :: nStage
    logical :: is_VarB_Used
    logical :: is_VStiffTerm_Implicit

    ! 実行文; Executable statement
    !

    timeIntMode = barocTimeIntMode
    nStage = nStage_BarocTimeInt
    is_VarB_Used = isVarBUsed_BarocTimeInt
    is_VStiffTerm_Implicit = SemiImplicitFlag

    if(CurrentTimeStep == 1) then
       ! For first time step, RK4 which  has an ability to self-start is used. 
       timeIntMode = timeIntMode_RK4; nStage=4; is_VarB_Used = .false.
       is_VStiffTerm_Implicit = .false.
    end if

    call HydroBouEqSolver_AdvanceTStep( &
         & DelTime, timeIntMode, nStage, is_VarB_Used )

  end subroutine GovernEqSolverDriver_AdvanceTStep

  !>
  !!
  !!
  subroutine GovernEqSolverDriver_Final()

    ! 実行文; Executable statements
    !

!!$    call HydroBouEqSolver_Final()
    call HydroBouEq_TimeInteg_Final()
    call EOSDriver_Final()

  end subroutine GovernEqSolverDriver_Final


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GovernEqSolverDriver_mod

