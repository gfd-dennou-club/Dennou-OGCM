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
  use HydroBouEqSolver_mod, only: &
       & HydroBouEqSolver_Init, HydroBouEqSolver_Final, &
       & HydroBouEqSolver_AdvanceTStep, &
       & timeIntMode_RK4, timeIntMode_LFTR, timeIntMode_LFAM3

  use EqState_JM95_mod, only: &
       & EqState_JM95_Init, EqState_JM95_Final

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
    call EqState_JM95_Init()

  end subroutine GovernEqSolverDriver_Init

  !> @brief 
  !!
  !!
  subroutine GovernEqSolverDriver_AdvanceTStep()
    
    ! モジュール引用; Use statement
    !
    use TemporalIntegSet_mod, only: &
         & CurrentTimeStep

    ! 宣言文; Declaration statement
    !
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    call HydroBouEqSolver_Init()

    if(CurrentTimeStep /= 1) then
       call HydroBouEqSolver_AdvanceTStep(timeIntMode_LFAM3)!TR)
    else
       ! For first time step, RK4 which is self starting is used. 
       call HydroBouEqSolver_AdvanceTStep(timeIntMode_RK4)
    end if

    call HydroBouEqSolver_Final()

  end subroutine GovernEqSolverDriver_AdvanceTStep

  !>
  !!
  !!
  subroutine GovernEqSolverDriver_Final()

    ! 実行文; Executable statements
    !

    call EqState_JM95_Final()

  end subroutine GovernEqSolverDriver_Final


end module GovernEqSolverDriver_mod

