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
  use PolyMesh_mod, only: &
       & PolyMesh
  
  use VariableSet_mod, only: &
       & VariableSet

  use HydroBouEqSolver_mod, only: &
       & HydroBouEqSolver_Init, HydroBouEqSolver_Final, &
       & HydroBouEqSolver_AdvanceTime

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: GovernEqSolverDriver_Init, GovernEqSolverDriver_Final
  public :: GovernEqSolverDriver_AdvanceTime

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


  end subroutine GovernEqSolverDriver_Init

  !> @brief 
  !!
  !!
  subroutine GovernEqSolverDriver_AdvanceTime(variable, mesh)
    
    ! 宣言文; Declaration statement
    !
    type(VariableSet), intent(inout) :: variable
    type(PolyMesh), intent(in) :: mesh
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    call HydroBouEqSolver_Init(mesh)
    call HydroBouEqSolver_AdvanceTime(variable)
    call HydroBouEqSolver_Final()

  end subroutine GovernEqSolverDriver_AdvanceTime

  !>
  !!
  !!
  subroutine GovernEqSolverDriver_Final()

    ! 実行文; Executable statements
    !

  end subroutine GovernEqSolverDriver_Final


end module GovernEqSolverDriver_mod

