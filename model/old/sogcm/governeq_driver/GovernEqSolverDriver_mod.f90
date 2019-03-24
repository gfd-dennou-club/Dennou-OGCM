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
  use DataFileSet_mod, only: &
       & DataFileSet

  use GovernEqSet_mod, only: &
       & DynEqType, EOSType, &
       & GOVERNEQSET_DYN_HYDROBOUSSINESQ, &
       & GOVERNEQSET_DYN_NONDYN_MIXEDLYR

  use HydroBoudEq_TimeInteg_v2_mod, only: &
       & HydroBouEq_TimeInteg_Init, HydroBouEq_TimeInteg_Final, &
       & HydroBouEqSolver_AdvanceTStep

  use NonDynMixedLyr_TimeInteg_mod, only: &
       & NonDynMixedLyr_TimeInteg_Init, NonDynMixedLyr_TimeInteg_Final, &
       & NonDynMixedLyr_AdvanceTStep
  
  use SeaiceEq_TimeInteg_mod, only: &
       & SeaiceEq_TimeInteg_Init, SeaiceEq_TimeInteg_Final, &
       & SeaiceEqSolver_AdvanceTStep
  
  use EOSDriver_mod, only: &
       & EOSDriver_Init, EOSDriver_Final



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

    ! 宣言文; Declaration statement
    !


    ! 実行文; Executable statements
    !

    call EOSDriver_Init(EOSType)

!!$    call HydroBouEqSolver_Init()
    select case (DynEqType)
    case (GOVERNEQSET_DYN_HYDROBOUSSINESQ)
       call HydroBouEq_TimeInteg_Init()
    case (GOVERNEQSET_DYN_NONDYN_MIXEDLYR)
       call NonDynMixedLyr_TimeInteg_Init()
    end select
          

  end subroutine GovernEqSolverDriver_Init

  !>
  !!
  !!
  subroutine GovernEqSolverDriver_Final()

    ! 実行文; Executable statements
    !

    select case (DynEqType)
    case (GOVERNEQSET_DYN_HYDROBOUSSINESQ)
       call HydroBouEq_TimeInteg_Final()
    case (GOVERNEQSET_DYN_NONDYN_MIXEDLYR)
       call NonDynMixedLyr_TimeInteg_Final()
    end select

    call EOSDriver_Final()

  end subroutine GovernEqSolverDriver_Final


  !> @brief 
  !!
  !!
  subroutine GovernEqSolverDriver_AdvanceTStep(isSelfStartSchemeUsed)
    
    ! モジュール引用; Use statement
    !
    
    use TemporalIntegUtil_mod, only: &
         & timeIntMode_RK4, timeIntMode_Euler

    use TemporalIntegSet_mod, only: &
         & CurrentTimeStep, DelTime, &
         & barocTimeIntMode, nStage_BarocTimeInt, isVarBUsed_BarocTimeInt

    use VarSetSeaice_mod, only: &
         & xy_SIceConA, xy_SIceConN, xy_Wice

    use BoundaryCondO_mod, only: &
         & BoundaryCondO_Update
    
    ! 宣言文; Declaration statement
    !
    logical, optional :: isSelfStartSchemeUsed
    
    ! 局所変数
    ! Local variables
    !
    integer :: timeIntMode
    integer :: nStage
    logical :: is_VarB_Used

    ! 実行文; Executable statement
    !

    timeIntMode = barocTimeIntMode
    nStage = nStage_BarocTimeInt
    is_VarB_Used = isVarBUsed_BarocTimeInt

    if(present(isSelfStartSchemeUsed) .and. isSelfStartSchemeUsed) then
       timeIntMode = timeIntMode_RK4; nStage=4; is_VarB_Used = .false.
    end if

    !
    call SeaiceEqSolver_AdvanceTStep( &
         & DelTime )

    !
    call BoundaryCondO_Update(xy_SIceConN, xy_Wice)

    select case (DynEqType)
    case (GOVERNEQSET_DYN_HYDROBOUSSINESQ)
       call HydroBouEqSolver_AdvanceTStep( &
            & DelTime, timeIntMode, nStage, is_VarB_Used )
    case (GOVERNEQSET_DYN_NONDYN_MIXEDLYR)
       call NonDynMixedLyr_AdvanceTStep( &
            & DelTime, timeIntMode, nStage, is_VarB_Used )
    end select
    

  end subroutine GovernEqSolverDriver_AdvanceTStep


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module GovernEqSolverDriver_mod

