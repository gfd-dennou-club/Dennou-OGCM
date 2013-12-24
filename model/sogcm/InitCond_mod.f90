!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module InitCond_mod 

  ! モジュール引用; Use statements
  !

  use VariableSet_mod, only: &
       & xyz_UN, xyz_VN, xyz_SigDot, &
       & xyz_PTempEddN, &
       & xy_WindStressU, xy_WindStressV

!!$  use Exp_W94_Case2_mod, only: &
!!$  use Exp_BarotRossbyWave_mod, only: &
!!$       & SetInitCondition 
  use Exp_InternalGravWave_mod, only: &
       & SetInitCondition
!!$  use Exp_WindDrivenCirculation_mod, only: &
!!$       & SetInitCondition

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: InitCond_Init, InitCond_Final
  public :: InitCond_Set

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'InitCond_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine InitCond_Init()

    ! 実行文; Executable statements
    !
    xyz_UN = 0d0; xyz_VN = 0d0; xyz_SigDot = 0d0
    xyz_PTempEddN = 0d0
    xy_WindStressU = 0d0; xy_WindStressV = 0d0;

  end subroutine InitCond_Init

  !> @brief 
  !!
  !!
  subroutine InitCond_Set()
    
    !
    !

    ! 宣言文; Declaration statement
    !
    
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !

    

    call SetInitCondition()

  end subroutine InitCond_Set

  !>
  !!
  !!
  subroutine InitCond_Final()

    ! 実行文; Executable statements
    !

  end subroutine InitCond_Final

end module InitCond_mod

