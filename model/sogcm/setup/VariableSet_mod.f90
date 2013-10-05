!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module VariableSet_mod 

  ! モジュール引用; Use statements
  !

  use dc_types, only: &
       & DP, STRING 

  use TemporalIntegSet_mod, only: &
       & nLongTimeLevel

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: VariableSet_Init, VariableSet_Final

  ! 公開変数
  ! Public variable
  !

  integer, public, save :: TracerNum
  integer, public, save :: SaltTracerID
  integer, public, save :: PTempTracerID
  real(DP), public, save :: refDens


  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'VariableSet_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine VariableSet_Init()

    ! モジュール引用; Use statements
    !
    use PolyMesh_mod, only: &
         & PolyMesh

    use GridSet_mod, only: &
         & kMax

    ! 宣言文; Declaration statement
    !

    ! 局所変数
    ! Local variable
    !
    integer :: n
    integer :: tl

    ! 実行文; Executable statements
    !


    !
    TracerNum = 2
    SaltTracerID = 1
    PTempTracerID = 2
    

  end subroutine VariableSet_Init

  !>
  !!
  !!
  subroutine VariableSet_Final()

    ! モジュール引用; Use statements
    !

    ! 宣言文; Declaration statement
    !


    ! 局所変数
    ! Local variable
    !
    integer :: n, tl
    
    ! 実行文; Executable statements
    !

  end subroutine VariableSet_Final

end module VariableSet_mod

