!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DGCoordConvert_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP
  use dc_message, only: &
       & MessageNotify

  use VectorSpace_mod
  use PolyMesh_mod
  use fvMeshInfo_mod
  use HexTriIcMesh_mod

  use SimParameters_mod, only: &
       & Radius

  use GridSet_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DGCoordConvert_Init, DGCoordConvert_Final

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DGCoordConvert_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DGCoordConvert_Init()

    ! 実行文; Executable statements
    !

  end subroutine DGCoordConvert_Init

  !>
  !!
  !!
  subroutine DGCoordConvert_Final()

    ! 実行文; Executable statements
    !

  end subroutine DGCoordConvert_Final

!!!!!!!!!!!!!!!!!!!


end module DGCoordConvert_mod

