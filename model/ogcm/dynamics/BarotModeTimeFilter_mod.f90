!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module BarotModeTimeFilter_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: BarotModeTimeFilter_Init, BarotModeTimeFilter_Final

  real(DP), public, save, allocatable :: Am(:), Bm(:)

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'BarotModeTimeFilter_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine BarotModeTimeFilter_Init()

    ! 実行文; Executable statements
    !

  end subroutine BarotModeTimeFilter_Init

  !>
  !!
  !!
  subroutine BarotModeTimeFilter_Final()

    ! 実行文; Executable statements
    !

  end subroutine BarotModeTimeFilter_Final

end module BarotModeTimeFilter_mod
