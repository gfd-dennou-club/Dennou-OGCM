!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module SGSVMixingConst_mod 

  ! モジュール引用; Use statements
  !

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SGSVMixingConst_Init, SGSVMixingConst_Final

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SGSVMixingConst_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine SGSVMixingConst_Init()

    ! 実行文; Executable statements
    !

  end subroutine SGSVMixingConst_Init

  !>
  !!
  !!
  subroutine SGSVMixingConst_Final()

    ! 実行文; Executable statements
    !

  end subroutine SGSVMixingConst_Final

  
end module SGSVMixingConst_mod

