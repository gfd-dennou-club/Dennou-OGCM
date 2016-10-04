!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module DOGCM_LPhys_DIFF_driver_mod 

  ! モジュール引用; Use statements
  !

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DOGCM_LPhys_DIFF_driver_Init, DOGCM_LPhys_DIFF_driver_Final

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_LPhys_DIFF_driver_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DOGCM_LPhys_DIFF_driver_Init()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_LPhys_DIFF_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_LPhys_DIFF_driver_Final()

    ! 実行文; Executable statements
    !

  end subroutine DOGCM_LPhys_DIFF_driver_Final

  !---------------------------------------------

  
end module DOGCM_LPhys_DIFF_driver_mod

