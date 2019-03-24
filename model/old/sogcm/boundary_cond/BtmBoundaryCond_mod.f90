!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module BtmBoundaryCond_mod 

  ! モジュール引用; Use statements
  !

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: BtmBoundaryCond_Init, BtmBoundaryCond_Final

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'BtmBoundaryCond_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine BtmBoundaryCond_Init()

    ! 実行文; Executable statements
    !

  end subroutine BtmBoundaryCond_Init

  !>
  !!
  !!
  subroutine BtmBoundaryCond_Final()

    ! 実行文; Executable statements
    !

  end subroutine BtmBoundaryCond_Final

end module BtmBoundaryCond_mod

