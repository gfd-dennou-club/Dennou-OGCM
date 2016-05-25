!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module LinearMapping_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP
  use VectorSpace_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: rotateZ

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'LinearMapping_mod' !< Module Name

contains

function rotateZ(vec, angle) result(rotVec)
  type(Vector3d), intent(in) :: vec
  real(DP), intent(in) :: angle
  type(Vector3d) :: rotVec

  rotVec%v_(1) = cos(angle)*vec%v_(1) - sin(angle)*vec%v_(2)
  rotVec%v_(2) = sin(angle)*vec%v_(1) + cos(angle)*vec%v_(2)
  rotVec%v_(3) = vec%v_(3)
end function rotateZ

end module LinearMapping_mod

