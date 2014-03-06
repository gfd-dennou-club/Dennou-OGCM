!-------------------------------------------------------------
! Copyright (c) 2013-2014 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module UnitConversion_mod 

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
  public :: bar2Pa, bar2dbar, degC2K

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'UnitConversion_mod' !< Module Name

contains

pure function  bar2Pa(bar) result(Pa)
  real(DP), intent(in) :: bar
  real(DP) :: Pa

  Pa = 1d05 * bar
end function bar2Pa

pure function  bar2dbar(bar) result(dbar)
  real(DP), intent(in) :: bar
  real(DP) :: dbar

  dbar = 10d0 * bar

end function bar2dbar

pure function  degC2K(degC) result(K)
  real(DP), intent(in) :: degC
  real(DP) :: K

  K = degC + 273.15d0

end function degC2K

end module UnitConversion_mod

