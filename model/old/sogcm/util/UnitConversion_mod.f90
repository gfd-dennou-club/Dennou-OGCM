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
  public :: bar2Pa, Pa2bar, bar2dbar, degC2K, K2degC

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'UnitConversion_mod' !< Module Name

contains

elemental function  bar2Pa(bar) result(Pa)
  real(DP), intent(in) :: bar
  real(DP) :: Pa

  Pa = 1d05 * bar
end function bar2Pa

elemental function  Pa2bar(Pa) result(bar)
  real(DP), intent(in) :: Pa
  real(DP) :: bar

  bar = 1d-05 * Pa
end function Pa2bar

elemental function  bar2dbar(bar) result(dbar)
  real(DP), intent(in) :: bar
  real(DP) :: dbar

  dbar = 10d0 * bar

end function bar2dbar

elemental function  degC2K(degC) result(K)
  real(DP), intent(in) :: degC
  real(DP) :: K

  K = degC + 273.15d0

end function degC2K

elemental function  K2degC(K) result(degC)
  real(DP), intent(in) :: K
  real(DP) :: degC

  degC = K - 273.15d0

end function K2degC

end module UnitConversion_mod

