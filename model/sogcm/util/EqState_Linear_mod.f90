!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module EqState_Linear_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: EqState_Linear_Init, EqState_Linear_Final
  public :: EqState_Linear_Eval

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'EqState_Linear_mod' !< Module Name
  real(DP), save :: RefDens = 1.023d03
  real(DP), save :: Cp0 = 3986d0
  real(DP), save :: BetaT = 1.67d-04
  real(DP), save :: T0 = 283d0

contains

  !>
  !!
  !!
  subroutine EqState_Linear_Init(refDens_, T0_, BetaT_, Cp0_)
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in), optional :: refDens_
    real(DP), intent(in), optional :: T0_
    real(DP), intent(in), optional :: BetaT_
    real(DP), intent(in), optional :: Cp0_

    ! 実行文; Executable statements
    !

    if(present(refDens_)) refDens = refDens_
    if(present(T0_)) T0 = T0_
    if(present(BetaT_)) BetaT = BetaT_
    if(present(Cp0_)) Cp0 = Cp0_

    call MessageNotify('M', module_name, &
         & "Set RefDens=%f, RefTemp=%f, ThermalExpanCoef=%f, Cp0=%f", d=(/refDens,T0,BetaT,Cp0/) )

  end subroutine EqState_Linear_Init

  !>
  !!
  !!
  subroutine EqState_Linear_Final()

    ! 実行文; Executable statements
    !

  end subroutine EqState_Linear_Final

  !> @brief 
  !!
  !!
  elemental function EqState_Linear_Eval(theta, s, p) result(densEdd)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: theta, s, p
    real(DP) :: densEdd

    ! 局所変数
    ! Local variables
    !

    !> Secan bulk modulus($K(S, \theta, p)$) 
    real(DP) :: K

    ! 実行文; Executable statement
    !
    densEdd = - refDens*( & 
         & BetaT*( theta*exp(-BetaT*p/(refDens*Cp0)) - T0 ) &
         & )

  end function EqState_Linear_Eval

end module EqState_Linear_mod
