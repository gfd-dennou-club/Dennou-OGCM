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
 
  real(DP) :: RefDens      !< Reference density
  real(DP) :: RefTemp      !< Reference temperature
  real(DP) :: PTempDepCoef
  real(DP) :: PressDepCoef

contains

  !>
  !!
  !!
  subroutine EqState_Linear_Init( &
       & refDens_, refTemp_, BetaT_, Cp0_, Cs0_  )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in), optional :: refDens_   !< Reference density
    real(DP), intent(in), optional :: refTemp_   !< Reference temperature
    real(DP), intent(in), optional :: BetaT_     !< (First) thermal expansion coffecient
    real(DP), intent(in), optional :: Cp0_       !< Specific heat capacity at constant pressure
    real(DP), intent(in), optional :: Cs0_       !< Reference sound speed

    ! 作業変数
    ! Work variables
    !
    real(DP) :: Cp0          !< Specific heat capacity at constant pressure
    real(DP) :: BetaT        !< (First) thermal expansion coffecient
    real(DP) :: Cs0          !< Reference sound speed

    ! 実行文; Executable statements
    !

    ! Set the default value of parameters used in linear EOS.  
    RefDens = 1.027d03  
    RefTemp = 283d0          
    Cp0 = 3986d0        
    BetaT = 1.67d-04    
    Cs0 = 1490d0  

    if(present(refDens_)) RefDens = refDens_
    if(present(refTemp_)) RefTemp = refTemp_
    if(present(BetaT_)) BetaT = BetaT_
    if(present(Cp0_)) Cp0 = Cp0_
    if(present(Cs0_)) Cs0 = Cs0_

    call MessageNotify('M', module_name, &
         & "Set RefDens=%f, RefTemp=%f, BetaT=%f, Cp0=%f, Cs0=%f", &
         & d=(/ RefDens, RefTemp, BetaT, Cp0, Cs0 /) )

    PTempDepCoef = - BetaT
    PressDepCoef = 1d0/(RefDens*Cs0**2) - RefTemp*BetaT**2/(RefDens*Cp0)


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
  elemental function EqState_Linear_Eval( theta, s, p ) result( densEdd )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: theta, s, p
    real(DP) :: densEdd

    ! 局所変数
    ! Local variables
    !

    ! 実行文; Executable statement
    !

!!$    densEdd = - refDens*( & 
!!$         & BetaT*( theta*exp(-BetaT*p/(refDens*Cp0)) - T0 ) &
!!$         & )

    densEdd = refDens*( & 
         &  PTempDepCoef*( theta - RefTemp )  + PressDepCoef*p &
         & )

  end function EqState_Linear_Eval

end module EqState_Linear_mod
