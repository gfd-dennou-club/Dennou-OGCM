!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module EqState_JM95_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  
  public :: EqState_JM95_Init, EqState_JM95_Final
  public :: EqState_JM95_Eval

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'EqState_JM95_mod' !< Module Name
  
  ! K(theta,S,p)
  real(DP), parameter :: K_0_0_0  =  1.965933d04
  real(DP), parameter :: K_1_0_0  =  1.444304d02
  real(DP), parameter :: K_2_0_0  = -1.706103d00
  real(DP), parameter :: K_3_0_0  =  9.648704d-03 
  real(DP), parameter :: K_4_0_0  = -4.190253d-05
  real(DP), parameter :: K_0_1_0  =  5.284855d01
  real(DP), parameter :: K_1_1_0  = -3.101089d-01
  real(DP), parameter :: K_2_1_0  =  6.283263d-03
  real(DP), parameter :: K_3_1_0  = -5.084188d-05
  real(DP), parameter :: K_0_15_0 =  3.886640d-01
  real(DP), parameter :: K_1_15_0 =  9.085835d-03
  real(DP), parameter :: K_2_15_0 = -4.619924d-04
  real(DP), parameter :: K_0_0_1  =  3.186519d00  
  real(DP), parameter :: K_1_0_1  =  2.212276d-02
  real(DP), parameter :: K_2_0_1  = -2.984642d-04
  real(DP), parameter :: K_3_0_1  =  1.956415d-06
  real(DP), parameter :: K_0_1_1  =  6.704388d-03
  real(DP), parameter :: K_1_1_1  = -1.847318d-04
  real(DP), parameter :: K_2_1_1  =  2.059331d-07
  real(DP), parameter :: K_0_15_1 =  1.480266d-04
  real(DP), parameter :: K_0_0_2  =  2.102898d-04
  real(DP), parameter :: K_1_0_2  = -1.202016d-05
  real(DP), parameter :: K_2_0_2  =  1.394680d-07
  real(DP), parameter :: K_0_1_2  = -2.040237d-06
  real(DP), parameter :: K_1_1_2  =  6.128773d-08
  real(DP), parameter :: K_2_1_2  =  6.207323d-10


contains

  !>
  !!
  !!
  subroutine EqState_JM95_Init()

    ! 実行文; Executable statements
    !

  end subroutine EqState_JM95_Init

  !>
  !!
  !!
  subroutine EqState_JM95_Final()

    ! 実行文; Executable statements
    !

  end subroutine EqState_JM95_Final

  !> @brief 
  !!
  !!
  elemental function EqState_JM95_Eval(theta, s, p) result(rho)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: theta, s, p
    real(DP) :: rho

    ! 局所変数
    ! Local variables
    !
    real(DP) :: rho0, rho00, A, B, C

    !> Secan bulk modulus($K(S, \theta, p)$) 
    real(DP) :: K

    ! 実行文; Executable statement
    !

    rho00 = 999.842594d0 + theta*( &
         &  6.793952d-02 + theta*(-9.095290d-03 + theta*(1.001685d-04 + theta*(-1.120083d-06 + theta*6.536332d-09))) &
         & )
    A = 8.24493d-01 + theta*(-4.0899d-03 + theta*(7.6438d-05 + theta*(-8.2467d-07 + theta*5.3875d-09)))
    B = -5.72466d-03 + theta*(1.0227d-04 - theta*1.6546d-06)
    C = 4.8314d-04
    rho0 = rho00 + S*(A + sqrt(S)*B + S*C)
 
    K = K_0_0_0 + &
         &   theta*(K_1_0_0 + theta*(K_2_0_0 + theta*(K_3_0_0 + theta*K_4_0_0))) &
         & + S*(K_0_1_0 + theta*(K_1_1_0 + theta*(K_2_1_0 + theta*K_3_1_0))) &
         & + S*sqrt(S)*(K_0_15_0 + theta*(K_1_15_0 + theta*K_2_15_0)) &
         & + p*(  K_0_0_1 + theta*(K_1_0_1 + theta*(K_2_0_1 + theta*K_3_0_1 + S*K_2_1_1) + S*K_1_1_1) &
         &      + S*( K_0_1_1 + sqrt(S)*K_0_15_1) ) &
         & + p**2*(K_0_0_2 + theta*(K_1_0_2 + theta*(K_2_0_2 + S*K_2_1_2) + S*K_1_1_2) + S*K_0_1_2)
 
    rho = rho0 / ( 1d0 - p/K )

  end function EqState_JM95_Eval

end module EqState_JM95_mod

