!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief This module provides a practical equation of state for seawater in Jackett and McDougall(1995). 
!! 
!! @author Kawai Yuta
!!
!! @detail
!> The EOS for seawater of JM95 is based on EOS-80 (UNESCO), 
!! but modified in order to use the potential temparature instead of in-situ 
!! temperature in calulating the density. 
!!
!! In this module, the units for temperature, salinity and pressure are 
!! [degree C], [psu] and [bar] respectively. 
!!
module EOS_JM95_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN

  use UnitConversion_mod, only: &
       & bar2dbar

  use SeaWaterProp_FM83_mod, only: &
       & SeaWaterProp_FM83_Temp2PTemp, &
       & SeaWaterProp_FM83_AdLapseRate, &
       & SeaWaterProp_FM83_HeatCapacity

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  
  public :: EOS_JM95_Init, EOS_JM95_Final
  public :: EOS_JM95_Eval
  public :: EOS_JM95_PTemp2Temp, EOS_JM95_Temp2PTemp
  public :: EOS_JM95_AdLapseRate
  public :: EOS_JM95_HeatCapacity
  public :: EOS_JM95_alpha_beta
  
  ! 公開変数
  ! Public variables
  !
  integer, public, parameter :: EOSTYPE_JM95   = 3
  character(TOKEN), public, parameter :: EOSTYPENAME_JM95   = 'EOS_JM95'

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'EOS_JM95_mod' !< Module Name
  
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
  subroutine EOS_JM95_Init()

    ! 実行文; Executable statements
    !

  end subroutine EOS_JM95_Init

  !>
  !!
  !!
  subroutine EOS_JM95_Final()

    ! 実行文; Executable statements
    !

  end subroutine EOS_JM95_Final

  !> @brief Calcualte the density from potential temparature, practical salnity and pressure. 
  !!
  !! @param [in] theta potential temperature [degree C]
  !! @param [in] s     practical salinity    [psu]
  !! @param [in] p     pressure              [bar]
  !
  elemental function EOS_JM95_Eval(theta, s, p) result(rho)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: theta
    real(DP), intent(in) :: s
    real(DP), intent(in) :: p
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

  end function EOS_JM95_Eval

  !> @brief 
  !!
  !!
  !! @param [in] P0     pressure              [bar]
  !! @param [in] S0     practical salinity    [psu]
  !! @param [in] T0     in-situ temperature [degree C]
  !! @param [in] Pref   Reference pressure    [bar]
  !! @return potential temperature [degree C]
  !!
  elemental function EOS_JM95_Temp2PTemp(P0, S0, T0, Pref) result(PTemp)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: P0
    real(DP), intent(in) :: S0
    real(DP), intent(in) :: T0
    real(DP), intent(in) :: Pref
    real(DP) :: PTemp

    ! 実行文; Executable statement
    !

    PTemp = SeaWaterProp_FM83_Temp2PTemp( &
         & P0=bar2dbar(P0), T0=T0, S0=S0, Pref=bar2dbar(Pref) )
    
  end function EOS_JM95_Temp2PTemp

  !> @brief 
  !!
  !!
  !! @param [in] P0     pressure              [bar]
  !! @param [in] S0     practical salinity    [psu]
  !! @param [in] T0     in-situ temperature [degree C]
  !! @param [in] Pref   Reference pressure    [bar]
  !! @return in-situ temperature [degree C]
  !!
  elemental function EOS_JM95_PTemp2Temp(P0, S0, PTemp, Pref) result(InSituTemp)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: P0
    real(DP), intent(in) :: S0
    real(DP), intent(in) :: PTemp
    real(DP), intent(in) :: Pref
    real(DP) :: InSituTemp

    ! 実行文; Executable statement
    !

    InSituTemp = SeaWaterProp_FM83_Temp2PTemp( &
         & P0=bar2dbar(Pref), T0=PTemp, S0=S0, Pref=bar2dbar(P0) )
    
  end function EOS_JM95_PTemp2Temp

  !> @brief 
  !!
  !! @param [in] P     pressure              [bar]
  !! @param [in] S     practical salinity    [psu]
  !! @param [in] T     in-situ temperature [degree C]
  !! @return 
  !!
  function EOS_JM95_HeatCapacity(p, S, T) result(cp)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: p
    real(DP), intent(in) :: S
    real(DP), intent(in) :: T
    real(DP) :: cp
    
    ! 実行文; Executable statement
    !

    cp = SeaWaterProp_FM83_HeatCapacity(S=S, t=T, p=bar2dbar(p))
    

  end function EOS_JM95_HeatCapacity

  !--------------------------------------------------------
  
  !> @brief 
  !!
  !! @param [in] P     pressure              [bar]
  !! @param [in] S     practical salinity    [psu]
  !! @param [in] T     in-situ temperature [degree C]
  !! @return 
  !!
  function EOS_JM95_AdLapseRate(p, S, T) result(adLapseRate)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: p
    real(DP), intent(in) :: S
    real(DP), intent(in) :: T
    real(DP) :: adLapseRate
    
    ! 実行文; Executable statement
    !
    
    adLapseRate = SeaWaterProp_FM83_AdLapseRate(S=S, t=T, p=bar2dbar(p))

  end function EOS_JM95_AdLapseRate

  !--------------------------------------------------------

  !> @brief 
  !!
  !! @param [out] alpha thermal expansion coefficient [degree C^-1]
  !! @param [out] beta  saline contraction coefficient [psu^-1]
  !! @param [in] P      pressure              [dbar]
  !! @param [in] S      practical salinity    [psu]
  !! @param [in] T      in-situ temperature [degree C]
  !!  
  elemental subroutine EOS_JM95_alpha_beta( alpha, beta, &
       & theta, S, p )

    real(DP), intent(out) :: alpha
    real(DP), intent(out) :: beta
    real(DP), intent(in) :: theta
    real(DP), intent(in) :: s
    real(DP), intent(in) :: p

    real(DP) :: alpha_div_beta

    real(DP), parameter :: ab000 =  0.665157d-1
    real(DP), parameter :: ab100 =  0.170907d-1
    real(DP), parameter :: ab200 = -0.203814d-3
    real(DP), parameter :: ab300 =  0.298357d-5
    real(DP), parameter :: ab400 = -0.255019d-7
    real(DP), parameter :: ab010 =  0.378110d-2
    real(DP), parameter :: ab110 = -0.846960d-4
    real(DP), parameter :: ab011 = -0.164759d-6
    real(DP), parameter :: ab012 = -0.251520d-11
    real(DP), parameter :: ab020 = -0.678662d-5
    real(DP), parameter :: ab001 =  0.380374d-4
    real(DP), parameter :: ab101 = -0.933746d-6
    real(DP), parameter :: ab201 =  0.791325d-8
    real(DP), parameter :: ab202 =  0.512857d-12
    real(DP), parameter :: ab003 = -0.302285d-13

    real(DP), parameter :: b000 =  0.785567d-3
    real(DP), parameter :: b100 = -0.301985d-5
    real(DP), parameter :: b200 =  0.555579d-7
    real(DP), parameter :: b300 = -0.415613d-9
    real(DP), parameter :: b010 = -0.356603d-6
    real(DP), parameter :: b110 =  0.788212d-8
    real(DP), parameter :: b011 =  0.408195d-10
    real(DP), parameter :: b012 = -0.602281d-15
    real(DP), parameter :: b020 =  0.515032d-8
    real(DP), parameter :: b001 = -0.121555d-7
    real(DP), parameter :: b101 =  0.192867d-9
    real(DP), parameter :: b201 = -0.213127d-11
    real(DP), parameter :: b002 =  0.176621d-12
    real(DP), parameter :: b102 = -0.175379d-14
    real(DP), parameter :: b003 =  0.121551d-17
    
    alpha_div_beta = ab000 &
         & +      theta*( ab100 + theta*(ab200 + theta*(ab300 + theta*ab400)) )    &
         & + (S - 35d0)*( ab010 + ab110*theta + p*(ab011 + ab012*p)                &
         & +            + (S - 35d0)*ab020                                    )    &
         & +          p*( ab001 + theta*(ab101 + theta*(ab201 + p*ab202))          &
         &                + p*p*ab003                                         )

    
    beta = b000 &
         & +    theta*(b100 + theta*(b200 + theta*b300))                 &
         & + (S-35d0)*(   b010 + b110*theta + p*(b011 + b012*p)          &
         & +            + (S-35d0)*b020                               )  &
         & +        p*(   b001 + theta*(b101 + theta*b201 + p*b102)      &
         &              + p*(b002 + p*b003)                           )

    alpha = alpha_div_beta * beta
    
  end subroutine EOS_JM95_alpha_beta

end module EOS_JM95_mod

