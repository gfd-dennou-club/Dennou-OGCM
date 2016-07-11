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

end module EOS_JM95_mod

