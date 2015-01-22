!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief This module provides a simple nonlinear equation of state for seawater. 
!! 
!! @author Kawai Yuta
!!
!! @detail
!> The simple nonlinear EOS for seawater implemented in this module is proposed by de Szoeke(2004).
!! This simply formulated EOS is used by,for example, Vallis(2006).
!!
!! In the module, the units for temperature, salinity and pressure are 
!! [K], [psu] and [Pa] respectively. 
!!
module EOS_SimpleNonLinear_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: EOS_SimpleNonLinear_Init, EOS_SimpleNonLinear_Final
  public :: EOS_SimpleNonLinear_Eval

  public :: EOS_SimpleNonLinear_PTemp2Temp, EOS_SimpleNonLinear_Temp2PTemp
  public :: EOS_SimpleNonLinear_AdLapseRate
  public :: EOS_SimpleNonLinear_HeatCapacity

  ! 公開変数
  ! Public variable
  !
  integer, public, parameter :: EOSTYPE_SIMPLENONLINEAR = 2
  character(TOKEN), public, parameter :: EOSTYPENAME_SIMPLENONLINEAR = 'EOS_SIMPLENONLINEAR'

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'EOS_SimpleNonLinear_mod' !< Module Name
 
  real(DP) :: RefDens       !< Reference density
  real(DP) :: RefTemp       !< Reference temperature
  real(DP) :: RefSalt       !< Reference salinity
  real(DP) :: Cp0           !< Specific heat capacity at constant pressure
  real(DP) :: BetaT1        !< First thermal expansion coefficient
  real(DP) :: BetaT2        !< Second thermal expansion coefficient
  real(DP) :: BetaS         !< Haline contraction coefficient
  real(DP) :: BetaP         !< Compressibility coefficient
  real(DP) :: Gamma         !< Thermobaric parameter
  real(DP) :: CpBetaS       !< Haline heat capcity coeficient
  real(DP) :: Cs0           !< Reference sound speed


  real(DP) :: PTempDepCoef1, PTempDepCoef2
  real(DP) :: PTempPressCoef
  real(DP) :: PressDepCoef
  real(DP) :: SaltDepCoef

contains

  !>
  !!
  !!
  subroutine EOS_SimpleNonLinear_Init( &
       & refDens_, refTemp_, refSalt_, BetaT1_, BetaT2_, BetaS_, BetaP_, Gamma_, Cp0_, CpBetaS_, Cs0_  )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in), optional :: refDens_    !< Reference density
    real(DP), intent(in), optional :: refTemp_    !< Reference temperature
    real(DP), intent(in), optional :: refSalt_    !< Reference salinity
    real(DP), intent(in), optional :: BetaT1_     !< (First) thermal expansion coffecient
    real(DP), intent(in), optional :: BetaT2_     !< (Second) thermal expansion coffecient
    real(DP), intent(in), optional :: BetaS_      !< Haline contraction cofficient
    real(DP), intent(in), optional :: BetaP_      !< Compressibility coefficient
    real(DP), intent(in), optional :: Gamma_      !< Thermobaric parameter
    real(DP), intent(in), optional :: Cp0_        !< Specific heat capacity at constant pressure
    real(DP), intent(in), optional :: CpBetaS_    !< Haline heat capcity coefficient
    real(DP), intent(in), optional :: Cs0_        !< Reference sound speed

    ! 作業変数
    ! Work variables
    !
    real(DP) :: GammaDash

    ! 実行文; Executable statements
    !

    ! Set the default value of parameters used in calculation of density.  
    ! These values are quoted from Table 1.2 in Vallis(2006). 
    RefDens = 1.027d03  
    RefTemp = 283d0
    RefSalt = 35d0 
    Cp0 = 3986d0        
    BetaT1 = 1.67d-04    
    BetaT2 = 1.00d-05
    BetaS = 0.78d-03
    CpBetaS = 1.5d-03
    BetaP = 4.39d-10
    Gamma = 1.1d-08
    Cs0 = 1490d0  

    if(present(refDens_)) RefDens = refDens_
    if(present(refTemp_)) RefTemp = refTemp_
    if(present(refSalt_)) RefSalt = refSalt_
    if(present(BetaT1_)) BetaT1 = BetaT1_
    if(present(BetaT2_)) BetaT2 = BetaT2_
    if(present(BetaS_)) BetaS = BetaS_
    if(present(Gamma_)) Gamma = Gamma_
    if(present(Cp0_)) Cp0 = Cp0_
    if(present(Cs0_)) Cs0 = Cs0_

    call MessageNotify('M', module_name, &
         & "Set RefDens=%f, RefTemp=%f, RefSalt=%f, BetaT1=%f, BetaT2=%f, BetaS=%f, Gamma=%f, Cp0=%f, Cs0=%f", &
         & d=(/ RefDens, RefTemp, RefSalt, BetaT1, BetaT2, BetaS, Gamma, Cp0, Cs0 /) )


    ! Calculate a coffecient dependent on potential temperature.
    PTempDepCoef1 = - BetaT1
    PTempDepCoef2 = - 0.5d0*BetaT2

    ! Calculate a coffecient dependent on potential temperature and pressure.
    PTempPressCoef = - BetaT1*(Gamma + RefTemp*BetaT2/(RefDens*Cp0))

    ! Calculate a coffecient dependent on salinity.
    SaltDepCoef = BetaS

    ! Calculate a coffecient dependent on pressure.
    ! The fist RHS term is inverse bulk modulus.  
    PressDepCoef = 1d0/(RefDens*Cs0**2) - RefTemp*BetaT1**2/(RefDens*Cp0)


  end subroutine EOS_SimpleNonLinear_Init

  !>
  !!
  !!
  subroutine EOS_SimpleNonLinear_Final()

    ! 実行文; Executable statements
    !

  end subroutine EOS_SimpleNonLinear_Final

  !> @brief 
  !!
  !!
  elemental function EOS_SimpleNonLinear_Eval( theta, s, p ) result( densEdd )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: theta, s, p
    real(DP) :: densEdd

    ! 局所変数
    ! Local variables
    !
    real(DP) :: thetaEdd

    ! 実行文; Executable statement
    !

    thetaEdd = theta - RefTemp

    densEdd = refDens*( & 
         &   (PTempDepCoef1 + PTempPressCoef*p + PTempDepCoef2*thetaEdd)*thetaEdd  &
         & + PressDepCoef*p &
         & + SaltDepCoef*(s - RefSalt) &
         & )

  end function EOS_SimpleNonLinear_Eval

  !> @brief 
  !!
  !!
  !! @param [in] P0     pressure              [Pa]
  !! @param [in] S0     practical salinity    [psu]
  !! @param [in] T0     in-situ temperature [K]
  !! @param [in] Pref   Reference pressure    [Pa]
  !! @return potential temperature [K]
  !!
  elemental function EOS_SimpleNonLinear_Temp2PTemp(P0, S0, T0, Pref) result(PTemp)

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: P0
    real(DP), intent(in) :: S0
    real(DP), intent(in) :: T0
    real(DP), intent(in) :: Pref
    real(DP) :: PTemp

    ! 実行文; Executable statement
    !

    PTemp = RefTemp + &
         & ( (T0 - RefTemp) - (RefTemp*BetaT1/(Cp0*RefDens)*P0*(1d0 + 0.5d0*Gamma*P0)) ) &
         &   / (1d0 + RefTemp*BetaT2/(Cp0*RefDens)*P0)

  end function EOS_SimpleNonLinear_Temp2PTemp

  !> @brief 
  !!
  !!
  !! @param [in] P0      pressure              [Pa]
  !! @param [in] S0      practical salinity    [psu]
  !! @param [in] PTemp   potential temperature [K]
  !! @param [in] Pref    Reference pressure    [Pa]
  !! @return in-situ temperature [K]
  !!
  elemental function EOS_SimpleNonLinear_PTemp2Temp(P0, S0, PTemp, Pref) result(InSituTemp)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: P0
    real(DP), intent(in) :: S0
    real(DP), intent(in) :: PTemp
    real(DP), intent(in) :: Pref
    real(DP) :: InSituTemp

    ! 実行文; Executable statement
    !
    
    InSituTemp = RefTemp + &
         &    RefTemp*BetaT1/(Cp0*RefDens)*P0*(1d0 + 0.5d0*Gamma*P0) &
         & + (PTemp - RefTemp)*(1d0 + RefTemp*BetaT2/(Cp0*RefDens)*P0)
    
  end function EOS_SimpleNonLinear_PTemp2Temp

  !> @brief 
  !!
  !! @return 
  !!
  function EOS_SimpleNonLinear_HeatCapacity(p, S, T) result(cp)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: p
    real(DP), intent(in) :: S
    real(DP), intent(in) :: T
    real(DP) :: cp
    
    ! 実行文; Executable statement
    !
    
    cp = Cp0*(1d0 - CpBetaS*(S - RefSalt)) - BetaT2/RefDens*p*T

  end function EOS_SimpleNonLinear_HeatCapacity

  !> @brief 
  !!
  !! @return 
  !!
  function EOS_SimpleNonLinear_AdLapseRate(p, S, T) result(adLapseRate)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: p
    real(DP), intent(in) :: S
    real(DP), intent(in) :: T
    real(DP) :: adLapseRate
    
    ! 作業変数
    ! Work variables
    !
    real(DP) :: Cp

    ! 実行文; Executable statement
    !
    
    Cp = EOS_SimpleNonLinear_HeatCapacity(p=p, S=S, T=T)
    adLapseRate = T/(Cp*RefDens) * (BetaT1*(1d0 + Gamma*p) + BetaT2*(T - RefTemp))

  end function EOS_SimpleNonLinear_AdLapseRate

end module EOS_SimpleNonLinear_mod
