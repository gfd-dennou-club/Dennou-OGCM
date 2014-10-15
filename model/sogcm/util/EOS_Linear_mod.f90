!-------------------------------------------------------------
! Copyright (c) 2013-2013 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!! @detail
!> The linear EOS for seawater implemented in this module is described in Vallis(2006).
!! In the module, the units for temperature, salinity and pressure are 
!! [K], [psu] and [Pa] respectively. 
!!
module EOS_Linear_mod 

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
  public :: EOS_Linear_Init, EOS_Linear_Final
  public :: EOS_Linear_Eval

  public :: EOS_Linear_PTemp2Temp, EOS_Linear_Temp2PTemp
  public :: EOS_Linear_AdLapseRate
  public :: EOS_Linear_HeatCapacity

  ! 公開変数
  ! Public variable
  !
  integer, public, parameter :: EOSTYPE_LINEAR = 1
  character(TOKEN), public, parameter :: EOSTYPENAME_LINEAR = 'EOS_LINEAR'

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'EOS_Linear_mod' !< Module Name
 
  real(DP) :: RefDens      !< Reference density
  real(DP) :: RefTemp      !< Reference temperature
  real(DP) :: RefSalt      !< Reference salinity
  real(DP) :: Cp0          !< Specific heat capacity at constant pressure
  real(DP) :: BetaT        !< (First) thermal expansion coefficient
  real(DP) :: BetaS        !< Haline contraction coefficient
  real(DP) :: CpBetaS      !< Haline heat capcity coeficient
  real(DP) :: Cs0          !< Reference sound speed

  real(DP) :: PTempDepCoef
  real(DP) :: PressDepCoef
  real(DP) :: SaltDepCoef

contains

  !>
  !!
  !!
  subroutine EOS_Linear_Init( &
       & refDens_, refTemp_, refSalt_, BetaT_, BetaS_, Cp0_, CpBetaS_, Cs0_  )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(in), optional :: refDens_     !< Reference density
    real(DP), intent(in), optional :: refTemp_     !< Reference temperature
    real(DP), intent(in), optional :: refSalt_     !< Reference salinity
    real(DP), intent(in), optional :: BetaT_       !< (First) thermal expansion coefficient
    real(DP), intent(in), optional :: BetaS_       !< Haline contraction cofficient
    real(DP), intent(in), optional :: Cp0_         !< Specific heat capacity at constant pressure
    real(DP), intent(in), optional :: CpBetaS_     !< Haline heat capcity coefficient
    real(DP), intent(in), optional :: Cs0_         !< Reference sound speed

    ! 作業変数
    ! Work variables
    !

    ! 実行文; Executable statements
    !

    ! Set the default value of parameters used in linear EOS.  
    ! These values are quoted from Table 1.2 in Vallis(2006). 
    RefDens = 1.027d03  
    RefTemp = 283d0          
    RefSalt = 35d0 
    Cp0 = 3986d0        
    BetaT = 1.67d-04
    BetaS = 0.78d-03
    CpBetaS = 1.5d-03
    Cs0 = 1490d0  

    if(present(refDens_)) RefDens = refDens_
    if(present(refTemp_)) RefTemp = refTemp_
    if(present(refSalt_)) RefSalt = refSalt_
    if(present(BetaT_)) BetaT = BetaT_
    if(present(BetaS_)) BetaS = BetaS_
    if(present(Cp0_)) Cp0 = Cp0_
    if(present(CpBetaS_)) CpBetaS = CpBetaS_
    if(present(Cs0_)) Cs0 = Cs0_

    call MessageNotify('M', module_name, &
         & "Set RefDens=%f, RefTemp=%f, RefSalt=%f, BetaT=%f, BetaS=%f, Cp0=%f, Cs0=%f", &
         & d=(/ RefDens, RefTemp, RefSalt, BetaT, BetaS, Cp0, Cs0 /) )

    PTempDepCoef = - BetaT
    PressDepCoef = 1d0/(RefDens*Cs0**2) - RefTemp*BetaT**2/(RefDens*Cp0)
    SaltDepCoef = BetaS

  end subroutine EOS_Linear_Init

  !>
  !!
  !!
  subroutine EOS_Linear_Final()

    ! 実行文; Executable statements
    !

  end subroutine EOS_Linear_Final

  !> @brief 
  !!
  !!
  elemental function EOS_Linear_Eval( theta, s, p ) result( densEdd )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: theta, s, p
    real(DP) :: densEdd

    ! 局所変数
    ! Local variables
    !

    ! 実行文; Executable statement
    !

    densEdd = refDens*( & 
         &  PTempDepCoef*(theta - RefTemp)  &
         & + PressDepCoef*p &
         & + SaltDepCoef*(s - RefSalt) &
         & )

  end function EOS_Linear_Eval

  !> @brief 
  !!
  !!
  !! @param [in] P0     pressure              [Pa]
  !! @param [in] S0     practical salinity    [psu]
  !! @param [in] T0     in-situ temperature [K]
  !! @param [in] Pref   Reference pressure    [Pa]
  !! @return potential temperature [K]
  !!
  elemental function EOS_Linear_Temp2PTemp(P0, S0, T0, Pref) result(PTemp)

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
         & ( (T0 - RefTemp) - RefTemp*BetaT/(Cp0*RefDens)*P0 )

  end function EOS_Linear_Temp2PTemp

  !> @brief 
  !!
  !!
  !! @param [in] P0      pressure              [Pa]
  !! @param [in] S0      practical salinity    [psu]
  !! @param [in] PTemp   potential temperature [K]
  !! @param [in] Pref    Reference pressure    [Pa]
  !! @return in-situ temperature [K]
  !!
  elemental function EOS_Linear_PTemp2Temp(P0, S0, PTemp, Pref) result(InSituTemp)
    
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
         &    RefTemp*BetaT/(Cp0*RefDens)*P0 &
         & + (PTemp - RefTemp)
    
  end function EOS_Linear_PTemp2Temp

  !> @brief 
  !!
  !! @return 
  !!
  function EOS_Linear_HeatCapacity(p, S, T) result(cp)
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: p
    real(DP), intent(in) :: S
    real(DP), intent(in) :: T
    real(DP) :: cp
    
    ! 実行文; Executable statement
    !
    
    cp = Cp0*(1d0 - CpBetaS*(S - RefSalt))

  end function EOS_Linear_HeatCapacity

  !> @brief 
  !!
  !! @return 
  !!
  function EOS_Linear_AdLapseRate(p, S, T) result(adLapseRate)
    
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
    
    Cp = EOS_Linear_HeatCapacity(p=p, S=S, T=T)
    adLapseRate = BetaT*T/(Cp*RefDens)

  end function EOS_Linear_AdLapseRate

end module EOS_Linear_mod
