!-------------------------------------------------------------
! Copyright (c) 2013-2015 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief Thermodynamics process of sea-ice model formulated by Winton(2000). 
!! @author Yuta Kawai
!!
!! This module provide some procedures to advance thermodynamics process of the sea-ice model.
!! The themodynamic model implemented here is formulated in Winton(2000). 
!!
!! To advance thermodynmics process you just should call four subroutines(update_effConductiveCoupling, calc_layersTemprature,
!! calc_SnowIceLyrMassChange and adjust_IceLyrInternal). If you consider stand-alone sea-ice thermodynamics
!! model run, and then the prescribed flux at boundary(i.e. We need not to calculate of the flux at boundary),
!! two interfaces(advance_SeaIceThermDynProc, update_SeaIceThermDynProcVar) will be convinence. 
!!
module SeaIceThermDyn_Winton2000_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN

  use SeaIceConstants_mod, only: &
       & SBConst, &
       & DensIce, DensSnow, DensSeaWater, &
       & CIce, LFreeze, &
       & KIce, KSnow, &
       & Mu, SaltSeaIce, FreezeTempSW, &
       & AlbedoSnow, AlbedoMeltSnow, AlbedoIce, &
       & I0
       
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: SeaIceThermDyn_Winton2000_Init, SeaIceThermDyn_Winton2000_Final

  ! Three subroutines which should be called in order to advance thermdynamics process
  ! of the sea-ice model.
  public :: update_effConductiveCoupling, calc_layersTemprature, calc_SnowIceLyrMassChange, adjust_IceLyrInternal
  
  ! Subroutines which are called in the stand-alone  sea-ice thermodynamics model run.
  public :: advance_SeaIceThermDynProc, update_SeaIceThermDynProcVar

  ! Service routines
  !
  public :: calc_E_IceLyr1, calc_E_IceLyr2
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 公開変数
  ! Public variable
  !

  !> A set of prognostic and diagnostic variables, and parameters.
  !! The derived type will be used to manage these variables in the
  !! stand-alone mode.
  !!
  type, public :: SeaIceThermDynVarSet
     real(DP) :: hsA, hsN
     real(DP) :: hiA, hiN
     real(DP) :: T1A, T1N
     real(DP) :: T2A, T2N
     real(DP) :: TsA, TsN

     real(DP) :: albedo_ice, albedo_snow
  end type SeaIceThermDynVarSet
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SeaIceThermDyn_Winton2000_mod' !< Module Name
  
contains

  !>
  !!
  !!
  subroutine SeaIceThermDyn_Winton2000_Init()

    ! 実行文; Executable statements
    !

  end subroutine SeaIceThermDyn_Winton2000_Init

  !>
  !!
  !!
  subroutine SeaIceThermDyn_Winton2000_Final()

    ! 実行文; Executable statements
    !

  end subroutine SeaIceThermDyn_Winton2000_Final

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_SeaIceThermDynProcVar(varset)

    ! 宣言文; Declaration statement
    !
    type(SeaIceThermDynVarSet), intent(inout) :: varset

    ! 実行文; Executable statement
    !    
    varset%hsN = varset%hsA
    varset%hiN = varset%hiA
    varset%TsN = varset%TsA
    varset%T1N = varset%T1A
    varset%T2N = varset%T2A
    
  end subroutine update_SeaIceThermDynProcVar

  subroutine advance_SeaIceThermDynProc(varset,  & ! (inout)
       & DelTime, Q_SI, Q_LI, SW, LW, Fb,        & ! (in)
       & excessMeltEnergy                        & ! (out)
       & )

    ! 宣言文; Declaration statement
    !    
    type(SeaIceThermDynVarSet), intent(inout) :: varset
    real(DP), intent(in) :: DelTime
    real(DP), intent(in) :: Q_SI, Q_LI, SW, LW, Fb
    real(DP), intent(out) :: excessMeltEnergy
    
    ! 局所変数
    ! Local variables
    !    
    real(DP) :: K_12, K_23
    real(DP) :: A, B
    real(DP) :: dhs, dh1, dh2
    real(DP) :: SurfMeltEn, BtmMeltEn
    real(DP) :: dummyTempO
    logical :: isSnowLyrExist, isSnowLyrMelt, isIceLyrMelt
    
    !> The flux of penetrating solar radiation(if there is no snow). 
    real(DP) :: I

    ! 実行文; Executable statements
    !
    
    !
    isSnowLyrMelt = .false.
    isIceLyrMelt = .false.
    if(varset%hsN <= 0d0) then
       isSnowLyrExist = .false.
    else
       isSnowLyrExist = .true.
    end if
    
    !
    call update_albedo(varset%albedo_snow, varset%albedo_ice, & ! (out)
         & varset%TsN, varset%hsN )

    !
    if(isSnowLyrExist) then
       call update_SurfaceUpwardFluxInfo(A, B, I, &  !(out)
            & varset%TsN, varset%albedo_snow, Q_SI, Q_LI, SW, LW, isSnowLyrExist &  !(in)
            & )
    else
       call update_SurfaceUpwardFluxInfo(A, B, I, &  !(out)
            & varset%TsN, varset%albedo_ice, Q_SI, Q_LI, SW, LW, isSnowLyrExist &  !(in)
            & )
    end if

!!$    write(*,*) "K_1/2=", K_12, ", K_3/2=", K_23
!!$    write(*,*) "A=", A, ", B=", B, ", I=", I
!!$    write(*,*) "K_1/2(T1 - Ts)=", K_12*(varset%T1N - varset%TsN)

    !
    call update_effConductiveCoupling(K_12, K_23, & ! (out)
         & varset%hsN, varset%hiN )

    BtmMeltEn = Fb
    call calc_layersTemprature( &
         & varset%TsA, varset%T1A, varset%T2A, SurfMeltEn,            &  ! (out)
         & BtmMeltEn,                                                 &  ! (inout)
         & varset%T1N, varset%T2N, varset%hiN, DelTime, A, B,         &  ! (in)
         & I, K_12, K_23, isSnowLyrExist )                               ! (in)

    !
    call calc_SnowIceLyrMassChange( &
         & dhs, dh1, dh2, excessMeltEnergy,                        & ! (out)
         & varset%T1A, varset%T2A, dummyTempO,                     & ! (inout)
         & varset%hsN, varset%hiN, varset%TsA, DelTime,            & ! (in)
         & SurfMeltEn, BtmMeltEn                                   & ! (in)
         & )

    !
    varset%hsA = varset%hsN + dhs
    call adjust_IceLyrInternal( &
         & varset%hsA, varset%hiA, varset%T1A, varset%T2A,                & !(inout)
         & 0.5d0*varset%hiN+dh1, 0.5d0*varset%hiN+dh2 ) !(in)    

  end subroutine advance_SeaIceThermDynProc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update_albedo(albedo_snow, albedo_ice, & ! (out)
       & Ts, hs )                                     ! (in)

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: albedo_snow, albedo_ice
    real(DP), intent(in) :: Ts, hs

    ! 実行文; Executable statements
    !
    
    if(Ts >= 0d0)then
       albedo_snow = AlbedoMeltSnow
    else
       albedo_snow = AlbedoSnow
    end if

    albedo_ice = AlbedoIce
    
  end subroutine update_albedo
  
  subroutine update_effConductiveCoupling(K_12, K_23, &  ! (out)
       & hs, hi                                       &  ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: K_12, K_23
    real(DP), intent(in) :: hs, hi

    ! 実行文; Executable statements
    !
    
    K_12 = 4d0*KIce*KSnow/(KSnow*hi + 4d0*KIce*hs)
    K_23 = 2d0*KIce/hi
    
  end subroutine update_effConductiveCoupling

  subroutine update_SurfaceUpwardFluxInfo(A, B, I, &  !(out)
       & Ts, albedo, Q_SI, Q_LI, SW, LW, isSnowLyrExist &  !(in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: A, B, I
    real(DP), intent(in) :: Ts, albedo, Q_SI, Q_LI, SW, LW
    logical, intent(in) :: isSnowLyrExist

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: Fs, dFsdTs, Ts_UnitK

    ! 実行文; Executable statements
    !
    
    Ts_UnitK = degC2K(Ts)
    Fs = Q_SI + Q_LI + LW + SBConst*(Ts_UnitK)**4
    if(isSnowLyrExist) then
       I = 0d0
       Fs = Fs + (1d0 - albedo)*SW
    else
       I = -I0*(1d0 - albedo)*SW
       Fs = Fs + (1d0 - I0)*(1d0 - albedo)*SW
    end if
    dFsdTs  = 4d0*SBConst*(Ts_UnitK)**3

!!$    write(*,*) "SI=", Q_SI, ", LI=", Q_LI, "SW=", SW, "LW=", LW, &
!!$         & "Sigma*T**4=", SBConst*(Ts_UnitK)**4, ", albedo=", albedo
!!$    write(*,*) "Fs=", Fs, ", dFsdTs=", dFsdTs, "I=", I, isSnowLyrExist
!!$    write(*,*) "A,B", A,b
    
    A = Fs - Ts*dFsdTs
    B = dFsdTs

  end subroutine update_SurfaceUpwardFluxInfo

  subroutine calc_layersTemprature( &
       & TsA, T1A, T2A, SurfMeltEn,                            &  ! (out)
       & BtmMeltEn,                                             &  ! (inout)
       & T1N, T2N, hi, dt, A, B, I, K_12, K_23, isSnowLyrExist &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: TsA, T1A, T2A
    real(DP), intent(out) :: SurfMeltEn, BtmMeltEn
    real(DP), intent(in) :: T1N, T2N, hi, I, dt
    real(DP), intent(in) :: A, B, K_12, K_23
    logical, intent(in) :: isSnowLyrExist

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: A1, B1, C1, A10, B10
    real(DP) :: HCIce, Work1
    logical :: MeltFlag

    ! 実行文; Executable statements
    !
    
    HCIce = DensIce*hi*CIce
    Work1 = 6d0*dt*K_23 + HCIce

    !
    A10 =     HCIce/(2d0*dt) &
         & + K_23*(4d0*dt*K_23 + HCIce)/Work1

    B10 =   - DensIce*hi/(2d0*dt)*(CIce*T1N - LFreeze*Mu*SaltSeaIce/T1N) &
         & - I                                                                        &
         & - K_23*(4d0*dt*K_23*FreezeTempSW + HCIce*T2N)/Work1

    A1 =   A10 + K_12*B/(K_12 + B)
    B1 =   B10 + K_12*A/(K_12 + B) 
    C1  = - DensIce*hi/(2d0*dt)*LFreeze*Mu*SaltSeaIce
    
    T1A = -(B1 + sqrt(B1**2 - 4d0*A1*C1))/(2d0*A1)
    T2A = (2d0*dt*K_23*(T1A + 2d0*FreezeTempSW) + HCIce*T2N)/Work1
    TsA = (K_12*T1A - A)/(K_12 + B)

    MeltFlag = .false.
    if(isSnowLyrExist .and. TsA > 0d0) then
       TsA = 0d0; MeltFlag = .true.
    else if(.not. isSnowLyrExist .and. TsA > -Mu*SaltSeaIce) then
       TsA = -Mu*SaltSeaIce; MeltFlag = .true.
    end if
    
    if(MeltFlag) then
       A1 = A10 + K_12
       B1 = B10 - K_12*TsA
       T1A = -(B1 + sqrt(B1**2 - 4d0*A1*C1))/(2d0*A1)
       T2A = (2d0*dt*K_23*(T1A + 2d0*FreezeTempSW) + HCIce*T2N)/Work1
    end if

    !
    SurfMeltEn = K_12*(T1A - TsA) - (A + B*TsA)
    BtmMeltEn  = BtmMeltEn - 4d0*KIce*(FreezeTempSW - T2A)/hi
    
  end subroutine calc_layersTemprature

  subroutine calc_SnowIceLyrMassChange( &
       & dhs, dh1, dh2, excessMeltEnergy,                 & ! (out)
       & T1, T2, TsOcean,                                 & ! (inout)
       & hsOld_, hiOld_, Ts, dt, SurfMeltEn, BtmMeltEn    & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: dhs, dh1, dh2, excessMeltEnergy
    real(DP), intent(inout) :: T1, T2, TsOcean
    real(DP), intent(in) :: hsOld_, hiOld_, Ts
    real(DP), intent(in) :: dt, SurfMeltEn, BtmMeltEn

    ! 局所変数
    ! Local variables
    !
    real(DP) :: E1, E2
    real(DP) :: hsOld, h1Old, h2Old
    real(DP) :: Work

    real(DP) :: dh2Freeze
    real(DP) :: dhsMelt, dh1Melt, dh2Melt

    ! 実行文; Executable statements
    !
    
    !
    hsOld = hsOld_
    h1Old = 0.5d0*hiOld_; h2Old = 0.5d0*hiOld_
    excessMeltEnergy = 0d0
    
    !* Freezing
    !
    
    dh2Freeze = 0d0
    if(BtmMeltEn <= 0d0) then
       dh2Freeze = BtmMeltEn*dt/(DensIce*calc_E_IceLyr2(FreezeTempSW, SaltSeaIce))
       T2 = (dh2Freeze*FreezeTempSW + h2Old*T2)/(h2Old + dh2Freeze)

       if(hiOld_ == 0d0) then
          if(dh2Freeze > 1d-4) then
             ! The first layer temperature for intialized sea ice is calculated by adjust_IceLyrInternal.
             ! Before calculating the temperature, T1 need to be set *non* zero. 
             T1 = - Mu*SaltSeaIce
          else
             dh2Freeze = 0d0
          end if
       end if

       if(dh2Freeze > 0d0) &
            & TsOcean = degC2K(FreezeTempSW)
    end if
    

    
    !* Melting
    !
    
    dhsMelt = 0d0; dh1Melt = 0d0; dh2Melt = 0d0
    if(SurfMeltEn > 0d0) then
       Work = SurfMeltEn*dt
       dhsMelt = - min(Work/(DensSnow*LFreeze), hsOld)

       E1 = calc_E_IceLyr1(T1, SaltSeaIce)
       Work = Work - DensSnow*LFreeze*hsOld
       dh1Melt = - min( max(-Work/(DensIce*E1), 0d0), h1Old)
       
       E2 = calc_E_IceLyr2(T2, SaltSeaIce)
       Work = Work + DensIce*E1*h1Old
       dh2Melt = - min(max(-Work/(DensIce*E2), 0d0), h2Old)

       excessMeltEnergy = max(Work + DensIce*E2*h2Old, 0d0)

       hsOld = hsOld + dhsMelt       
       h1Old = h1Old + dh1Melt
       h2Old = h2Old + dh2Melt       
    end if
    
    if(BtmMeltEn > 0d0) then
       E2 = calc_E_IceLyr2(T2, SaltSeaIce)
       Work =BtmMeltEn*dt
       dh2Melt = dh2Melt - min(-Work/(DensIce*E2), h2Old)
       
       E1 = calc_E_IceLyr1(T1, SaltSeaIce)
       Work = Work + DensIce*E2*h2Old
       dh1Melt = dh1Melt - min( max(-Work/(DensIce*E1), 0d0), h1Old )

       Work = Work + DensIce*E1*h1Old
       dhsMelt = dhsMelt - min( max(Work/(DensSnow*LFreeze), 0d0), hsOld )

       excessMeltEnergy = excessMeltEnergy + max(Work - DensSnow*LFreeze*hsOld, 0d0)

    end if

    dhs = dhsMelt
    dh1 = dh1Melt
    dh2 = dh2Freeze + dh2Melt
    
!!$    write(*,*) "dhs=", dhs, ", dh1=", dh1, ", dh2=", dh2
!!$    write(*,*) "T1New=", T1, ", T2New=", T2
!!$    write(*,*) "excessMeltEnergy=", excessMeltEnergy

  end subroutine calc_SnowIceLyrMassChange

  subroutine adjust_IceLyrInternal( &
       & hs, hi, T1, T2,                & !(inout)
       & h1, h2 ) !(in)

    real(DP), intent(inout) :: hs, hi, T1, T2
    real(DP), intent(in) :: h1, h2

    ! 局所変数
    ! Local variables
    !    
    Real(DP) :: dhs, dh1, h1_
    real(DP) :: coef, dummy, T1Tmp, T2Tmp, h1Tmp
    real(DP) :: extraSE, dh

    
    !
    !
    h1_ = h1
    hi = h1_ + h2    
    
    !
    !
    if(hs > 0d0) then
       coef = (hs - (DensSeaWater - DensIce)/DensSnow*hi)/DensSeaWater
       dhs = - max(coef*DensIce, 0d0)
       dh1 =   max(coef*DensSeaWater, 0d0)

       !
       if(dh1 > 0d0) then
          hs = hs + dhs; h1_ = h1 + dh1
          call calc_NewTemp(T1Tmp, dummy, &
               & T1, -Mu*SaltSeaIce, h1/h1_)
          T1 = T1Tmp
       end if
!!$    write(*,*) "Adjust the snow layer below waterline dhs=", dhs, ", dh1=", dh1, ", hi=", hi

    end if


    
    !
    !


    if(h1_ < h2) then
!       write(*,*) "h1 < h2"
       call calc_NewTemp(T1Tmp, dummy, &
            & T1, T2, h1_/(0.5d0*hi))
       T1 = T1Tmp
    else
!       write(*,*) "h1 >= h2", h1Tmp, h2
       call calc_NewTemp(dummy, T2Tmp, &
            & T1, T2, h1_/(0.5d0*hi)-1d0)
       T2 = T2Tmp
       
       if(T2 > - Mu*SaltSeaIce) then

          T2 = - Mu*SaltSeaIce
          extraSE = 0.5d0*hi*( &
               & calc_E_IceLyr2(T2, SaltSeaIce) - calc_E_IceLyr2(-Mu*SaltSeaIce,SaltSeaIce) &
               & )
          dh = 2d0*extraSE/(calc_E_IceLyr1(T1,SaltSeaIce) + calc_E_IceLyr2(-Mu*SaltSeaIce,SaltSeaIce))
          hi = hi + dh
       end if
    end if


!!$    write(*,*) "hs=", hs, "hi=", hi
!!$    write(*,*) "Adjust temp T1=", T1, ", T2=", T2, "dummy=", dummy

  contains
    
    subroutine calc_NewTemp(T1New, T2New, &
         & T1, T2, f1 )
      ! 宣言文; Declaration statement
      !      
      real(DP), intent(out) :: T1New, T2New
      real(DP), intent(in) :: T1, T2, f1

      ! 実行文; Executable statements
      !
      
      T2New =   f1*(T1 - LFreeze/CIce*Mu*SaltSeaIce/T1) &
           &  + (1d0 - f1)*T2
      T1New = 0.5d0*(T2New - sqrt(T2New**2 + 4d0*Mu*SaltSeaIce*LFreeze/CIce))
    end subroutine calc_NewTemp

  end subroutine adjust_IceLyrInternal

!!!!!!!!!!!!!
  
  function calc_E_IceLyr1(T, S) result(E)
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: T, S
    real(DP) :: E

    ! 実行文; Executable statements
    !
    
    E =     CIce*(T + Mu*S) &
         & - LFreeze*(1d0 + Mu*S/T)
    
  end function calc_E_IceLyr1
  
  function calc_E_IceLyr2(T, S) result(E2)
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: T, S
    real(DP) :: E2

    ! 実行文; Executable statements
    !
    
    E2 = CIce*(T + Mu*S) - LFreeze
  end function calc_E_IceLyr2
  
  function degC2K(degC) result(K)
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: degC
    real(DP) :: K

    ! 実行文; Executable statements
    !

    K = 273.15d0 + degC
  end function degC2K
  
  
end module SeaIceThermDyn_Winton2000_mod

