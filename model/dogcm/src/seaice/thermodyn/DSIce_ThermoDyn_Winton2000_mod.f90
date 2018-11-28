!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief A module to construct a 3 layer thermodynamical sea-ice model by Winton(2000).
!! 
!! @author Kawai Yuta
!!
!!
module DSIce_ThermoDyn_Winton2000_mod
     
  ! モジュール引用; Use statements
  !
 
  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN

  use dc_message, only: &
       & MessageNotify
  
  !* Dennou-SIce

  use DOGCM_Admin_Constants_mod, only: &
       & LEvap => LatentHeat
  
  use DSIce_Admin_Constants_mod, only: &
       & SBConst,                               &
       & DensIce, DensSnow, DensSeaWater,       &
       & CIce, LFreeze,                         &
       & KIce, KSnow,                           &
       & Mu, SaltSeaIce, FreezeTempSW,          &
       & AlbedoSnow, AlbedoMeltSnow, AlbedoIce, &
       & I0,                                    &
       & IceThickMin

#include "../DSIceDef.h"
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DSIce_ThermoDyn_Winton2000_Init, DSIce_ThermoDyn_Winton2000_Final

  public :: DSIce_ThermoDyn_Winton2000_CalcIceTemp
  public :: DSIce_ThermoDyn_Winton2000_CalcLyrMassChange
  public :: DSIce_ThermoDyn_Winton2000_AdjustLyrInternal

  public :: calc_E_IceLyr1
  public :: calc_E_IceLyr2
  public :: calc_Temp_IceLyr1
  public :: calc_Temp_IceLyr2
  public :: calc_SIceTotEn
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_ThermoDyn_Winton2000_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DSIce_ThermoDyn_Winton2000_Init()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_ThermoDyn_Winton2000_Init

  !>
  !!
  !!
  subroutine DSIce_ThermoDyn_Winton2000_Final()

    ! 実行文; Executable statements
    !

  end subroutine DSIce_ThermoDyn_Winton2000_Final

  !----------------------------------------------

  subroutine DSIce_ThermoDyn_Winton2000_CalcIceTemp( &
       & SfcTempA, z_IceTempA, SfcResHFlx, BtmResHFlx,   &  ! (out)
       & SfcTempN, z_IceTempN, IceThickN, SnowThickN,    &  ! (in)
       & SfcHFlx, DSfcHFlxDTs, PenSWRFlx, SfcFrzTemp,    &  ! (in)
       & BtmHFlx, OcnFrzTemp,                            &  ! (in)
       & dt,                                             &  ! (in)
       & debugFlag                                       &  ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: SfcTempA
    real(DP), intent(out) :: z_IceTempA(2)
    real(DP), intent(out) :: SfcResHFlx
    real(DP), intent(out) :: BtmResHFlx
    real(DP), intent(in) :: SfcTempN
    real(DP), intent(in) :: z_IceTempN(2)
    real(DP), intent(in) :: IceThickN
    real(DP), intent(in) :: SnowThickN
    real(DP), intent(in) :: SfcHFlx
    real(DP), intent(in) :: DSfcHFlxDTs
    real(DP), intent(in) :: PenSWRFlx
    real(DP), intent(in) :: SfcFrzTemp
    real(DP), intent(in) :: BtmHFlx    
    real(DP), intent(in) :: OcnFrzTemp
    real(DP), intent(in) :: dt
    logical, intent(in) :: debugFlag
       
    ! 局所変数
    ! Local variables
    !
    real(DP) :: K_12
    real(DP) :: K_23
    real(DP) :: A1
    real(DP) :: A10
    real(DP) :: B1
    real(DP) :: B10
    real(DP) :: C1
    real(DP) :: HCIce
    real(DP) :: Work1
    real(DP) :: A
    real(DP) :: B
    real(DP) :: SIceMeltTemp
    
#ifdef DEBUG_SEAICE
    real(DP) :: E1
    real(DP) :: E2
#endif
    
    ! 実行文; Executable statements
    !

    SIceMeltTemp = - Mu*SaltSeaIce
    
    K_12 = 4d0*KIce*KSnow/(KSnow*IceThickN + 4d0*KIce*SnowThickN)
    K_23 = 2d0*KIce/IceThickN
      
    HCIce = DensIce * IceThickN * CIce
    Work1 = 6d0 * dt * K_23 + HCIce
 
    B = DSfcHFlxDTs
    A = SfcHFlx - SfcTempN*DSfcHFlxDTs

    !--------------------------------------------------
    
    A10 =     HCIce/(2d0*dt) &
         & + K_23*(4d0*dt*K_23 + HCIce)/Work1

    B10 =  - DensIce*IceThickN/(2d0*dt) &
         &      *(CIce*z_IceTempN(1) + LFreeze*SIceMeltTemp/z_IceTempN(1))            &
         & - PenSWRFlx                                                                &
         & - K_23*(4d0*dt*K_23*OcnFrzTemp + HCIce*z_IceTempN(2))/Work1

    A1 =   A10 + K_12*B/(K_12 + B)
    B1 =   B10 + K_12*A/(K_12 + B) 
    C1  =  DensIce*IceThickN/(2d0*dt)*LFreeze*SIceMeltTemp
     
    z_IceTempA(1) = -(B1 + sqrt(B1**2 - 4d0*A1*C1))/(2d0*A1)
    z_IceTempA(2) = (2d0*dt*K_23*(z_IceTempA(1) + 2d0*OcnFrzTemp) + HCIce*z_IceTempN(2))/Work1
    SfcTempA = (K_12*z_IceTempA(1) - A)/(K_12 + B)

    if (SfcTempA > SfcFrzTemp) then
       SfcTempA = SfcFrzTemp
  
       A1 = A10 + K_12
       B1 = B10 - K_12*SfcTempA
       z_IceTempA(1) = -(B1 + sqrt(B1**2 - 4d0*A1*C1))/(2d0*A1)
       z_IceTempA(2) = (2d0*dt*K_23*(z_IceTempA(1) + 2d0*OcnFrzTemp) + HCIce*z_IceTempN(2))/Work1
    end if
 
    !
    SfcResHFlx = K_12*(z_IceTempA(1) - SfcTempA) - (A + B*SfcTempA)
    BtmResHFlx = BtmHFlx - 4d0*KIce*(OcnFrzTemp - z_IceTempA(2))/IceThickN

    !------------------------------------------
           

#ifdef DEBUG_SEAICE
    if ( debugFlag ) then
!!$       write(*,*) "Lyr1: LHS=", DensIce*IceThickN/(2d0*dt)*(CIce + LFreeze*Mu*SaltSeaIce/(z_IceTempA(1)*z_IceTempN(1))) &
!!$            & *(z_IceTempA(1) - z_IceTempN(1))
!!$       write(*,*) "Lyr1: RHS=", K_12*(SfcTempA - z_IceTempA(1)) + K_23*(z_IceTempA(2) - z_IceTempA(1)) + PenSWRFlx
!!$
!!$       write(*,*) "Lyr2: LHS=", DensIce*IceThickN/(2d0*dt)*CIce*(z_IceTempA(2) - z_IceTempN(2))
!!$       write(*,*) "Lyr2: RHS=", K_23*(z_IceTempA(1) - z_IceTempA(2)) + 2d0*K_23*(OcnFrzTemp - z_IceTempA(2))
!!$       
!!$       write(*,*) "CalcIceTemp (Check HFlx):", &
!!$            & "SfcIceHFlxA=", K_12*(z_IceTempA(1) - SfcTempA), &
!!$            & "BtmIceHFlx=", 4d0*KIce*(OcnFrzTemp - z_IceTempA(2))/IceThickN
!!$            
!!$       write(*,*) "CalcIceTemp (Check IceEngy): (SfcHFlx-BtmHFl)*dtx=", &
!!$            & dt*( - K_12*(z_IceTempA(1) - SfcTempA)                     &
!!$            &      + 4d0*KIce*(OcnFrzTemp - z_IceTempA(2))/IceThickN     &
!!$            & )
       E1 = (calc_E_IceLyr1(z_IceTempA(1), SaltSeaIce) - calc_E_IceLyr1(z_IceTempN(1), SaltSeaIce))
       E2 = (calc_E_IceLyr2(z_IceTempA(2), SaltSeaIce) - calc_E_IceLyr2(z_IceTempN(2), SaltSeaIce))
       write(*,*) "CalcIceTemp check:",              &
            & 0.5d0*IceThickN*DensIce*(E1 + E2),                  &
            & - dt*(   ((A + B*SfcTempA) - PenSWRFlx - BtmHFlx)   &
            &        + ( BtmResHFlx + SfcResHFlx  )               &
            & )
    end if
#endif

  end subroutine DSIce_ThermoDyn_Winton2000_CalcIceTemp

  !-----------------------------------------------

  subroutine DSIce_ThermoDyn_Winton2000_CalcLyrMassChange( &
       & dSnowThick, dIceThick1, dIceThick2, wice, excessMeltEn,   & ! (out)
       & z_IceTemp,                                                & ! (inout)
       & SnowThick0, IceThick0,                                    & ! (in)
       & SfcResHFlx, BtmResHFlx, OcnFrzTemp,                       & ! (in)
       & RainFall, SnowFall, Evap,                                 & ! (in)
       & dt,                                                       & ! (in)
       & debugFlag                                                 & ! (in)
       & )

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: dSnowThick
    real(DP), intent(out) :: dIceThick1
    real(DP), intent(out) :: dIceThick2
    real(DP), intent(out) :: wice
    real(DP), intent(out) :: excessMeltEn
    real(DP), intent(inout) :: z_IceTemp(2)    
    real(DP), intent(in) :: SnowThick0
    real(DP), intent(in) :: IceThick0
    real(DP), intent(in) :: SfcResHFlx
    real(DP), intent(in) :: BtmResHFlx !< [W/m2]
    real(DP), intent(in) :: OcnFrzTemp !< [degC]
    real(DP), intent(in) :: RainFall   !< [kg/(m2.s)]
    real(DP), intent(in) :: SnowFall   !< [kg/(m2.s)]
    real(DP), intent(in) :: Evap       !< [kg/(m2.s)]    
    real(DP), intent(in) :: dt
    logical, intent(in) :: debugFlag

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: E1
    real(DP) :: E2
    real(DP) :: E1Old
    real(DP) :: E2Old
    real(DP) :: hEtmp
    
    real(DP) :: hsOld
    real(DP) :: h1Old
    real(DP) :: h2Old

    real(DP) :: work
    real(DP) :: dh2Frz
    real(DP) :: dhsMelt
    real(DP) :: dh1Melt
    real(DP) :: dh2Melt
    real(DP) :: DensQTmp
       
    real(DP) :: dhsEvap
    real(DP) :: dh1Evap
    real(DP) :: dh2Evap

    
    ! 実行文; Executable statements
    !
    
    hsOld = SnowThick0
    h1Old = 0.5d0*IceThick0
    h2Old = 0.5d0*IceThick0
    E1Old = calc_E_IceLyr1(z_IceTemp(1), SaltSeaIce)
    E2Old = calc_E_IceLyr2(z_IceTemp(2), SaltSeaIce)

    dSnowThick = 0d0
    dIceThick1 = 0d0
    dIceThick2 = 0d0
    wice = 0d0
    excessMeltEn = 0d0
    
    ! Freezing -----------------------------------------------------------------
    
    if ( BtmResHFlx < 0d0 ) then
       dh2Frz = BtmResHFlx * dt / (DensIce*calc_E_IceLyr2(OcnFrzTemp, SaltSeaIce)) ! > 0
       dIceThick2 = dIceThick2 + dh2Frz
       wice = wice + DensIce*dh2Frz / dt
       z_IceTemp(2) = (dh2Frz*OcnFrzTemp + h2Old*z_IceTemp(2))/(h2Old + dh2Frz)
    end if
       
    ! Melting (added evaporation on the surface) -------------------------------
 
    E1 = calc_E_IceLyr1(z_IceTemp(1), SaltSeaIce)
    E2 = calc_E_IceLyr2(z_IceTemp(2), SaltSeaIce)

    if ( SfcResHFlx > 0d0 ) then

       work = SfcResHFlx*dt
       dhsMelt = 0d0
       dh1Melt = 0d0
       dh2Melt = 0d0

       ! Melt the snow layer
       DensQTmp = DensSnow*LFreeze
       dhsMelt = min(work/DensQTmp, hsOld + dSnowThick)
       dSnowThick = dSnowThick - dhsMelt
       work = work - DensQTmp*dhsMelt

       ! Melt the upper ice layer
       if (work > 0d0) then
          DensQTmp = - DensIce*E1
          if (work < (h1Old + dIceThick1)*DensQTmp) then
             dh1Melt = work/DensQTmp
          else
             dh1Melt = h1Old + dIceThick1
          end if
!!$       dh1Melt = min(max(-work/(DensIce*E1), 0d0), h1Old+dIceThick1)
          dIceThick1 = dIceThick1 - dh1Melt          
          work = work - DensQTmp*dh1Melt
       end if
       
       ! Melt the lower ice layer
       if (work > 0d0) then
          DensQTmp = - DensIce*E2
          if (work < (h2Old + dIceThick2)*DensQTmp) then
             dh2Melt = work/DensQTmp
          else
             dh2Melt = h2Old + dIceThick2
          end if
!!$       dh2Melt = min( max(-work/(DensIce*E2), 0d0), h2Old+dIceThick2)
          dIceThick2 = dIceThick2 - dh2Melt
          work = work - DensQTmp*dh2Melt
       end if

       !
       excessMeltEn =   excessMeltEn + work 
       wice =   wice &
            & - (DensSnow*dhsMelt + DensIce*(dh1Melt + dh2Melt)) / dt
    end if

#ifdef DEBUG_SEAICE
    if (debugFlag) then
       write(*,*) "(Surface) dhsMelt=", dhsMelt, ", dh1Melt=", dh1Melt, ", dh2Melt=", dh2Melt, &
            & "melt excessMeltEn=", excessMeltEn
       dhsMelt = 0d0; dh1Melt = 0d0; dh2Melt = 0d0
    end if
#endif ! <- DEBUG_SEAICE
    
    if ( BtmResHFlx > 0d0 ) then

       work = BtmResHFlx*dt
       dhsMelt = 0d0
       dh1Melt = 0d0
       dh2Melt = 0d0

       ! Melt the lower ice layer
       DensQTmp = - DensIce*E2
       if (work < (h2Old + dIceThick2)*DensQTmp) then
          dh2Melt = work/DensQTmp
       else
          dh2Melt = h2Old + dIceThick2
       end if
!!$       dh2Melt = min( max(-work/(DensIce*E2), 0d0), h2Old+dIceThick2)
       dIceThick2 = dIceThick2 - dh2Melt
       work = work - DensQTmp*dh2Melt

       ! Melt the upper  ice layer
       if (work > 0d0) then
          DensQTmp = - DensIce*E1
          if (work < (h1Old + dIceThick1)*DensQTmp) then
             dh1Melt = work/DensQTmp
          else
             dh1Melt = h1Old+dIceThick1
          end if
       end if
!!$       dh1Melt = min( max(-work/(DensIce*E1), 0d0), h1Old+dIceThick1)
       dIceThick1 = dIceThick1 - dh1Melt
       work = work - DensQTmp*dh1Melt

       ! Melt the snow layer
       if (work > 0d0) then
          DensQTmp = DensSnow*LFreeze
          if (work < (hsOld + dSnowThick)*DensQTmp) then
             dhsMelt = work/DensQTmp
          else
             dhsMelt = hsOld + dSnowThick
          end if
       end if
!!$       dhsMelt = min( max(work/(DensSnow*LFreeze), 0d0), hsOld+dSnowThick)
       dSnowThick = dSnowThick - dhsMelt
       work = work - DensQTmp*dhsMelt
  
       excessMeltEn =   excessMeltEn + work
       wice =   wice &
            & - (DensSnow*dhsMelt + DensIce*(dh1Melt + dh2Melt)) / dt
    end if
    
#ifdef DEBUG_SEAICE
    if (debugFlag) then
       write(*,*) "(bottom) dhsMelt=", dhsMelt, ", dh1Melt=", dh1Melt, ", dh2Melt=", dh2Melt, &
            & "melt excessMeltEn=", excessMeltEn
    end if
#endif ! <- DEBUG_SEAICE

    !-- SnowFall & Evaporation ------------------------------------------------------

    !* snowfall
    if (SnowFall > 0d0) then
       dSnowThick = dSnowThick + SnowFall*dt/DensSnow
       work = 0d0
    else
       ! If the sign of snowfall is negative because of 2nd-order conservative interpolation scheme,
       ! the negative snowfall is treated as the evaporation of snow-layer. 
       work = -SnowFall*dt
    end if
#ifdef DEBUG_SEAICE
    if (debugFlag) then
       write(*,*) "(Surface) SnowFall [kg]:", SnowFall*dt
    end if
#endif ! <- DEBUG_SEAICE

    !* evaporation
    work = work + Evap*dt
    dhsEvap = min(work/DensSnow, hsOld+dSnowThick)
    dSnowThick = dSnowThick - dhsEvap
    
    work = work - DensSnow*dhsEvap
    dh1Evap = min(max(work/DensIce, 0d0), h1Old+dIceThick1)
    dIceThick1 = dIceThick1 - dh1Evap
    if (h1Old + dIceThick1 > 0d0 .and. dh1Evap > 1d-20) then
       hEtmp = (h1Old + dIceThick1 + dh1Evap)*E1 + LFreeze*dh1Evap               ! E1 < 0
       E1 = hEtmp/(h1Old + dIceThick1)
       z_IceTemp(1) = calc_Temp_IceLyr1(E1, SaltSeaIce)
    end if
     
    work = work - DensIce*dh1Evap
    dh2Evap = min( max(work/DensIce, 0d0), h2Old+dIceThick2)
    dIceThick2 = dIceThick2 - dh2Evap
    if (h2Old + dIceThick2 > 0d0 .and. dh2Evap > 1d-20) then
       hEtmp = (h2Old + dIceThick2 + dh2Evap)*E2 + LFreeze*dh2Evap                ! E2 < 0
       E2 = hEtmp/(h2Old + dIceThick2)
       z_IceTemp(2) = calc_Temp_IceLyr2(E2, SaltSeaIce)
    end if
        
    work = work - DensIce*dh2Evap
    if( work > 1d-20) then
       excessMeltEn = excessMeltEn + LFreeze*work
       wice =   wice + work/dt
    end if
    
#ifdef DEBUG_SEAICE
    if (debugFlag) then

       write(*,*) "dhsEvap=", dhsEvap, ", dh1Evap=", dh1Evap, ", dh2Evap=", dh2Evap, &
            & "evap excessMeltEn=", excessMeltEn, LEvap*work, ", evap [kg]:", Evap*dt

       write(*,*) " (SfcResHFlx+BtmResHFlx)*dt=", (SfcResHFlx + BtmResHFlx)*dt
       write(*,*) "Evap+Snow", -LFreeze*dt*(SnowFall - Evap)
       write(*,*) " New - Old =", &
            & - DensSnow*LFreeze*dSnowThick                                        &
            & + DensIce*(E1*(h1Old + dIceThick1) - E1Old*h1Old + (E2*(h2Old + dIceThick2) - E2Old*h2Old))
       write(*,*) "check:", - DensSnow*LFreeze*dSnowThick                          &
            & + DensIce*(E1*(h1Old + dIceThick1) - E1Old*h1Old + (E2*(h2Old + dIceThick2) - E2Old*h2Old))  &
            & - ((SfcResHFlx + BtmResHFlx - LFreeze*(SnowFall - Evap))*dt - excessMeltEn)
       write(*,*) " ExcessMeltEn =", ExcessMeltEn
       
    end if
#endif ! <- DEBUG_SEAICE

       
  end subroutine DSIce_ThermoDyn_Winton2000_CalcLyrMassChange

  !-----------------------------------------------
   
  subroutine DSIce_ThermoDyn_Winton2000_AdjustLyrInternal( &
       & SnowThick, IceThick,                                & ! (out)
       & z_IceTemp, Wice, excessMeltEn,                      & ! (inout)
       & SnowThick0, dSnowThick,                             & ! (in)
       & IceThick0, dIceThick1, dIceThick2,                  & ! (in)
       & dt, debugFlag                                       & ! (in)
       & )

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: SnowThick
    real(DP), intent(out) :: IceThick
    real(DP), intent(inout) :: z_IceTemp(2)
    real(DP), intent(inout) :: Wice
    real(DP), intent(inout) :: excessMeltEn
    real(DP), intent(in) :: SnowThick0
    real(DP), intent(in) :: dSnowThick
    real(DP), intent(in) :: IceThick0
    real(DP), intent(in) :: dIceThick1
    real(DP), intent(in) :: dIceThick2
    real(DP), intent(in) :: dt
    logical, intent(in) :: debugFlag
    
    ! 局所変数
    ! Local variables
    !    
    real(DP) :: h1
    real(DP) :: h2
    real(DP) :: dh1
    real(DP) :: dhs
    real(DP) :: coef

    real(DP) :: T1Tmp
    real(DP) :: T2Tmp
    real(DP) :: dummy
    real(DP) :: extraSHEn ! extra sensible heat energy
    real(DP) :: SIceTotEn
    real(DP) :: dIceMelt
    
    real(DP) :: IceMeltTemp
    
    ! 実行文; Executable statements
    !
    
    SnowThick = SnowThick0 + dSnowThick
    h1 = 0.5d0*IceThick0 + dIceThick1
    h2 = 0.5d0*IceThick0 + dIceThick2
    IceThick = h1 + h2

#ifdef DEBUG_SEAICE
    if (debugFlag) then
       write(*,*) "AdjustLyr internal before:", &
            & "SIceEn:", calc_SIceTotEn( SnowThick, h1, z_IceTemp(1), h2, z_IceTemp(2) )
!!$       write(*,*) "AdjustSIceLyr IceLyrAdj: TotMass=", (DensSnow*SnowThick + DensIce*IceThick)
    end if
#endif
    
    IceMeltTemp = - Mu*SaltSeaIce
        
    if (SnowThick > 0d0) then
       coef =  (SnowThick - (DensSeaWater - DensIce)/DensSnow*IceThick) &
            &  /DensSeaWater
       dhs = - max(coef*DensIce, 0d0)
       dh1 =   max(coef*DensSnow, 0d0)

       if ( dh1 > 0d0 ) then
          call calc_NewIceTemp( T1Tmp, dummy,                & ! (out)
               & z_IceTemp(1), IceMeltTemp, h1/(h1 + dh1)    & ! (in)
               & )
          z_IceTemp(1) = T1Tmp
          
          SnowThick = SnowThick + dhs
          IceThick = IceThick + dh1
          h1 = h1 + dh1
       end if
    end if

!!$    if (debugFlag) then
!!$       write(*,*) "AdjustSIceLyr Snow2Ice: TotMass=", (DensSnow*SnowThick + DensIce*IceThick)
!!$    end if
    
    if (h1 < h2) then
       call calc_NewIceTemp( T1Tmp, dummy,                    & ! (out)
            & z_IceTemp(1), z_IceTemp(2), h1/(0.5d0*IceThick) & ! (in)
            & )
       z_IceTemp(1) = T1Tmp
    else
       call calc_NewIceTemp( dummy, T2Tmp,                          & ! (out)
            & z_IceTemp(1), z_IceTemp(2), h1/(0.5d0*IceThick) - 1d0 & ! (in)
            & )
       z_IceTemp(2) = T2Tmp

       if ( z_IceTemp(2) > IceMeltTemp ) then
          ! Note: extraSHEn > 0          
          extraSHEn = LFreeze + calc_E_IceLyr2(z_IceTemp(2),SaltSeaIce) ! When the temperature of lower ice layer equals to freezing point, 
                                                                        ! the entalphy is defined to be -LFreeze in Winton(2002). 

!!$          if (debugFlag) then
!!$             write(*,*) "Temp of lower layer > -Mu*S: extraSHEn=",  extraSHEn, "Temp=",  z_IceTemp(1:2)
!!$          end if
          z_IceTemp(2) = IceMeltTemp
          z_IceTemp(1) = calc_Temp_IceLyr1( calc_E_IceLyr1(z_IceTemp(1),SaltSeaIce) + extraSHEn, SaltSeaIce )
          
       end if
    end if

#ifdef DEBUG_SEAICE
    if (debugFlag) then
       write(*,*) "AdjustLyr internal after:", &
            & "SIceEn=", calc_SIceTotEn( SnowThick, 0.5d0*IceThick, z_IceTemp(1), 0.5d0*IceThick, z_IceTemp(2) )
!!$       write(*,*) "AdjustSIceLyr IceLyrAdj: TotMass=", (DensSnow*SnowThick + DensIce*IceThick)
    end if
#endif
    
  contains
    subroutine calc_NewIceTemp( &
         & T1New, T2New,        & ! (out)
         & T1, T2, f1           & ! (in)
         & )

      ! 宣言文; Declaration statement
      !          
      real(DP), intent(out) :: T1New
      real(DP), intent(out) :: T2New
      real(DP), intent(in) :: T1
      real(DP), intent(in) :: T2
      real(DP), intent(in) :: f1

      ! 実行文; Executable statements
      !
      
      T2New =   f1*(T1 + LFreeze/CIce * IceMeltTemp/T1) &
           &  + (1d0 -f1)*T2
      
      T1New = 0.5d0*(T2New - sqrt(T2New**2 - 4d0*IceMeltTemp*LFreeze/CIce))
    end subroutine calc_NewIceTemp
  end subroutine DSIce_ThermoDyn_Winton2000_AdjustLyrInternal
     
  !--------------------------------------------------------------------------------

  elemental function calc_Temp_IceLyr1( &
       & E1, S ) result(T1)

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(in) :: E1, S
    real(DP) :: T1
    real(DP) :: a, b, c

    ! 実行文; Executable statements
    !
    
    a = CIce; b = CIce*Mu*s - LFreeze - E1; c = - LFreeze*Mu*S
    T1 = 0.5d0*(- b - sqrt(b**2 - 4d0*a*c))/a
    
  end function calc_Temp_IceLyr1

  !--------------------------------------------------------------------------------
  
  elemental function calc_Temp_IceLyr2( &
       & E2, S ) result(T2)
    
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: E2, S
    real(DP) :: T2

    ! 実行文; Executable statements
    !
    
    T2 = (E2 + LFreeze)/CIce - Mu*S
    
  end function calc_Temp_IceLyr2

  !--------------------------------------------------------------------------------
  
  elemental function calc_E_IceLyr1(T, S) result(E)
    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: T, S
    real(DP) :: E

    ! 実行文; Executable statements
    !
    
    E =     CIce*(T + Mu*S) &
         & - LFreeze*(1d0 + Mu*S/T)

  end function calc_E_IceLyr1

  !--------------------------------------------------------------------------------
  
  elemental function calc_E_IceLyr2(T, S) result(E2)
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: T, S
    real(DP) :: E2

    ! 実行文; Executable statements
    !
    
    E2 = CIce*(T + Mu*S) - LFreeze
  end function calc_E_IceLyr2

  !--------------------------------------------------------------------------------

  elemental function calc_SIceTotEn( SnowThick, IceThick1, SIceTemp1, IceThick2, SIceTemp2 ) &
       & result( SIceTotEn )
    ! 宣言文; Declaration statement
    !
    real(DP), intent(in) :: SnowThick
    real(DP), intent(in) :: IceThick1
    real(DP), intent(in) :: SIceTemp1
    real(DP), intent(in) :: IceThick2
    real(DP), intent(in) :: SIceTemp2
    real(DP) :: SIceTotEn

    ! 実行文; Executable statements
    !

    SIceTotEn =   DensSnow * SnowThick * LFreeze               &
         &      - DensIce * (                                  &
         &           IceThick1*calc_E_IceLyr1( SIceTemp1, SaltSeaIce )   &
         &         + IceThick2*calc_E_IceLyr2( SIceTemp2, SaltSeaIce )   &
         &        )
    
  end function calc_SIceTotEn
  
end module DSIce_ThermoDyn_Winton2000_mod

