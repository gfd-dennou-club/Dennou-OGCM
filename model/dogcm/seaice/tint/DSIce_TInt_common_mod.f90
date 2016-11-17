!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_TInt_common_mod
  
  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use UnitConversion_mod, only: &
       & degC2K, K2degC

  use DOGCM_Admin_Constants_mod, only: &
       & UNDEFVAL,                     &
       & CpOcn => Cp0,                 &
       & DensSeaWater => RefDens
  
  use DSIce_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM
        
  use DSIce_Admin_Constants_mod, only: &
       & DensFreshWater,               &
       & DensIce, DensSnow,            &
       & Mu, SaltSeaIce,               &
       & LFreeze,                      &
       & IceMaskMin, IceThickMin

  use DSIce_Admin_TInteg_mod, only: &
       & CurrentTime
  
  use DSIce_Boundary_vars_mod, only: &
       & xy_SfcHFlxAI, xy_DSfcHFlxAIDTs,   &
       & xy_PenSDRFlx,                     &
       & xy_SeaSfcTemp, xy_SeaSfcSalt,     &
       & xy_OcnFrzTemp, xy_BtmHFlxIO,      &
       & xy_OcnMixLyrDepth,                &
       & xy_RainFall, xy_SnowFall, xy_Evap
       

    
#include "../DSIceDef.h"

  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DSIce_TInt_common_Init, DSIce_TInt_common_Final

  public :: DSIce_TInt_common_advance_Dyn
  public :: DSIce_TInt_common_advance_ThermoDyn
  public :: DSIce_TInt_common_advance_Phys
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DSIce_TInt_common_mod' !< Module Name
  integer :: j_dbg 

contains

  !>
  !!
  !!
  Subroutine DSIce_TInt_common_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

    
    !    call read_nmlData(configNmlName)

#ifdef DEBUG_SEAICE    
    j_dbg = DEBUG_SEAICE_j_dbg 
#else
    j_dbg = JS - 1
#endif

  end subroutine DSIce_TInt_common_Init

  !>
  !!
  !!
  subroutine DSIce_TInt_common_Final()

    ! 宣言文; Declaration statement
    !
    
    ! 実行文; Executable statements
    !


  end subroutine DSIce_TInt_common_Final

  !-----------------------------------------------------------------------------------

  subroutine DSIce_TInt_common_advance_Phys( &
       & dt                                                      & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(in) :: dt

    ! 実行文; Executable statements
    !
    
    != Physical process ======================================================

    
  end subroutine DSIce_TInt_common_advance_Phys

  !----------------------------------------------------------------------
  
  subroutine DSIce_TInt_common_advance_Dyn(  & 
       & xy_SIceConA, xy_IceThickA, xy_SnowThickA,                          & ! (out)
       & xy_SIceSfcTempA,                                                   & ! (inout)
       & xyz_SIceTempA, xyz_SIceEnA,                                        & ! (out)
       & xy_Wice,                                                           & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,                          & ! (out)
       & xy_SIceSfcTemp0, xyz_SIceTemp0, xyz_SIceEn0,                       & ! (in)
       & xy_SIceCon_thm, xy_IceThick_thm, xy_SnowThick_thm, xyz_SIceEn_thm, & ! (in) 
       & dt                                                                 & ! (in)
       & )


    ! モジュール引用; Use statements
    !    

    use DSIce_ThermoDyn_Winton2000_mod, only: &
        & calc_Temp_IceLyr1, calc_Temp_IceLyr2
    
    use DSIce_Dyn_driver_mod, only: &
       & DSIce_Dyn_driver_ADVRHS
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xy_SIceConA(IA,JA)
    real(DP), intent(out) :: xy_IceThickA(IA,JA)
    real(DP), intent(out) :: xy_SnowThickA(IA,JA)
    real(DP), intent(inout) :: xy_SIceSfcTempA(IA,JA)
    real(DP), intent(out) :: xyz_SIceTempA(IA,JA,KA)
    real(DP), intent(out) :: xyz_SIceEnA(IA,JA,KA)
    real(DP), intent(out) :: xy_Wice(IA,JA)
    real(DP), intent(in) :: xy_SIceCon0(IA,JA)
    real(DP), intent(in) :: xy_IceThick0(IA,JA)
    real(DP), intent(in) :: xy_SnowThick0(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp0(IA,JA)
    real(DP), intent(in) :: xyz_SIceTemp0(IA,JA,KA)
    real(DP), intent(in) :: xyz_SIceEn0(IA,JA,KA)
    real(DP), intent(in) :: xy_SIceCon_thm(IA,JA)
    real(DP), intent(in) :: xy_IceThick_thm(IA,JA)
    real(DP), intent(in) :: xy_SnowThick_thm(IA,JA)
    real(DP), intent(in) :: xyz_SIceEn_thm(IA,JA,KA)    
    real(DP), intent(in) :: dt
    
    
    ! 作業変数
    ! Work variables
    
    integer :: i
    integer :: j

    real(DP) :: xy_SIceCon_RHS(IA,JA)
    real(DP) :: xy_IceThick_RHS(IA,JA)
    real(DP) :: xy_SnowThick_RHS(IA,JA)
    real(DP) :: xyz_SIceEn_RHS(IA,JA,KA)
    real(DP) :: xy_SIceSfcTemp_RHS(IA,JA)
    
    real(DP) :: xy_SIceU(IA,JA)
    real(DP) :: xy_SIceV(IA,JA)
    real(DP) :: IceMass

    real(DP) :: xy_SIceSfcTempADV(IA,JA)
    
    real(DP), parameter :: EPS_ICETHICK = 1d-13

    
    ! 実行文; Executable statements
    !

    where( xy_SIceSfcTempADV > 0d0)
       xy_SIceSfcTempADV(:,:) = 0d0
    elsewhere
       xy_SIceSfcTempADV(:,:) = xy_SIceSfcTemp0
    end where
    
    call DSIce_Dyn_driver_ADVRHS(                               & 
       & xy_SIceCon_RHS, xy_IceThick_RHS, xy_SnowThick_RHS,     & ! (out)
       & xyz_SIceEn_RHS, xy_SIceSfcTemp_RHS,                    & ! (out)
       & xy_SIceU, xy_SIceV,                                    & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,              & ! (in)
       & xyz_SIceEn0, xy_SIceSfcTempADV                         & ! (in)
       & )    
!!$    xy_SIceCon_RHS = 0d0
!!$    xy_IceThick_RHS = 0d0
!!$    xy_SnowThick_RHS = 0d0
!!$    xyz_SIceEn_RHS = 0d0
    
    !$omp parallel do private(IceMass) collapse(2)
    do j = JS, JE
       do i = IS, IE

          ! * Update the tendency of ice thickness and the entalphy, and snow thickness.  
          !

          xy_IceThickA(i,j)  = xy_IceThick0(i,j)                                     &
               & + dt*(xy_IceThick_RHS(i,j) + xy_IceThick_thm(i,j))

          xy_SnowThickA(i,j) = xy_SnowThick0(i,j)                                    &
               & + dt*(xy_SnowThick_RHS(i,j) + xy_SnowThick_thm(i,j))

          xyz_SIceEnA(i,j,KS:KS+1) = xyz_SIceEn0(i,j,KS:KS+1)                        &
               & + dt*(xyz_SIceEn_RHS(i,j,KS:KS+1) + xyz_SIceEn_thm(i,j,KS:KS+1))

          if ( xy_IceThickA(i,j) > 0d0 ) then
             
             IceMass = DensIce*0.5d0*xy_IceThickA(i,j)
             xyz_SIceTempA(i,j,KS  ) = calc_Temp_IceLyr1(xyz_SIceEnA(i,j,KS  )/IceMass, SaltSeaIce)
             xyz_SIceTempA(i,j,KS+1) = calc_Temp_IceLyr2(xyz_SIceEnA(i,j,KS+1)/IceMass, SaltSeaIce)
             xy_SIceConA(i,j)        = 1d0

             ! * Set the sea-ice surface temperature to the freezing point, if  no sea-ice exists at time N,
             !   or sea-ice is melted completely in thermodynamics process, but then the sea-ice is advected
             !   into the grid cell by the horizontal transport.
             !
             if ( xy_SIceSfcTempA(i,j) == UNDEFVAL ) then
                xy_SIceSfcTempA(i,j) = - Mu*SaltSeaIce
             end if
             
          else
             xy_SIceConA(i,j)           = 0d0
             xy_IceThickA(i,j)          = 0d0
             xy_SnowThickA(i,j)         = 0d0
             xyz_SIceEnA(i,j,KS:KS+1)   = 0d0
             xyz_SIceTempA(i,j,KS:KS+1) = UNDEFVAL
             xy_SIceSfcTempA(i,j)       = UNDEFVAL
          end if
          
       end do
    end do

#ifdef DEBUG_SEAICE
    do j = JS, JE
       do i = IS, IE
          if (i==IS .and. j==j_dbg) then
             write(*,*) "sice dynamics summary: (i,j)=", i, j, "time=", CurrentTime, "-----------------"
             write(*,*) &
                  & "hiA=", xy_IceThickA(i,j), "hsA=", xy_SnowThickA(i,j),             &
                  & "TsfcA=", xy_SIceSfcTempA(i,j), "TA=", xyz_SIceTempA(i,j,KS:KE),   &
                  & "SIceConA=", xy_SIceConA(i,j)
             write(*,*) "SIceEnA=", xyz_SIceEnA(i,j,KS:KE)

             write(*,*) &
                  & "hi0=", xy_IceThick0(i,j), "hs0=", xy_SnowThick0(i,j),             &
                  & "Tsfc0=", xy_SIceSfcTemp0(i,j), "T0=", xyz_SIceTemp0(i,j,KS:KE),   &
                  & "SIceCon0=", xy_SIceCon0(i,j)
             write(*,*) "SIceEn0=", xyz_SIceEn0(i,j,KS:KE)
             
             write(*,*) "BtmHFlxIO=", xy_BtmHFlxIO(i,j), "FreshWtFlxS=", -xy_Wice(i,j)/DensFreshWater
             write(*,*) "***********************************************************************************"
             if (isNan(xy_BtmHFlxIO(i,j))) then
                write(*,*) "Detect Nan. Stop!"
                stop
             end if
          end if
       end do
    end do
#endif !<- DEBUG_SEAICE
    
    !-- Final part of dynamical process -------------------------
    
    call DSIce_TInt_common_MeltThinMinIce( &
       & xy_SnowThickA, xy_IceThickA, xy_SIceConA,     & ! (inout)
       & xy_SIceSfcTempA, xyz_SIceTempA,  xyz_SIceEnA, & ! (inout)
       & xy_BtmHFlxIO, xy_Wice,                        & ! (inout)
       & dt )

    
  end subroutine DSIce_TInt_common_advance_Dyn

  !----------------------------------------------------------------------
  
  subroutine DSIce_TInt_common_advance_ThermoDyn(  &
       & xy_SIceCon_RHS, xy_IceThick_RHS, xy_SnowThick_RHS, & ! (out)
       & xyz_SIceEn_RHS, xy_SIceSfcTempA,                   & ! (out) 
       & xy_Wice,                                           & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,          & ! (out)
       & xy_SIceSfcTemp0, xyz_SIceTemp0, xyz_SIceEn0,       & ! (in)
       & dt                                                 & ! (in)
       & )

    ! モジュール引用; Use statements
    !        
    use DSIce_ThermoDyn_Winton2000_mod, only: &
        & DSIce_ThermoDyn_Winton2000_CalcIceTemp,       &
        & DSIce_ThermoDyn_Winton2000_CalcLyrMassChange, &
        & DSIce_ThermoDyn_Winton2000_AdjustLyrInternal, &
        & calc_E_IceLyr1, calc_E_IceLyr2,               &
        & calc_SIceTotEn
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xy_SIceCon_RHS(IA,JA)
    real(DP), intent(out) :: xy_IceThick_RHS(IA,JA)
    real(DP), intent(out) :: xy_SnowThick_RHS(IA,JA)
    real(DP), intent(out) :: xyz_SIceEn_RHS(IA,JA,KA)
    real(DP), intent(out) :: xy_SIceSfcTempA(IA,JA)
    real(DP), intent(out) :: xy_Wice(IA,JA)
    real(DP), intent(in) :: xy_SIceCon0(IA,JA)
    real(DP), intent(in) :: xy_IceThick0(IA,JA)
    real(DP), intent(in) :: xy_SnowThick0(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp0(IA,JA)
    real(DP), intent(in) :: xyz_SIceTemp0(IA,JA,KA)
    real(DP), intent(in) :: xyz_SIceEn0(IA,JA,KA)
    real(DP), intent(in) :: dt
    
    ! 作業変数
    ! Work variables
    
    integer :: i
    integer :: j

    real(DP) :: SIceConA
    real(DP) :: IceThickA    
    real(DP) :: SnowThickA
    real(DP) :: z_SIceEnA(KS:KE)
    real(DP) :: z_SIceTempA(KS:KE)
    real(DP) :: z_SIceEn0(KS:KE)
    real(DP) :: z_SIceTemp0(KS:KE)
    
    real(DP) :: SfcResHFlx
    real(DP) :: BtmResHFlx
    real(DP) :: SfcFrzTemp
    real(DP) :: dSnowThick
    real(DP) :: dIceThick1
    real(DP) :: dIceThick2
    real(DP) :: FrzPot
    real(DP) :: excessMeltEn
    real(DP) :: OcnFrzTemp

    
    ! 実行文; Executable statements
    !

    
    !$omp parallel do private( &
    !$omp   i, SfcResHFlx, BtmResHFlx, SfcFrzTemp, OcnFrzTemp, FrzPot,             &
    !$omp   SIceConA, IceThickA, SnowThickA, z_SIceEnA, z_SIceTempA,               &
    !$omp   z_SIceEn0, z_SIceTemp0,                                                &
    !$omp   dSnowThick, dIceThick1, dIceThick2, excessMeltEn                       &
    !$omp ) schedule(guided)
    do j = JS, JE
       do i = IS, IE

          z_SIceEn0(KS:KE) = xyz_SIceEn0(i,j,KS:KE)
          z_SIceTemp0(KS:KE) = xyz_SIceTemp0(i,j,KS:KE)
          
          SfcFrzTemp = 0d0
          if (xy_SnowThick0(i,j) == 0d0) SfcFrzTemp = - Mu*SaltSeaIce

          OcnFrzTemp = K2degC( xy_OcnFrzTemp(i,j) )
          FrzPot = UNDEFVAL
                
          !---------------------------

          if ( xy_SIceCon0(i,j) >= IceMaskMin ) then

#ifdef DEBUG_SEAICE             
             if (i==IS .and. j==j_dbg) then
                write(*,*) "(i,j)=", i, j, "time=", CurrentTime, "*****************************"
                write(*,*) "- Sfc & Btm heat flux ---------"
                write(*,*) "SfcHFlx=", xy_SfcHFlxAI(i,j), "DSfcHFlxDTs=", xy_DSfcHFlxAIDTs(i,j)
                write(*,*) "PenSDRFlx=", xy_PenSDRFlx(i,j), "BtmHFlx=", xy_BtmHFlxIO(i,j)
             end if
#endif  ! <- DEBUG_SEAICE             

             call DSIce_ThermoDyn_Winton2000_CalcIceTemp( &
                  & xy_SIceSfcTempA(i,j), z_SIceTempA(KS:KE),            & ! (out)
                  & SfcResHFlx, BtmResHFlx,                              & ! (out)
                  & xy_SIceSfcTemp0(i,j), z_SIceTemp0(KS:KE),            & ! (in)
                  & xy_IceThick0(i,j), xy_SnowThick0(i,j),               & ! (in)
                  & xy_SfcHFlxAI(i,j), xy_DSfcHFlxAIDTs(i,j),            & ! (in)
                  & xy_PenSDRFlx(i,j), SfcFrzTemp,                       & ! (in)
                  & xy_BtmHFlxIO(i,j), OcnFrzTemp,                       & ! (in)
                  & dt,  (i==IS .and. j==j_dbg) )

#ifdef DEBUG_SEAICE             
             if (i==IS .and. j==j_dbg) then
                write(*,*) "-calcTemp----------------------"
                write(*,*) "SfcResHFlx=", SfcResHFlx, "BtmResHFlx=", BtmResHFlx
                write(*,*) "Temp0=", xy_SIceSfcTemp0(i,j), xyz_SIceTemp0(i,j,KS:KE)
                write(*,*) "TempA=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)
                write(*,*) "SIceEn*-SIceEn0=", &
                     &   calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS+1) ) &
                     & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS+1) )
                write(*,*) "FinalNetEnBudget=", dt*( &
                     &   (xy_SfcHFlxAI(i,j) + xy_DSfcHFlxAIDTs(i,j)*(xy_SIceSfcTempA(i,j) - xy_SIceSfcTemp0(i,j))) &
                     & - xy_BtmHFlxIO(i,j)                                                                         &
                     & )
                write(*,*) "(SfcResHFlx + BtmResHFlx)*dt=", dt*(SfcResHFlx + BtmResHFlx)
             end if
#endif ! <- DEBUG_SEAICE
             
             !---------------------------
          
             call DSIce_ThermoDyn_Winton2000_CalcLyrMassChange( &
                  & dSnowThick, dIceThick1, dIceThick2,               & ! (out)
                  & xy_Wice(i,j), excessMeltEn,                       & ! (out)
                  & z_SIceTempA(KS:KE),                               & ! (inout)
                  & xy_SnowThick0(i,j), xy_IceThick0(i,j),            & ! (in)
                  & SfcResHFlx, BtmResHFlx, OcnFrzTemp,               & ! (in)
                  & xy_RainFall(i,j), xy_SnowFall(i,j), xy_Evap(i,j), & ! (in)
                  & dt,  (i==IS .and. j==j_dbg) )                       ! (in)


             xy_BtmHFlxIO(i,j) = xy_BtmHFlxIO(i,j) - excessMeltEn/dt
             
#ifdef DEBUG_SEAICE                          
             if (i==IS .and. j==j_dbg) then
                write(*,*) "-calcLyrMassChange----------------------"
                write(*,*) "dhs=", dSnowThick, "dhi1=", dIceThick1, "dhi2=", dIceThick2
                write(*,*) "Temp=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)                
                write(*,*) "excessMeltEn=", excessMeltEn, "BtmHFlxIO=", xy_BtmHFlxIO(i,j)
                write(*,*) "SnowFall=", xy_SnowFall(i,j), "Wice=", xy_Wice(i,j)
                write(*,*) "TotMassCange=", &
                     & DensSnow*dSnowThick + DensIce*(dIceThick1 + dIceThick2), &
                     & "SnowFall+Wice=", (xy_Wice(i,j) + xy_SnowFall(i,j))*dt
                write(*,*) "SIceEn**-SIceEn*=", &
                     &   calc_SIceTotEn( xy_SnowThick0(i,j)+dSnowThick, 0.5d0*xy_IceThick0(i,j)+dIceThick1,     &
                     &                   z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j)+dIceThick2, z_SIceTempA(KS+1) )                                   &
                     & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS+1) ) &
                     & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j))*dt                                                                               &
                     & + excessMeltEn
                write(*,*) "-------------------------"
             end if
#endif ! <- DEBUG_SEAICE
             
             !---------------------------

             call DSIce_ThermoDyn_Winton2000_AdjustLyrInternal( &
                  & SnowThickA, IceThickA,                               & ! (out)
                  & z_SIceTempA(KS:KE), xy_Wice(i,j),                    & ! (inout)
                  & xy_SnowThick0(i,j), dSnowThick,                      & ! (in)
                  & xy_IceThick0(i,j), dIceThick1, dIceThick2,           & ! (in)
                  & dt, (i==IS .and. j==j_dbg) )                           ! (in)

             
#ifdef DEBUG_SEAICE             
             if (i==IS .and. j==j_dbg) then
                write(*,*) "-adjustLyrInternal ----------------------"
                write(*,*) "hsA=", SnowThickA, "hiA=", IceThickA, "TempA", z_SIceTempA(KS:KE)
                write(*,*) "hs0=", xy_SnowThick0(i,j), "hi0=", xy_IceThick0(i,j)
                write(*,*) "MassTotA-MassTotB=", &
                     &    DensSnow*SnowThickA         + DensIce*IceThickA           &
                     & - (DensSnow*xy_SnowThick0(i,j) + DensIce*xy_IceThick0(i,j)), &
                     & "SnowFall+Wice=", (xy_Wice(i,j) + xy_SnowFall(i,j))*dt
                write(*,*) "SIceEnA-SIceEn0=", &
                     &   calc_SIceTotEn( SnowThickA, 0.5d0*IceThickA, z_SIceTempA(KS), 0.5d0*IceThickA, z_SIceTempA(KS+1) )                         &
                     & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS+1) ) &
                     & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j))*dt                                                                               &
                     & + excessMeltEn
                write(*,*) "-------------------------"
             end if
#endif ! <- DEBUG_SEAICE
             
             !---------------------------

             SIceConA = 1d0  ! Note: After implmentining fractional sea-ice, this code need to be modified.
             if ( IceThickA <= 0d0 ) then
                SIceConA = 0d0
                xy_SIceSfcTempA(i,j) = UNDEFVAL
                z_SIceTempA(KS:KE)   = UNDEFVAL
             end if

             !---------------------------
             
             z_SIceEnA(KS  ) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr1(z_SIceTempA(KS  ),SaltSeaIce)
             z_SIceEnA(KS+1) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr2(z_SIceTempA(KS+1),SaltSeaIce)             
             
          else
             !- The case when there is no sea-ice at the grid point -----------------------------------------
             
             FrzPot =    CpOcn*DensSeaWater*xy_OcnMixLyrDepth(i,j)/dt         &
                  &    * ( OcnFrzTemp - K2degC(xy_SeaSfcTemp(i,j)) )            ! [J/(K.kg)*(kg.m-3)*m.s-1.K = J/(m2.s)]
             
             if ( FrzPot > 0d0 ) then
                !- Initalize sea ice ------------------------------------------------                

                z_SIceTempA(KS  )    = OcnFrzTemp
                z_SIceTempA(KS+1)    = OcnFrzTemp
                xy_SIceSfcTempA(i,j) = OcnFrzTemp
                
                SIceConA   = 1d0
                SnowThickA = 0d0
                
                IceThickA  = FrzPot * dt * &
                     &       2d0 / ( - DensIce*calc_E_IceLyr1(z_SIceTempA(KS  ),SaltSeaIce)   &
                     &               - DensIce*calc_E_IceLyr2(z_SIceTempA(KS+1),SaltSeaIce) )  

                z_SIceEnA(KS  ) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr1(z_SIceTempA(KS  ),SaltSeaIce)
                z_SIceEnA(KS+1) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr2(z_SIceTempA(KS+1),SaltSeaIce)             
                
                xy_Wice(i,j) = DensIce*IceThickA/dt
                xy_BtmHFlxIO(i,j) = - FrzPot &
                     &              + LFreeze * xy_SnowFall(i,j) ! [J/kg * kg/(m2.s)]
             else
                !- No sea ice -------------------------------------------------------
                SIceConA     = 0d0
                SnowThickA   = 0d0
                IceThickA    = 0d0
                z_SIceEnA(:) = 0d0
                xy_SIceSfcTempA(i,j) = UNDEFVAL

                xy_Wice(i,j)      = 0d0
                xy_BtmHFlxIO(i,j) = LFreeze * xy_SnowFall(i,j) ! [J/kg * kg/(m2.s)]
             end if

          end if
             
          !---------------------------------------------------------------------
                          
          ! Calculate tendency due to sea-ice thermodynamics
          xy_SIceCon_RHS(i,j)   = (SIceConA - xy_SIceCon0(i,j))/dt
          xy_IceThick_RHS(i,j)  = (IceThickA - xy_IceThick0(i,j))/dt
          xy_SnowThick_RHS(i,j) = (SnowThickA - xy_SnowThick0(i,j))/dt
          xyz_SIceEn_RHS(i,j,KS:KS+1) = (z_SIceEnA(KS:KS+1) - z_SIceEn0(KS:KS+1))/dt
             

#ifdef DEBUG_SEAICE
          if (i==IS .and. j==j_dbg) then
             write(*,*) "sice thermodynamics summary: (i,j)=", i, j, "time=", CurrentTime, "-----------------"
             write(*,*) &
                  & "hiA=", IceThickA, "hsA=", SnowThickA,                                                     &
                  & "TsfcA=", xy_SIceSfcTempA(i,j), "TA=", z_SIceTempA(KS:KE), "EnA=", z_SIceEnA(KS:KE),       &
                  & "SIceConA=", SIceConA
             write(*,*) &
                  & "hi0=", xy_IceThick0(i,j), "hs0=", xy_SnowThick0(i,j),                                     &
                  & "Tsfc0=", xy_SIceSfcTemp0(i,j), "T0=", xyz_SIceTemp0(i,j,KS:KE), "En0=", z_SIceEn0(KS:KE), &
                  & "SIceCon0=", xy_SIceCon0(i,j)
             
             write(*,*) "SST=", xy_SeaSfcTemp(i,j), "OcnFrzTemp=", OcnFrzTemp, &
                  & "FrzPot=", FrzPot, "dhi=", FrzPot/( - DensIce*calc_E_IceLyr2(OcnFrzTemp,SaltSeaIce) )*dt

             write(*,*) "BtmHFlxIO=", xy_BtmHFlxIO(i,j), "FreshWtFlxS=", -xy_Wice(i,j)/DensFreshWater
             write(*,*) "***********************************************************************************"
          end if
#endif !<- DEBUG_SEAICE
          
          
       end do
    end do

    !-- Final part of sea-ice thermodynamics process  -----------------------

!!$
!!$
!!$    call DSIce_TInt_common_MeltThinMinIce( &
!!$       & xy_SnowThickA, xy_IceThickA, xy_SIceConA,     & ! (inout)
!!$       & xy_SIceSfcTempA, xyz_SIceTempA,  xyz_SIceEnA, & ! (inout)
!!$       & xy_BtmHFlxIO, xy_Wice,                        & ! (inout)
!!$       & dt )
    

          
#ifdef DEBUG_SEAICE
    write(*,*) "----------------------------------"
#endif !<- DEBUG_SEAICE
    
  end subroutine DSIce_TInt_common_advance_ThermoDyn

  subroutine DSIce_TInt_common_MeltThinMinIce( &
       & xy_SnowThick, xy_IceThick, xy_SIceCon,     & ! (inout)
       & xy_SIceSfcTemp, xyz_SIceTemp,  xyz_SIceEn, & ! (inout)
       & xy_BtmHFlxIO, xy_Wice,                     & ! (inout)
       & dt )                                         ! (in)

    ! モジュール引用; Use statements
    !        
    use DSIce_ThermoDyn_Winton2000_mod, only: &
        & calc_E_IceLyr1, calc_E_IceLyr2,               &
        & calc_SIceTotEn
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xy_SnowThick(IA,JA)
    real(DP), intent(inout) :: xy_IceThick(IA,JA)
    real(DP), intent(inout) :: xy_SIceCon(IA,JA)
    real(DP), intent(inout) :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(inout) :: xyz_SIceTemp(IA,JA,KA)
    real(DP), intent(inout) :: xyz_SIceEn(IA,JA,KA)
    real(DP), intent(inout) :: xy_BtmHFlxIO(IA,JA)
    real(DP), intent(inout) :: xy_Wice(IA,JA)
    real(DP), intent(in) :: dt
    
    ! 作業変数
    ! Work variables
    
    integer :: i
    integer :: j


    ! 実行文; Executable statements
    !
    
    !$omp parallel do private(i)
    do j = JS, JE
       do i = IS, IE
          
          if (        (xy_IceThick(i,j) > 0d0 .and. xy_IceThick(i,j) < IceThickMin) &  ! <- [i]  the case of too thin ice
               & .or. ( xyz_SIceTemp(i,j,KS+1) >  - Mu*SaltSeaIce                 ) &  ! <- [ii] the case where the temperature of lower ice layer
               & ) then                                                                !         exceeds the freezing point

             ! Melt too thin ice in order to avoid numerical instability
             
             xy_BtmHFlxIO(i,j) = xy_BtmHFlxIO(i,j) + ( &
                  &         - (xyz_SIceEn(i,j,KS) + xyz_SIceEn(i,j,KS+1))        & ! <- The sign is > 0
                  &         + DensSnow*xy_SnowThick(i,j)*LFreeze                 & ! <- The sign is > 0
                  &         )/dt

             xy_Wice(i,j) =   xy_Wice(i,j) - (                                     &
                  &  DensIce*xy_IceThick(i,j) + DensSnow*xy_SnowThick(i,j)         &
                  &  )/dt

#ifdef DEBUG_SEAICE
          if (i==IS .and. j==j_dbg) then
             write(*,*) "sice MeltThinIce summary: (i,j)=", i, j, "time=", CurrentTime, "-----------------"
             write(*,*) "SIceEnOri=", xyz_SIceEn(i,j,KS:KE)
             write(*,*) "BtmHFlxIO=", xy_BtmHFlxIO(i,j), "FreshWtFlxS=", -xy_Wice(i,j)/DensFreshWater
             write(*,*) "***********************************************************************************"
          end if
#endif !<- DEBUG_SEAICE
             
             xy_SIceCon(i,j)         = 0d0
             xy_IceThick(i,j)        = 0d0
             xy_SnowThick(i,j)       = 0d0
             xyz_SIceEn(i,j,KS:KE)   = 0d0 
             xyz_SIceTemp(i,j,KS:KE) = UNDEFVAL
             xy_SIceSfcTemp(i,j)     = UNDEFVAL

             
          end if
          
       end do
    end do
    
  end subroutine DSIce_TInt_common_MeltThinMinIce
  
end module DSIce_TInt_common_mod

