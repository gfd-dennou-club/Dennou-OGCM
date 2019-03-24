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
       & DensSeaWater => RefDens,      &
       & LatentHeatVap => LatentHeat
  
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

  use DOGCM_Admin_Constants_mod, only: &
       & LEvap => LatentHeat
  
  use DSIce_Admin_TInteg_mod, only: &
       & CurrentTime
  
  use DSIce_Boundary_vars_mod, only: &
       & xy_SfcHFlxAI, xy_DSfcHFlxAIDTs, xy_DelSfcHFlxAI,  &
       & xy_PenSDRFlx,                                     &
       & xy_SeaSfcTemp, xy_SeaSfcSalt,                     &
       & xy_OcnFrzTemp, xy_SfcHFlxAO, xy_BtmHFlxIO,        &
       & xy_OcnMixLyrDepth,                                &
       & xy_RainFall, xy_SnowFall, xy_Evap
       
  use DSIce_IO_History_mod, only: &
       & DSIce_IO_History_RegistVar, &
       & DSIce_IO_History_HistPut
    
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
  public :: DSIce_TInt_common_advance_ThermoDyn2
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
       & xy_SIceUA, xy_SIceVA,                                              & ! (out)
       & xy_Wice,                                                           & ! (inout)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,                          & ! (out)
       & xy_SIceSfcTemp0, xyz_SIceEn0,                                      & ! (in)
       & dt                                                                 & ! (in)
       & )


    ! モジュール引用; Use statements
    !    
      
    use DSIce_ThermoDyn_Winton2000_mod, only: &
        & calc_Temp_IceLyr1, calc_Temp_IceLyr2
    
    use DSIce_Dyn_driver_mod, only: &
       & DSIce_Dyn_driver_ADVRHS

    use SpmlUtil_mod
 
     
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xy_SIceConA(IA,JA)
    real(DP), intent(out) :: xy_IceThickA(IA,JA)
    real(DP), intent(out) :: xy_SnowThickA(IA,JA)
    real(DP), intent(inout) :: xy_SIceSfcTempA(IA,JA)
    real(DP), intent(out) :: xyz_SIceTempA(IA,JA,KA)
    real(DP), intent(out) :: xyz_SIceEnA(IA,JA,KA)
    real(DP), intent(out) :: xy_SIceUA(IA,JA)
    real(DP), intent(out) :: xy_SIceVA(IA,JA)
    real(DP), intent(inout) :: xy_Wice(IA,JA)
    real(DP), intent(in) :: xy_SIceCon0(IA,JA)
    real(DP), intent(in) :: xy_IceThick0(IA,JA)
    real(DP), intent(in) :: xy_SnowThick0(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp0(IA,JA)
    real(DP), intent(in) :: xyz_SIceEn0(IA,JA,KA)
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
    real(DP) :: xyz_SIceTemp0(IA,JA,KA)
    
    real(DP) :: IceMass
    real(DP) :: SIceConOld
  
    real(DP) :: xy_Tmp1(IA,JA)
    real(DP) :: xy_Tmp2(IA,JA)
    real(DP) :: xy_Tmp(IA,JA)
    real(DP) :: tmp
    
    ! 実行文; Executable statements
    !


    !$omp parallel do collapse(2) private(IceMass)
    do j = JS, JE
    do i = IS, IE
       IceMass = DensIce*0.5d0*xy_IceThick0(i,j)
       xyz_SIceTemp0(i,j,KS  ) = calc_Temp_IceLyr1(xyz_SIceEn0(i,j,KS  )/IceMass, SaltSeaIce)
       xyz_SIceTemp0(i,j,KS+1) = calc_Temp_IceLyr2(xyz_SIceEn0(i,j,KS+1)/IceMass, SaltSeaIce)
    end do
    end do
    
    call DSIce_Dyn_driver_ADVRHS(                               & 
       & xy_SIceCon_RHS, xy_IceThick_RHS, xy_SnowThick_RHS,     & ! (out)
       & xyz_SIceEn_RHS, xy_SIceSfcTemp_RHS,                    & ! (out)
       & xy_SIceUA, xy_SIceVA,                                  & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,              & ! (in)
       & xyz_SIceEn0                                            & ! (in)
       & )
    
!!$    call DSIce_Dyn_driver_ADVRHS(                               & 
!!$       & xy_SIceCon_RHS, xy_IceThick_RHS, xy_SnowThick_RHS,     & ! (out)
!!$       & xyz_SIceEn_RHS, xy_SIceSfcTemp_RHS,                    & ! (out)
!!$       & xy_SIceUA, xy_SIceVA,                                  & ! (out)
!!$       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,              & ! (in)
!!$       & xyz_SIceEn0                                            & ! (in)
!!$       & )
    
    !$omp parallel do private(IceMass, SIceConOld) collapse(2)
    do j = JS, JE
       do i = IS, IE

          ! * Update the tendency of ice thickness and the entalphy, and snow thickness.  
          !

          SIceConOld = xy_SiceCon0(i,j)
          xy_SIceConA(i,j) = SIceConOld + dt*xy_SIceCon_RHS(i,j)
!!$          if (xy_IceThickA(i,j) + dt*xy_IceThick_RHS(i,j) > 0d0) then
!!$             xy_SIceConA(i,j) = 1d0
!!$          else
!!$             xy_SIceConA(i,j) = 0d0
!!$          end if
          
          if (xy_SIceConA(i,j) > 0d0) then
             xy_IceThickA(i,j) = (SIceConOld*xy_IceThick0(i,j) + dt*xy_IceThick_RHS(i,j))/xy_SIceConA(i,j)
             xy_SnowThickA(i,j) = (SIceConOld*xy_SnowThick0(i,j) + dt*xy_SnowThick_RHS(i,j))/xy_SIceConA(i,j)
             xyz_SIceEnA(i,j,KS:KS+1) = (SIceConOld*xyz_SIceEn0(i,j,KS:KS+1) + dt*xyz_SIceEn_RHS(i,j,KS:KS+1))/xy_SIceConA(i,j)

             IceMass = DensIce*0.5d0*xy_IceThickA(i,j)
             xyz_SIceTempA(i,j,KS  ) = calc_Temp_IceLyr1(xyz_SIceEnA(i,j,KS  )/IceMass, SaltSeaIce)
             xyz_SIceTempA(i,j,KS+1) = calc_Temp_IceLyr2(xyz_SIceEnA(i,j,KS+1)/IceMass, SaltSeaIce)

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

#ifdef DEBUG_SEAICE
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
          end if
#endif
           
       end do
    end do

    !-- Final part of dynamical process -------------------------
    
    call DSIce_TInt_common_MeltThinMinIce( &
       & xy_SnowThickA, xy_IceThickA, xy_SIceConA,     & ! (inout)
       & xy_SIceSfcTempA, xyz_SIceTempA,  xyz_SIceEnA, & ! (inout)
       & xy_BtmHFlxIO, xy_Wice,                        & ! (inout)
       & dt ) 

#ifdef DEBUG_SEAICE
    xy_Tmp1 = xy_SnowFall + xy_RainFall - xy_Evap
    write(*,*) "Top", AvrLonLat_xy(xy_Tmp1(IS:IE,JS:JE))
    
    where(xy_IceThick0 > 0d0)
       xy_Tmp2 =   xy_SIceCon0*(-xy_Wice + xy_RainFall)                                  &
            &    + (1d0 - xy_SIceCon0)*(-xy_Wice + xy_SnowFall + xy_RainFall - xy_Evap)
    elsewhere
       xy_Tmp2 = -xy_Wice + xy_SnowFall + xy_RainFall - xy_Evap
    end where
    write(*,*) "Btm", AvrLonLat_xy(xy_Tmp2(IS:IE,JS:JE))

    xy_Tmp = xy_SIceConA*(DensSnow*xy_SnowThickA + DensIce*xy_IceThickA) &
         & - xy_SIceCon0*(DensSnow*xy_SnowThick0 + DensIce*xy_IceThick0)
!!$    xy_Tmp = DensSnow*(xy_SnowThickA - xy_SnowThick0) + DensIce*(xy_IceThickA - xy_IceThick0)

    write(*,*) "IceMassChenge:", AvrLonLat_xy(xy_Tmp(IS:IE,JS:JE))/dt, &
         & AvrLonLat_xy(xy_Tmp1(IS:IE,JS:JE) - xy_Tmp2(IS:IE,JS:JE))
#endif
 
      
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

    use SpmlUtil_mod, only: AvrLonLat_xy

    
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

!!$    real(DP) :: xy_Tmp(IA,JA)
!!$    real(DP) :: xy_Tmp00(IA,JA)
!!$    real(DP) :: xy_Tmp0(IA,JA)
!!$    real(DP) :: xy_Tmp1(IA,JA)
!!$    real(DP) :: xy_Tmp2(IA,JA)
!!$    real(DP) :: xy_Tmp3(IA,JA)
!!$    real(DP) :: xy_Tmp4(IA,JA)
!!$    real(DP) :: xy_Tmp5(IA,JA)
!!$    real(DP) :: xy_TmpNetHFlx(IA,JA)
    
    ! 実行文; Executable statements
    !
    
!!$    xy_Tmp = 0d0
!!$    xy_Tmp00 = 0d0
!!$    xy_Tmp0 = 0d0
!!$    xy_Tmp1 = 0d0
!!$    xy_Tmp2 = 0d0
!!$    xy_Tmp3 = 0d0
!!$    xy_Tmp4 = 0d0
!!$    xy_Tmp5 = 0d0
     
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
       if (xy_SnowThick0(i,j) <= 1d-20) SfcFrzTemp = - Mu*SaltSeaIce

       xy_DelSfcHFlxAI(i,j) = 0d0
       OcnFrzTemp   = K2degC( xy_OcnFrzTemp(i,j) )
       FrzPot       = UNDEFVAL
       excessMeltEn = 0d0
       
       !---------------------------

       if ( xy_SIceCon0(i,j) >= IceMaskMin ) then
 
#ifdef DEBUG_SEAICE             
          if (i==IS .and. j==j_dbg) then
             write(*,*) "(i,j)=", i, j, "time=", CurrentTime, "*****************************"
             write(*,*) "- Sfc & Btm heat flux ---------"
             write(*,*) "SfcHFlx=", xy_SfcHFlxAI(i,j), "DSfcHFlxDTs=", xy_DSfcHFlxAIDTs(i,j)
             write(*,*) "PenSDRFlx=", xy_PenSDRFlx(i,j), "BtmHFlx0=", xy_BtmHFlxIO(i,j)
          end if
#endif             

          call DSIce_ThermoDyn_Winton2000_CalcIceTemp( &
               & xy_SIceSfcTempA(i,j), z_SIceTempA(KS:KE),            & ! (out)
               & SfcResHFlx, BtmResHFlx,                              & ! (out)
               & xy_SIceSfcTemp0(i,j), z_SIceTemp0(KS:KE),            & ! (in)
               & xy_IceThick0(i,j), xy_SnowThick0(i,j),               & ! (in)
               & xy_SfcHFlxAI(i,j), xy_DSfcHFlxAIDTs(i,j),            & ! (in)
               & xy_PenSDRFlx(i,j), SfcFrzTemp,                       & ! (in)
               & xy_BtmHFlxIO(i,j), OcnFrzTemp,                       & ! (in)
               & dt,  (i==IS .and. j==j_dbg) )
   
          xy_DelSfcHFlxAI(i,j) = xy_DSfcHFlxAIDTs(i,j)*(xy_SIceSfcTempA(i,j) - xy_SIceSfcTemp0(i,j))
           
#ifdef DEBUG_SEAICE             
          if (i==IS .and. j==j_dbg) then
             write(*,*) "-calcTemp----------------------"
             write(*,*) "SfcResHFlx=", SfcResHFlx, "BtmResHFlx=", BtmResHFlx, "DelSfcHFlxAI=", xy_DelSfcHFlxAI(i,j)
             write(*,*) "Temp0=", xy_SIceSfcTemp0(i,j), xyz_SIceTemp0(i,j,KS:KE)
             write(*,*) "TempA=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)
             write(*,*) "SIceEn*-SIceEn0=", &
                  &   calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS+1) ) &
                  & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS+1) )
             write(*,*) "FinalNetEnBudget=", dt*( &
                  &   (xy_SfcHFlxAI(i,j) + xy_DSfcHFlxAIDTs(i,j)*(xy_SIceSfcTempA(i,j) - xy_SIceSfcTemp0(i,j))) &
                  & - xy_BtmHFlxIO(i,j)                                                                         &
                  & - xy_PenSDRFlx(i,j)                                                                         &
                  & )
             write(*,*) "(SfcResHFlx + BtmResHFlx)*dt=", dt*(SfcResHFlx + BtmResHFlx), dt*SfcResHFlx, dt*BtmResHFlx
          end if
#endif
!!$          xy_TmpNetHFlx(i,j) = dt*( &
!!$               &   (xy_SfcHFlxAI(i,j) + xy_DSfcHFlxAIDTs(i,j)*(xy_SIceSfcTempA(i,j) - xy_SIceSfcTemp0(i,j))) &
!!$               & - xy_BtmHFlxIO(i,j)                                                                         &
!!$               & - xy_PenSDRFlx(i,j)                                                                         &
!!$               & )
!!$          xy_Tmp00(i,j) = calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), &
!!$               &                   z_SIceTemp0(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS+1) )
!!$          xy_Tmp0(i,j) = calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), &
!!$               &                   z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS+1) )
            
          !---------------------------
        
          call DSIce_ThermoDyn_Winton2000_CalcLyrMassChange( &
               & dSnowThick, dIceThick1, dIceThick2,                 & ! (out)
               & xy_Wice(i,j), excessMeltEn,                         & ! (out)
               & z_SIceTempA(KS:KE),                                 & ! (inout)
               & xy_SnowThick0(i,j), xy_IceThick0(i,j),              & ! (in)
               & SfcResHFlx, BtmResHFlx, OcnFrzTemp,                 & ! (in)
               & xy_RainFall(i,j), xy_SnowFall(i,j), xy_Evap(i,j),   & ! (in)
               & dt,  (i==IS .and. j==j_dbg) )                         ! (in)

!!$          xy_Tmp1(i,j) = calc_SIceTotEn( (xy_SnowThick0(i,j)+dSnowThick), 0.5d0*xy_IceThick0(i,j)+dIceThick1, &
!!$               &                   z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j)+dIceThick2, z_SIceTempA(KS+1) )
!!$          xy_Tmp2(i,j) = -SfcResHFlx*dt
!!$          xy_Tmp3(i,j) = -BtmResHFlx*dt
!!$          xy_Tmp4(i,j) = +LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j))*dt
!!$          xy_Tmp(i,j)  = +excessMeltEn

#ifdef DEBUG_SEAICE                          
          if (i==IS .and. j==j_dbg) then
             write(*,*) "-calcLyrMassChange----------------------"
             write(*,*) "dhs=", dSnowThick, "dhi1=", dIceThick1, "dhi2=", dIceThick2
             write(*,*) "Temp=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)                
             write(*,*) "excessMeltEn=", excessMeltEn
             write(*,*) "SnowFall=", xy_SnowFall(i,j), "Wice=", xy_Wice(i,j)
             write(*,*) "TotMassCange=", &
                  & DensSnow*dSnowThick + DensIce*(dIceThick1 + dIceThick2), &
                  & "SnowFall+Wice=", (xy_Wice(i,j) + xy_SnowFall(i,j))*dt, "Evap=", xy_Evap(i,j)
             write(*,*) "SIceEn**-SIceEn*=", &
                  &   calc_SIceTotEn( xy_SnowThick0(i,j)+dSnowThick, 0.5d0*xy_IceThick0(i,j)+dIceThick1,     &
                  &                   z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j)+dIceThick2, z_SIceTempA(KS+1) )                                   &
                  & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS+1) ) &
                  & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j))*dt                                                                               &
                  & + excessMeltEn
             write(*,*) "-------------------------"
          end if
#endif

          !---------------------------

          call DSIce_ThermoDyn_Winton2000_AdjustLyrInternal( &
               & SnowThickA, IceThickA,                               & ! (out)
               & z_SIceTempA(KS:KE), xy_Wice(i,j), excessMeltEn,      & ! (inout)
               & xy_SnowThick0(i,j), dSnowThick,                      & ! (in)
               & xy_IceThick0(i,j), dIceThick1, dIceThick2,           & ! (in)
               & dt, (i==IS .and. j==j_dbg) )                           ! (in)


#ifdef DEBUG_SEAICE             
          if (i==IS .and. j==j_dbg) then
             write(*,*) "-adjustLyrInternal ----------------------"
             write(*,*) "hsA=", SnowThickA, "hiA=", IceThickA, "TempA", z_SIceTempA(KS:KE)
             write(*,*) "hs0=", xy_SnowThick0(i,j), "hi0=", xy_IceThick0(i,j)
             write(*,*) "excessMeltEn=", excessMeltEn
             write(*,*) "MassTotA-MassTotB=", &
                  &    DensSnow*SnowThickA         + DensIce*IceThickA           &
                  & - (DensSnow*xy_SnowThick0(i,j) + DensIce*xy_IceThick0(i,j)), &
                  & "SnowFall+Wice=", (xy_Wice(i,j) + (xy_SnowFall(i,j) - xy_Evap(i,j)))*dt
             write(*,*) "SIceEnA-SIceEn0=", &
                  &   calc_SIceTotEn( SnowThickA, 0.5d0*IceThickA, z_SIceTempA(KS), 0.5d0*IceThickA, z_SIceTempA(KS+1) )                         &
                  & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS+1) ) &
                  & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j))*dt                                                                               &
                  & + excessMeltEn
             write(*,*) "-------------------------"
          end if
#endif

          !---------------------------

          SIceConA = 1d0  ! Note: After implmentining fractional sea-ice, this code need to be modified.
          if ( IceThickA <= 0d0 ) then
             SIceConA = 0d0
             xy_SIceSfcTempA(i,j) = UNDEFVAL
             z_SIceTempA(KS:KE)   = UNDEFVAL
          end if

          !---------------------------

          xy_BtmHFlxIO(i,j) = xy_BtmHFlxIO(i,j) - excessMeltEn/dt             
          z_SIceEnA(KS  ) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr1(z_SIceTempA(KS  ),SaltSeaIce)
          z_SIceEnA(KS+1) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr2(z_SIceTempA(KS+1),SaltSeaIce)             
!!$          xy_Tmp(i,j) = excessMeltEn
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

             IceThickA  = FrzPot * dt *                                                    &
                  &       2d0 / ( - DensIce*calc_E_IceLyr1(z_SIceTempA(KS  ),SaltSeaIce)   &
                  &               - DensIce*calc_E_IceLyr2(z_SIceTempA(KS+1),SaltSeaIce) )  

             z_SIceEnA(KS  ) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr1(z_SIceTempA(KS  ),SaltSeaIce)
             z_SIceEnA(KS+1) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr2(z_SIceTempA(KS+1),SaltSeaIce)             

             xy_Wice(i,j) = DensIce*IceThickA/dt
             xy_BtmHFlxIO(i,j) = - FrzPot                                                  &
                  &              + LFreeze * xy_SnowFall(i,j) ! [J/kg * kg/(m2.s)]
          else
             !- No sea ice -------------------------------------------------------
             SIceConA     = 0d0
             SnowThickA   = 0d0
             IceThickA    = 0d0
             z_SIceEnA(:) = 0d0
             z_SIceTempA(:)       = UNDEFVAL
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

          write(*,*) "BtmHFlxIOA=", xy_BtmHFlxIO(i,j), "FreshWtFlxS=", -xy_Wice(i,j)/DensFreshWater
          write(*,*) "***********************************************************************************"
       end if
#endif

    end do
    end do
      
!!$    write(*,*) "SIceEn0:", AvrLonLat_xy(xy_Tmp00(IS:IE,JS:JE))
!!$    write(*,*) "SIceEnDel(CalTemp):", AvrLonLat_xy(xy_Tmp0(IS:IE,JS:JE) - xy_Tmp00(IS:IE,JS:JE))
!!$    write(*,*) "-"
!!$    write(*,*) "SIceEnDel(MassChange):", AvrLonLat_xy(xy_Tmp1(IS:IE,JS:JE) - xy_Tmp0(IS:IE,JS:JE)), &
!!$         & "check:", AvrLonLat_xy(  xy_Tmp2(IS:IE,JS:JE) + xy_Tmp3(IS:IE,JS:JE) + xy_Tmp4(IS:IE,JS:JE)  &
!!$         &                        + xy_Tmp(IS:IE,JS:JE))
!!$    write(*,*) "SfcBtmRes*:", AvrLonLat_xy(xy_Tmp2(IS:IE,JS:JE)+xy_Tmp3(IS:IE,JS:JE))
!!$    write(*,*) "SfcRes*:", AvrLonLat_xy(xy_Tmp2(IS:IE,JS:JE))
!!$    write(*,*) "BtmRes*:", AvrLonLat_xy(xy_Tmp3(IS:IE,JS:JE))
!!$    write(*,*) "NetHDelEn*:", AvrLonLat_xy(xy_TmpNetHFlx(IS:IE,JS:JE))
!!$    write(*,*) "Excess*:", AvrLonLat_xy(xy_Tmp(IS:IE,JS:JE))
!!$    write(*,*) "SnowEvap*:", AvrLonLat_xy(xy_Tmp4(IS:IE,JS:JE))
!!$    write(*,*) "CheckDel:", (xy_Tmp1(IS:IE,JS:JE) - xy_Tmp0(IS:IE,JS:JE)) &
!!$         & - (  xy_Tmp2(IS:IE,JS:JE) + xy_Tmp3(IS:IE,JS:JE) + xy_Tmp4(IS:IE,JS:JE)  &
!!$         &                        + xy_Tmp(IS:IE,JS:JE))
     
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
#endif

  end subroutine DSIce_TInt_common_advance_ThermoDyn

  subroutine DSIce_TInt_common_advance_ThermoDyn2(  &
       & xy_SIceConA, xy_IceThickA, xy_SnowThickA,                & ! (out)
       & xyz_SIceEnA, xy_SIceSfcTempA,                            & ! (out) 
       & xy_Wice,                                                 & ! (out)
       & xy_SIceCon0, xy_IceThick0, xy_SnowThick0,                & ! (out)
       & xy_SIceSfcTemp0, xyz_SIceTemp0, xyz_SIceEn0,             & ! (in)
       & dt                                                       & ! (in)
       & )

    ! モジュール引用; Use statements
    !        
    use DSIce_ThermoDyn_Winton2000_mod, only: &
        & DSIce_ThermoDyn_Winton2000_CalcIceTemp,       &
        & DSIce_ThermoDyn_Winton2000_CalcLyrMassChange, &
        & DSIce_ThermoDyn_Winton2000_AdjustLyrInternal, &
        & calc_E_IceLyr1, calc_E_IceLyr2,               &
        & calc_SIceTotEn

    use SpmlUtil_mod, only: AvrLonLat_xy

    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(out) :: xy_SIceConA(IA,JA)
    real(DP), intent(out) :: xy_IceThickA(IA,JA)
    real(DP), intent(out) :: xy_SnowThickA(IA,JA)
    real(DP), intent(out) :: xyz_SIceEnA(IA,JA,KA)
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
    real(DP) :: SIceCon0
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

    real(DP) :: z_SIceEnNew(KS:KE)
    real(DP) :: IceThickNew
    real(DP) :: IceTempNew
    real(DP) :: Fw
    real(DP) :: Fi
    real(DP) :: Fs
    real(DP) :: z_FiceEn(KS:KE)
    real(DP) :: HFlxAI

    real(DP) :: xy_Tmp1(IA,JA)
    real(DP) :: xy_Tmp2(IA,JA)
    
    ! 実行文; Executable statements
    !
    
    !$omp parallel do private( &
    !$omp   i, SfcResHFlx, BtmResHFlx, SfcFrzTemp, OcnFrzTemp, FrzPot,             &
    !$omp   SIceConA, IceThickA, SnowThickA, z_SIceEnA, z_SIceTempA,               &
    !$omp   z_SIceEn0, z_SIceTemp0, SIceCon0,                                      &
    !$omp   dSnowThick, dIceThick1, dIceThick2, excessMeltEn,                      &
    !$omp   Fw, Fi, Fs,z_FiceEn,z_SIceEnNew, IceThickNew, HFlxAI                   &  
    !$omp ) schedule(guided)
    do j = JS, JE
    do i = IS, IE

       z_SIceEn0(KS:KE) = xyz_SIceEn0(i,j,KS:KE)
       z_SIceTemp0(KS:KE) = xyz_SIceTemp0(i,j,KS:KE)
       SIceCon0 = xy_SIceCon0(i,j)
       
       SfcFrzTemp = 0d0
       if (xy_SnowThick0(i,j) <= 1d-20) SfcFrzTemp = - Mu*SaltSeaIce

       xy_DelSfcHFlxAI(i,j) = 0d0
       OcnFrzTemp   = K2degC( xy_OcnFrzTemp(i,j) )
       FrzPot       = UNDEFVAL
       excessMeltEn = 0d0

       xy_Wice(i,j) = 0d0
       
       !---------------------------

       if ( xy_SIceCon0(i,j) >= IceMaskMin ) then
 
#ifdef DEBUG_SEAICE             
          if (i==IS .and. j==j_dbg) then
             write(*,*) "(i,j)=", i, j, "time=", CurrentTime, "*****************************"
             write(*,*) "- Sfc & Btm heat flux ---------"
             write(*,*) "SfcHFlx=", xy_SfcHFlxAI(i,j), "DSfcHFlxDTs=", xy_DSfcHFlxAIDTs(i,j)
             write(*,*) "PenSDRFlx=", xy_PenSDRFlx(i,j), "BtmHFlx0=", xy_BtmHFlxIO(i,j)
          end if
#endif             
          
          call DSIce_ThermoDyn_Winton2000_CalcIceTemp( &
               & xy_SIceSfcTempA(i,j), z_SIceTempA(KS:KE),            & ! (out)
               & SfcResHFlx, BtmResHFlx,                              & ! (out)
               & xy_SIceSfcTemp0(i,j), z_SIceTemp0(KS:KE),            & ! (in)
               & xy_IceThick0(i,j), xy_SnowThick0(i,j),               & ! (in)
               & xy_SfcHFlxAI(i,j), xy_DSfcHFlxAIDTs(i,j),            & ! (in)
               & xy_PenSDRFlx(i,j), SfcFrzTemp,                       & ! (in)
               & xy_BtmHFlxIO(i,j), OcnFrzTemp,                       & ! (in)
               & dt,  (i==IS .and. j==j_dbg) )
    
          xy_DelSfcHFlxAI(i,j) = xy_DSfcHFlxAIDTs(i,j)*(xy_SIceSfcTempA(i,j) - xy_SIceSfcTemp0(i,j)) 
             
#ifdef DEBUG_SEAICE             
          if (i==IS .and. j==j_dbg) then
             write(*,*) "-calcTemp----------------------"
             write(*,*) "SfcResHFlx=", SfcResHFlx, "BtmResHFlx=", BtmResHFlx, "DelSfcHFlxAI=", xy_DelSfcHFlxAI(i,j)
             write(*,*) "Temp0=", xy_SIceSfcTemp0(i,j), xyz_SIceTemp0(i,j,KS:KE)
             write(*,*) "TempA=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)
             write(*,*) "SIceEn*-SIceEn0=", &
                  &   calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS+1) ) &
                  & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS+1) )
             write(*,*) "FinalNetEnBudget=", dt*( &
                  &   (xy_SfcHFlxAI(i,j) + xy_DSfcHFlxAIDTs(i,j)*(xy_SIceSfcTempA(i,j) - xy_SIceSfcTemp0(i,j))) &
                  & - xy_BtmHFlxIO(i,j)                                                                         &
                  & - xy_PenSDRFlx(i,j)                                                                         &
                  & )
             write(*,*) "(SfcResHFlx + BtmResHFlx)*dt=", dt*(SfcResHFlx + BtmResHFlx), dt*SfcResHFlx, dt*BtmResHFlx
          end if
#endif
           
          !---------------------------
        
          call DSIce_ThermoDyn_Winton2000_CalcLyrMassChange( &
               & dSnowThick, dIceThick1, dIceThick2,                 & ! (out)
               & xy_Wice(i,j), excessMeltEn,                         & ! (out)
               & z_SIceTempA(KS:KE),                                 & ! (inout)
               & xy_SnowThick0(i,j), xy_IceThick0(i,j),              & ! (in)
               & SfcResHFlx, BtmResHFlx, OcnFrzTemp,                 & ! (in)
               & xy_RainFall(i,j), xy_SnowFall(i,j), xy_Evap(i,j),   & ! (in)
               & dt,  (i==IS .and. j==j_dbg) )                         ! (in)

#ifdef DEBUG_SEAICE                          
          if (i==IS .and. j==j_dbg) then
             write(*,*) "-calcLyrMassChange----------------------"
             write(*,*) "dhs=", dSnowThick, "dhi1=", dIceThick1, "dhi2=", dIceThick2
             write(*,*) "Temp=", xy_SIceSfcTempA(i,j), z_SIceTempA(:)                
             write(*,*) "excessMeltEn=", excessMeltEn
             write(*,*) "SnowFall=", xy_SnowFall(i,j), "Wice=", xy_Wice(i,j)
             write(*,*) "TotMassCange=", &
                  & DensSnow*dSnowThick + DensIce*(dIceThick1 + dIceThick2), &
                  & "SnowFall+Wice=", (xy_Wice(i,j) + xy_SnowFall(i,j))*dt, "Evap=", xy_Evap(i,j)
             write(*,*) "SIceEn**-SIceEn*=", &
                  &   calc_SIceTotEn( xy_SnowThick0(i,j)+dSnowThick, 0.5d0*xy_IceThick0(i,j)+dIceThick1,     &
                  &                   z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j)+dIceThick2, z_SIceTempA(KS+1) )                                   &
                  & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTempA(KS+1) ) &
                  & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j))*dt                                                                               &
                  & + excessMeltEn
             write(*,*) "-------------------------"
          end if
#endif

          !---------------------------

          call DSIce_ThermoDyn_Winton2000_AdjustLyrInternal( &
               & SnowThickA, IceThickA,                               & ! (out)
               & z_SIceTempA(KS:KE), xy_Wice(i,j), excessMeltEn,      & ! (inout)
               & xy_SnowThick0(i,j), dSnowThick,                      & ! (in)
               & xy_IceThick0(i,j), dIceThick1, dIceThick2,           & ! (in)
               & dt, (i==IS .and. j==j_dbg) )                           ! (in)


#ifdef DEBUG_SEAICE             
          if (i==IS .and. j==j_dbg) then
             write(*,*) "-adjustLyrInternal ----------------------"
             write(*,*) "hsA=", SnowThickA, "hiA=", IceThickA, "TempA", z_SIceTempA(KS:KE)
             write(*,*) "hs0=", xy_SnowThick0(i,j), "hi0=", xy_IceThick0(i,j)
             write(*,*) "excessMeltEn=", excessMeltEn
             write(*,*) "MassTotA-MassTotB=", &
                  &    DensSnow*SnowThickA         + DensIce*IceThickA           &
                  & - (DensSnow*xy_SnowThick0(i,j) + DensIce*xy_IceThick0(i,j)), &
                  & "SnowFall+Wice=", (xy_Wice(i,j) + (xy_SnowFall(i,j) - xy_Evap(i,j)))*dt
             write(*,*) "SIceEnA-SIceEn0=", &
                  &   calc_SIceTotEn( SnowThickA, 0.5d0*IceThickA, z_SIceTempA(KS), 0.5d0*IceThickA, z_SIceTempA(KS+1) )                         &
                  & - calc_SIceTotEn( xy_SnowThick0(i,j), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS), 0.5d0*xy_IceThick0(i,j), z_SIceTemp0(KS+1) ) &
                  & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j))*dt                                                                               &
                  & + excessMeltEn
             write(*,*) "-------------------------"
          end if
#endif

  
          z_SIceEnA(KS  ) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr1(z_SIceTempA(KS  ),SaltSeaIce)
          z_SIceEnA(KS+1) = 0.5d0*IceThickA*DensIce*calc_E_IceLyr2(z_SIceTempA(KS+1),SaltSeaIce)             
          
          Fi  = (IceThickA - xy_IceThick0(i,j))/dt
          Fs = (SnowThickA - xy_SnowThick0(i,j))/dt
          z_FiceEn(KS:KS+1) = (z_SIceEnA(KS:KS+1) - xyz_SIceEn0(i,j,KS:KS+1))/dt
          
       else
          Fi = 0d0
          Fs = 0d0
          z_FiceEn(:) = 0d0
       end if


       FrzPot =  CpOcn*DensSeaWater*xy_OcnMixLyrDepth(i,j)/dt         &
            &    * ( OcnFrzTemp - K2degC(xy_SeaSfcTemp(i,j)) )            ! [J/(K.kg)*(kg.m-3)*m.s-1.K = J/(m2.s)]


       xy_SfcHFlxAO(i,j) = xy_SfcHFlxAO(i,j)
       if ( (FrzPot + xy_SfcHFlxAO(i,j)) > 0d0 ) then
          !- Initalize sea ice ------------------------------------------------                
!!$          xy_SIceSfcTempA(i,j) = OcnFrzTemp
          IceThickNew  = (FrzPot + xy_SfcHFlxAO(i,j))* dt *                          &
                  &       2d0 / ( - DensIce*calc_E_IceLyr1(OcnFrzTemp, SaltSeaIce)   &
                  &               - DensIce*calc_E_IceLyr2(OcnFrzTemp, SaltSeaIce) )

          z_SIceEnNew(KS  ) = 0.5d0*IceThickNew*DensIce*calc_E_IceLyr1(OcnFrzTemp, SaltSeaIce)
          z_SIceEnNew(KS+1) = 0.5d0*IceThickNew*DensIce*calc_E_IceLyr2(OcnFrzTemp, SaltSeaIce)          
          Fw = IceThickNew / dt
          if (SIceCon0 < IceMaskMin) xy_SIceSfcTempA(i,j) = OcnFrzTemp          
       else 
          Fw = 0d0
          FrzPot = 0d0
          z_SIceEnNew(:) = 0d0
       end if
        
       !------------------------------------------------------------------------------------          
       
       
       SIceConA = (SIceCon0 + dt * Fw/0.2d0) &
            & /( 1d0 + dt*(Fw/0.2d0 - min(0d0,Fi)/(2d0*SIceCon0*(xy_IceThick0(i,j) + Fi*dt) + 1d-12)) )
       xy_BtmHFlxIO(i,j) =   (1d0 - SIceCon0) * (min(xy_SfcHFlxAO(i,j),-FrzPot) + LFreeze*xy_SnowFall(i,j))  & 
            &              + SIceCon0*(xy_BtmHFlxIO(i,j) - excessMeltEn/dt)
       
       if (SIceCon0 < IceMaskMin) then
          xy_SfcHFlxAO(i,j) = LFreeze*xy_SnowFall(i,j)
          xy_BtmHFlxIO(i,j) = xy_SfcHFlxAO(i,j)
       end if
       
       xy_Wice(i,j) = (1d0 - SIceCon0)*DensIce*Fw + SIceCon0*xy_Wice(i,j)

       ! Resize
       if (SIceConA < IceMaskMin) then
       end if
         
       if (SIceConA > 1d-10) then
          IceThickA  = (SIceCon0*xy_IceThick0(i,j)  + dt*((1d0 - SIceCon0)*Fw + SIceCon0*Fi)   )/SIceConA
          SnowThickA = (SIceCon0*xy_SnowThick0(i,j) + dt*(                    + SIceCon0*Fs)   )/SIceConA
          z_SIceEnA(KS:KE) = (SIceCon0*z_SIceEn0(KS:KE) + (1d0 - SIceCon0)*z_SIceEnNew(KS:KE) &
               &                                        + dt*SIceCon0*z_FiceEn(KS:KE)          )/SIceConA
       else
          !- In the case of no sea-ice -------------------------------------------------------
          SIceConA     = 0d0
          SnowThickA   = 0d0
          IceThickA    = 0d0
          z_SIceEnA(:) = 0d0
          z_SIceTempA(:)       = UNDEFVAL
          xy_SIceSfcTempA(i,j) = UNDEFVAL

!!$          if (SIceCon0 < IceMaskMin) then
!!$             xy_Wice(i,j)      = 0d0
!!$             xy_BtmHFlxIO(i,j) = LFreeze * xy_SnowFall(i,j) ! [J/kg * kg/(m2.s)]
!!$          end if
       end if

!!$       write(*,*) "CheckEn:", j, &
!!$            & (  SIceConA*(sum(z_SIceEnA(KS:KE)) - LFreeze*SnowThickA*DensSnow) &
!!$            &  - SIceCon0*(sum(z_SIceEn0(KS:KE)) - LFreeze*xy_SnowThick0(i,j)*DensSnow))/dt, &
!!$            & SIceCon0*(sum(z_FiceEn(KS:KS+1)) - LFreeze*Fs*DensSnow), &
!!$            & SIceCon0*( &
!!$            & - (xy_SfcHFlxAI(i,j) + xy_DelSfcHFlxAI(i,j) - xy_PenSDRFlx(i,j))   &
!!$            & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j)))                       &
!!$            & + xy_BtmHFlxIO(i,j) - (1d0 - SIceCon0)*(xy_SfcHFlxAO(i,j)),        &
!!$            & ":", (1d0 - SIceCon0)*sum(z_SIceEnNew(KS:KE)), (1d0 - SIceCon0)*(-FrzPot)

!!$       if (SIceCon0 >= IceMaskMin) then
!!$          xy_Tmp1(i,j) = SIceCon0*( &
!!$               & - (xy_SfcHFlxAI(i,j) + xy_DelSfcHFlxAI(i,j) - xy_PenSDRFlx(i,j))   &
!!$               & - LFreeze*(xy_SnowFall(i,j) - xy_Evap(i,j)))                       &
!!$               & + xy_BtmHFlxIO(i,j) - (1d0 - SIceCon0)*(xy_SfcHFlxAO(i,j))
!!$       else
!!$          xy_Tmp1(i,j) = 0d0
!!$       end if
       
       !---------------------------------------------------------------------

       ! Calculate tendency due to sea-ice thermodynamics
       xy_SIceConA(i,j) = SIceConA
       xy_IceThickA(i,j) = IceThickA
       xy_SnowThickA(i,j) = SnowThickA
       xyz_SIceEnA(i,j,KS:KE) = z_SIceEnA(KS:KE)
       
#ifdef DEBUG_SEAICE
       if (i==IS .and. j==j_dbg) then
          FrzPot = CpOcn*DensSeaWater*xy_OcnMixLyrDepth(i,j)/dt         &
               &    * ( OcnFrzTemp - K2degC(xy_SeaSfcTemp(i,j)) )            ! [J/(K.kg)*(kg.m-3)*m.s-1.K = J/(m2.s)]
          
          write(*,*) "sice thermodynamics summary: (i,j)=", i, j, "time=", CurrentTime, "-----------------"
          write(*,*) &
               & "SIceConA=",SIceConA, "hiA=", IceThickA, "hsA=", SnowThickA,                                                     &
               & "TsfcA=", xy_SIceSfcTempA(i,j), "TA=", z_SIceTempA(KS:KE), "EnA=", z_SIceEnA(KS:KE),       &
               & "SIceConA=", SIceConA
          write(*,*) &
               & "SIceCon0=", xy_SIceCon0(i,j), "hi0=", xy_IceThick0(i,j), "hs0=", xy_SnowThick0(i,j),                                     &
               & "Tsfc0=", xy_SIceSfcTemp0(i,j), "T0=", xyz_SIceTemp0(i,j,KS:KE), "En0=", z_SIceEn0(KS:KE), &
               & "SIceCon0=", xy_SIceCon0(i,j)

          write(*,*) "SST=", xy_SeaSfcTemp(i,j), "OcnFrzTemp=", OcnFrzTemp, &
               & "FrzPot=", FrzPot, "dhi=", FrzPot/( - DensIce*calc_E_IceLyr2(OcnFrzTemp,SaltSeaIce) )*dt

          write(*,*) "BtmHFlxIOA=", xy_BtmHFlxIO(i,j), "FreshWtFlxS=", -xy_Wice(i,j)/DensFreshWater
          write(*,*) "***********************************************************************************"
       end if
#endif
  
    end do
    end do
           
    !-- Final part of sea-ice thermodynamics process  -----------------------

#ifdef DEBUG_SEAICE
    write(*,*) "----------------------------------"
#endif

  end subroutine DSIce_TInt_common_advance_ThermoDyn2
  
  !----------------------------------------------------------------------
  
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
    real(DP) :: SIceConOld
    real(DP) :: iceVol
    real(DP) :: fac
    
    ! 実行文; Executable statements
    !

    !$omp parallel do private(i, SIceConOld, iceVol, fac)
    do j = JS, JE
    do i = IS, IE

       SIceConOld = xy_SIceCon(i,j)
       if ( SIceConOld > 1d0 ) then
          xy_SIceCon(i,j) = 1d0
          xy_IceThick(i,j) = xy_IceThick(i,j)*SIceConOld
          xy_SnowThick(i,j) = xy_SnowThick(i,j)*SIceConOld
          xyz_SIceEn(i,j,KS:KS+1) = xyz_SIceEn(i,j,KS:KS+1)*SIceConOld
       else if( SIceConOld <= IceMaskMin .and. SIceConOld > 0d0) then
          xy_SIceCon(i,j) = IceMaskMin
          fac = SIceConOld/IceMaskMin
          xy_IceThick(i,j) = fac*xy_IceThick(i,j)
          xy_SnowThick(i,j) = fac*xy_SnowThick(i,j)
          xyz_SIceEn(i,j,KS:KS+1) = fac*xyz_SIceEn(i,j,KS:KS+1)
       end if
          
       if ( (xy_IceThick(i,j) > 0d0 .and. xy_IceThick(i,j) < IceThickMin) ) then
          iceVol = xy_SIceCon(i,j)*xy_IceThick(i,j)
          if (iceVol >= IceMaskMin*IceThickMin) then
             xy_SIceCon(i,j) = iceVol/IceThickMin
             fac = IceThickMin/xy_IceThick(i,j)
             xy_SnowThick(i,j) = fac*xy_SnowThick(i,j)
             xyz_SIceEn(i,j,KS:KS+1) = fac*xyz_SIceEn(i,j,KS:KS+1)
             xy_IceThick(i,j) = IceThickMin
          end if 
       end if 
 
       if (        (xy_IceThick(i,j) > 0d0 .and. xy_IceThick(i,j) < IceThickMin) &  ! <- [i]  the case of too thin ice
            & .or. ( xyz_SIceTemp(i,j,KS+1) >  - Mu*SaltSeaIce                 ) &  ! <- [ii] the case where the temperature of lower ice layer
            & ) then                                                                !         exceeds the freezing point

          ! Melt too thin ice in order to avoid numerical instability
          
          xy_BtmHFlxIO(i,j) = xy_BtmHFlxIO(i,j) + xy_SIceCon(i,j)*( &
               &         - (xyz_SIceEn(i,j,KS) + xyz_SIceEn(i,j,KS+1))        & ! <- The sign is > 0
               &         + DensSnow*xy_SnowThick(i,j)*LFreeze                 & ! <- The sign is > 0
               &         )/dt

          xy_Wice(i,j) =   xy_Wice(i,j) - xy_SIceCon(i,j)*(                     &
               &  DensIce*xy_IceThick(i,j) + DensSnow*xy_SnowThick(i,j)         &
               &  )/dt
 
#ifdef DEBUG_SEAICE
          if (i==IS .and. j==j_dbg) then
             write(*,*) "sice MeltThinIce summary: (i,j)=", i, j, "time=", CurrentTime, "-----------------"
             write(*,*) "SIceEnOri=", xyz_SIceEn(i,j,KS:KE)
             write(*,*) "BtmHFlxIO=", xy_BtmHFlxIO(i,j), "FreshWtFlxS=", -xy_Wice(i,j)/DensFreshWater
             write(*,*) "***********************************************************************************"
          end if
#endif
             
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

