!-------------------------------------------------------------
! Copyright (c) 2015-2017 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief A module to set some variables (e.g., surface and bottom flux) in order to satisfy boundary conditions
!! 
!! @author Yuta Kawai
!!
!!
module DSIce_Boundary_common_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM / SeaIce

  use DOGCM_Admin_Constants_mod, only: &
       & UNDEFVAL,                     &
       & CpOcn => Cp0,                 &
       & DensSeaWater => RefDens,      &
       & LatentHeatVap => LatentHeat

  use UnitConversion_mod, only: &
       & degC2K, K2degC

  use DSIce_Admin_Constants_mod, only: &
       & SBConst, FreezeTempWater,                           &
       & DensFreshWater,                                     &
       & Mu,                                                 &
       & AlbedoOcean, AlbedoIce, AlbedoSnow, AlbedoMeltSnow, &
       & EmissivOcean, EmissivSnow, EmissivIce,              &
       & I0,                                                 &
       & BaseMeltHeatTransCoef,                              &
       & IceMaskMin,                                         &
       & LatentHeatFusion => LFreeze
  
  use DSIce_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  use DSIce_Admin_TInteg_mod, only: &
       & DelTimeOcn

  use DSIce_Admin_Variable_mod, only: &
       & xy_Wice

  use DOGCM_Admin_BC_mod, only: &       
       &  ThermBCTYPE_PresFlux_Han1984Method, ThermBC_Surface

  use DOGCM_Boundary_Vars_mod, only: &
       & xy_SeaSfcTemp0
  
  use DSIce_Boundary_vars_mod, only: &
       & xy_SfcHFlxAI, xy_DSfcHFlxAIDTs, xy_SfcHFlxAO,        &
       & xy_SfcHFlxAI0_ns, xy_SfcHFlxAI0_sr,                      &
       & xy_SfcHFlxAO0,                        &       
       & xy_DelSfcHFlxAI,                                     &
       & xy_SUwRFlx, xy_SDwRFlx, xy_LUwRFlx, xy_LDwRFlx, xy_LatHFlx, xy_SenHFlx, &
       & xy_SnowFall, xy_RainFall, xy_Evap,                   &
       & xy_DLatSenHFlxDTs,                                   &
       & xy_PenSDRFlx,                                        &
       & xy_BtmHFlxIO,                                        &
       & xy_WindStressUAI, xy_WindStressVAI,                  &
       & xy_WindStressUIO, xy_WindStressVIO,                  &
       & xy_FreshWtFlxS, xy_FreshWtFlx,                       &
       & xy_SfcAlbedoAI,                                      &
       & xy_SeaSfcTemp, xy_SeaSfcSalt,                        &
       & xy_SIceSfcTemp0,                                     &
       & xy_SeaSfcU, xy_SeaSfcV,                              &
       & xy_OcnFrzTemp, xy_OcnMixLyrDepth

  use DOGCM_Boundary_Vars_mod, only: &
       & xy_DSfcHFlxAODTs => xy_DSfcHFlxDTs
  
  use DSIce_Boundary_SfcAlbedo_mod, only: &
       & DSIce_Boundary_SfcAlbedo_Init, DSIce_Boundary_SfcAlbedo_Final, &
       & DSIce_Boundary_SfcAlbedo_Get

#include "../DSIceDef.h"

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DSIce_Boundary_common_Init, DSIce_Boundary_common_Final

  public :: DSIce_Boundary_common_UpdateBeforeTstep
  public :: DSIce_Boundary_common_UpdateAfterTstep
  
  public :: DSIce_Boundary_common_CalcSfcHFlx
  public :: DSIce_Boundary_common_CalcBtmHFlx
  
  ! 公開変数
  ! Public variable
  !
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DSIce_Boundary_common_mod' !< Module Name

  integer :: j_dbg

  real(DP), parameter :: IceThickMax_PenSDRFlx = 10d0 ! [m]
  ! This parameter is nessesary to prevent the sea ice temperature from exceeding
  ! the melting point due to penetrative shortwave radiation if the ice layer represented
  ! by two-layer model is too thick. (In such case, however, we should increase the number of ice layers.)
    
contains

  !>
  !!
  !!
  Subroutine DSIce_Boundary_common_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !


!    call read_nmlData(configNmlName)
    call DSIce_Boundary_SfcAlbedo_Init( configNmlName )

#ifdef DEBUG_SEAICE    
    j_dbg = DEBUG_SEAICE_j_dbg 
#else
    j_dbg = JS - 1
#endif
    
  end subroutine DSIce_Boundary_common_Init

  !>
  !!
  !!
  subroutine DSIce_Boundary_common_Final()

    ! 実行文; Executable statements
    !

    call DSIce_Boundary_SfcAlbedo_Final()
    
  end subroutine DSIce_Boundary_common_Final

  !-----------------------------------------

  subroutine DSIce_Boundary_common_UpdateBeforeTstep( &
       & xy_SIceCon, xy_IceThick, xy_SnowThick, xy_SIceSfcTemp, xyz_SIceTemp & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_IceThick(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(in) :: xyz_SIceTemp(IA,JA,KA)

    ! 実行文; Executable statements
    !
 
    !-------------------------------------------------

    !$omp parallel
    !$omp workshare
    
    xy_OcnFrzTemp(:,:) = 273.15d0 - 1.8d0  !- Mu*xy_SeaSfcSalt
    
    where (xy_SIceCon(:,:) >= IceMaskMin)
       xy_FreshWtFlxS = ( &
            &   xy_RainFall                                &
            & + (1d0 - xy_SIceCon)*(xy_SnowFall - xy_Evap) &
            & )/ DensFreshWater
    elsewhere
       xy_FreshWtFlxS = 0d0
    end where

    !$omp end workshare
    !$omp end parallel

    !--------------------------------
    
!!$    call DSIce_Boundary_common_CalcSfcHFlx( &
!!$         & xy_SfcHFlxAI, xy_SfcHFlxAO, xy_PenSDRFlx,                                    & ! (out)
!!$         & xy_DSfcHFlxAIDTs, xy_SfcAlbedoAI,                                            & ! (out)
!!$         & xy_SUwRFlx, xy_SDwRFlx, xy_LUwRFlx, xy_LDwRFlx, xy_LatHFlx, xy_SenHFlx,      & ! (in)
!!$         & xy_DLatSenHFlxDTs,                                                           & ! (in)
!!$         & xy_SIceCon, xy_IceThick, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp,        & ! (in)
!!$         & xy_SfcAlbedoAI                                                               & ! (in)
!!$         & )
    call DSIce_Boundary_common_CalcSfcHFlx_v2( &
       & xy_SfcHFlxAI, xy_SfcHFlxAO, xy_PenSDRFlx,                               & ! (out)
       & xy_DSfcHFlxAIDTs,                                                       & ! (out)
       & xy_SIceCon, xy_IceThick, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp,   & ! (in)
       & xy_SfcAlbedoAI                                                          & ! (in)
       & )

    !--------------------------------
    
    call DSIce_Boundary_common_CalcBtmHFlx( &
         & xy_BtmHFlxIO,                                        & ! (out)
         & xy_SIceCon, xy_SeaSfcU, xy_SeaSfcV, xy_SeaSfcTemp,   & ! (in)
         & xy_OcnFrzTemp, xy_OcnMixLyrDepth,                    & ! (in)
         & DelTimeOcn                                           & ! (in)
         & )
      
    !--------------------------------
    
  end subroutine DSIce_Boundary_common_UpdateBeforeTstep


  subroutine DSIce_Boundary_common_UpdateAfterTstep( &
       & xy_SIceCon, xy_IceThick, xy_SnowThick, xy_SIceSfcTemp, xyz_SIceTemp & ! (in)
       & )
 
    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_IceThick(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)
    real(DP), intent(in) :: xyz_SIceTemp(IA,JA,KA)


    ! 実行文; Executable statements
    !
 
    !----------------------

    xy_FreshWtFlxS(:,:) =   xy_FreshWtFlxS                              &
         &                - xy_Wice / DensFreshWater

    !----------------------

    call DSIce_Boundary_SfcAlbedo_Get( &
         & xy_SfcAlbedoAI,                                                      & ! (out)
         & xy_SIceCon, xy_SnowThick, xy_IceThick, xy_SIceSfcTemp, xy_SeaSfcTemp & ! (in)
         & )
 
    !----------------------
    
  end subroutine DSIce_Boundary_common_UpdateAfterTstep
  
  !-----------------------------------------  
  
  subroutine DSIce_Boundary_common_CalcSfcHFlx( &
       & xy_SfcHFlxAI, xy_SfcHFlxAO, xy_PenSWRFlxAI,                             & ! (out)
       & xy_DSfcHFlxDTsAI, xy_SfcAlbedoAI,                                       & ! (out)
       & xy_SUWRFlx, xy_SDWRFlx, xy_LUwRFlx, xy_LDWRFlx, xy_LatHFlx, xy_SenHFlx, & ! (in)
       & xy_DLatSensHFlxDTs,                                                     & ! (in)
       & xy_SIceCon, xy_IceThick, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp,   & ! (in)
       & xy_SIceAlbedoAI                                                         & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SfcHFlxAI(IA,JA)
    real(DP), intent(out) :: xy_SfcHFlxAO(IA,JA)
    real(DP), intent(out) :: xy_PenSWRFlxAI(IA,JA)
    real(DP), intent(out) :: xy_DSfcHFlxDTsAI(IA,JA)
    real(DP), intent(out) :: xy_SfcAlbedoAI(IA,JA)
    real(DP), intent(inout) :: xy_SUWRFlx(IA,JA)
    real(DP), intent(in) :: xy_SDWRFlx(IA,JA)
    real(DP), intent(inout) :: xy_LUWRFlx(IA,JA)
    real(DP), intent(in) :: xy_LDWRFlx(IA,JA)
    real(DP), intent(in) :: xy_LatHFlx(IA,JA)
    real(DP), intent(in) :: xy_SenHFlx(IA,JA)
    real(DP), intent(in) :: xy_DLatSensHFlxDTs(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_IceThick(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)  !< degC
    real(DP), intent(in) :: xy_SeaSfcTemp(IA,JA)   !< K
    real(DP), intent(in) :: xy_SIceAlbedoAI(IA,JA)
    
    ! 局所変数
    ! Local variables
    !    
    integer :: i
    integer :: j
    real(DP) :: emisiv
    real(DP) :: SIceSfcTempK
    real(DP) :: SIceSfcAlbedo
    real(DP) :: SIceCellAvgTempK
    real(DP) :: LUwRFlxCellAvg
    real(DP) :: DLUwRFlxDTs
    real(DP) :: SLR
    real(DP) :: SSR
    real(DP) :: xy_Check(IA,JA)
    
    real(DP), parameter :: EPS = 1d-20
    
    ! 実行文; Executable statements
    !
      
    !$omp parallel do private( emisiv, SIceSfcTempK, SIceSfcAlbedo, SIceCellAvgTempK, &
    !$omp                     LUwRFlxCellAvg, DLUwRFlxDTs, SLR, SSR &
    !$omp ) collapse(2)
    do j = JS, JE
    do i = IS, IE
          
       if ( xy_SIceCon(i,j) >= IceMaskMin ) then
 
          SIceSfcAlbedo = xy_SIceAlbedoAI(i,j)
          SIceSfcTempK = degC2K( xy_SIceSfcTemp(i,j) )
          SIceCellAvgTempK =   (1d0 - xy_SIceCon(i,j))*xy_SeaSfcTemp(i,j) &
               &             + xy_SIceCon(i,j)*SIceSfcTempK
          
          if ( xy_SnowThick(i,j) <= EPS ) then
             emisiv = EmissivIce
          else
             emisiv = EmissivSnow
          end if
  
          if (ThermBC_Surface == ThermBCTYPE_PresFlux_Han1984Method) then
             LUwRFlxCellAvg = SBConst*SIceCellAvgTempK**4
             DLUwRFlxDTs = 4d0*SBConst*SIceCellAvgTempK**3
             SLR = emisiv*( &
                  &   LUwRFlxCellAvg + DLUwRFlxDTs*(SIceSfcTempK - SIceCellAvgTempK)                   &
                  &    + 4d0*SBConst*xy_SIceSfcTemp0(i,j)**3*(SIceCellAvgTempK - xy_SIceSfcTemp0(i,j)) &
                  & - xy_LDWRFlx(i,j) )
             SSR = - (1d0 - SIceSfcAlbedo)*xy_SDWRFlx(i,j)
          else
             LUwRFlxCellAvg = xy_LUwRFlx(i,j)! SBConst*SIceCellAvgTempK**4
             xy_Check(i,j) = LUwRFlxCellAvg
             DLUwRFlxDTs = 4d0*SBConst*SIceCellAvgTempK**3
             SLR = emisiv*( &
                  &   LUwRFlxCellAvg + DLUwRFlxDTs*(SIceSfcTempK - SIceCellAvgTempK)                   &
                  & - xy_LDWRFlx(i,j) )
!!$             SSR = - (1d0 - SIceSfcAlbedo)*xy_SDWRFlx(i,j)
             SSR = xy_SUWRFlx(i,j) - xy_SDWRFlx(i,j) 
          end if

          if ( xy_SnowThick(i,j) <= EPS .and. xy_IceThick(i,j) < IceThickMax_PenSDRFlx) then
             xy_PenSWRFlxAI(i,j) = I0*SSR
!!$                xy_PenSWRFlxAI(i,j) = - I0*(1d0 - SIceSfcAlbedo)*xy_SDWRFlx(i,j) ! Note: the sign is positive.
          else
             xy_PenSWRFlxAI(i,j) = 0d0
          end if
                    
          xy_SfcHFlxAI(i,j) = &
               &   xy_SenHFlx(i,j)                                 &
               & + (LatentHeatVap + LatentHeatFusion)*xy_Evap(i,j) &
               & + SSR + SLR + xy_PenSWRFlxAI(i,j)                             

#ifdef DEBUG_SEAICE
          if(i==IS .and. j==j_dbg)then
             write(*,*) "SfcHFlxAI(LatH,SensH,LDwR,SDwR,LUwR):",         &
                  & xy_LatHFlx(i,j), xy_SenHFlx(i,j), xy_LDwRFlx(i,j),   &
                  & xy_SDwRFlx(i,j),  emisiv*(SBConst*SIceSfcTempK**4),  &
                  & "SfcHFlxAI0", xy_SfcHFlxAI(i,j),                     &
                  & "SIceSfcTemp:", SIceSfcTempK, "SfcAlbedo:", SIceSfcAlbedo
             write(*,*) "SUwR", xy_SUwRFlx(i,j), "LUwR", xy_LUwRFlx(i,j)
          end if
#endif              
           
          if (ThermBC_Surface == ThermBCTYPE_PresFlux_Han1984Method) then
             ! When we consider the change of latent and sensible heat flux due to
             ! surface temperature in the manner based on Han (1984), we assume that
             ! xy_SIceSfcTemp0 is set to a basic surface temperature (units: K)
             ! at each  grid point.  
          
             xy_SfcHFlxAI(i,j) = xy_SfcHFlxAI(i,j) + &
                  & xy_DLatSenHFlxDTs(i,j)*(SIceCellAvgTempK - xy_SIceSfcTemp0(i,j))
!                  & xy_DLatSenHFlxDTs(i,j)*(SIceSfcTempK - xy_SIceSfcTemp0(i,j))
          end if
           
#ifdef DEBUG_SEAICE
          if(i==IS .and. j==j_dbg)then
             write(*,*) "SfcHFlxModAI:", &
                  & xy_DLatSenHFlxDTs(i,j)*(SIceSfcTempK - xy_SIceSfcTemp0(i,j)), &
                  & SIceSfcTempK, xy_SIceSfcTemp0(i,j)
          end if
#endif             
          xy_DSfcHFlxDTsAI(i,j) = &
               &   xy_DLatSensHFlxDTs(i,j)                      &
               & + emisiv*DLUwRFlxDTs

              
       else ! grid cell not covered by sea ice
          SIceSfcAlbedo       = UNDEFVAL
          xy_SfcHFlxAI(i,j)   = UNDEFVAL
       end if


       if (ThermBC_Surface == ThermBCTYPE_PresFlux_Han1984Method) then
          SLR = EmissivOcean*( &
               &   LUwRFlxCellAvg + DLUwRFlxDTs*(xy_SeaSfcTemp(i,j) - SIceCellAvgTempK)              &
               &    + 4d0*SBConst*xy_SIceSfcTemp0(i,j)**3*(SIceCellAvgTempK - xy_SIceSfcTemp0(i,j))  &
               & - xy_LDWRFlx(i,j) )
          SSR = - (1d0 - AlbedoOcean)*xy_SDWRFlx(i,j)
       else
          SLR = emisiv*( &
               &   LUwRFlxCellAvg + DLUwRFlxDTs*(SIceSfcTempK - SIceCellAvgTempK)                   &
               & - xy_LDWRFlx(i,j) )
          SSR = - (1d0 - AlbedoOcean)*xy_SDWRFlx(i,j)
       end if
       
       xy_SfcHFlxAO(i,j) = &
            &   xy_LatHFlx(i,j) + xy_SenHFlx(i,j)            &
            & + SLR + SSR

       if (ThermBC_Surface == ThermBCTYPE_PresFlux_Han1984Method) then          
          xy_SfcHFlxAO(i,j) = xy_SfcHFlxAO(i,j) + &
               &  xy_DLatSenHFlxDTs(i,j)*(SIceCellAvgTempK - xy_SIceSfcTemp0(i,j))
!!$               &  xy_DLatSenHFlxDTs(i,j)*(xy_SeaSfcTemp(i,j) - xy_SeaSfcTemp0(i,j))
       end if
         
#ifdef DEBUG_SEAICE
       if(i==IS .and. j==j_dbg)then
          write(*,*) "SfcHFlxAO:", xy_SfcHFlxAO(i,j)
       end if
#endif       

    end do
    end do

!!$    write(*,*) "LUwRF:", xy_LUwRFlx(IS,JS:JE)
!!$    write(*,*) "Check:", xy_Check(IS,JS:JE)

  end subroutine DSIce_Boundary_common_CalcSfcHFlx

  subroutine DSIce_Boundary_common_CalcSfcHFlx_v2( &
       & xy_SfcHFlxAI, xy_SfcHFlxAO, xy_PenSWRFlxAI,                             & ! (out)
       & xy_DSfcHFlxDTsAI,                                                       & ! (out)
       & xy_SIceCon, xy_IceThick, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp,   & ! (in)
       & xy_SIceAlbedoAI                                                         & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SfcHFlxAI(IA,JA)
    real(DP), intent(out) :: xy_SfcHFlxAO(IA,JA)
    real(DP), intent(out) :: xy_PenSWRFlxAI(IA,JA)
    real(DP), intent(inout) :: xy_DSfcHFlxDTsAI(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_IceThick(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)  !< degC
    real(DP), intent(in) :: xy_SeaSfcTemp(IA,JA)   !< K
    real(DP), intent(in) :: xy_SIceAlbedoAI(IA,JA)
    
    ! 局所変数
    ! Local variables
    !    
    integer :: i
    integer :: j
    real(DP) :: emisiv
    real(DP) :: SIceSfcTempK
    real(DP) :: SIceSfcAlbedo
    real(DP) :: SIceCellAvgTempK
    real(DP) :: LUwRFlxCellAvg
    real(DP) :: DLUwRFlxDTs
    real(DP) :: SLR
    real(DP) :: SSR
    real(DP) :: xy_Check(IA,JA)
    
    real(DP), parameter :: EPS = 1d-20
    
    ! 実行文; Executable statements
    !
      
    !$omp parallel do private( emisiv, SIceSfcTempK, SIceSfcAlbedo, SIceCellAvgTempK, &
    !$omp                     LUwRFlxCellAvg, DLUwRFlxDTs, SLR, SSR &
    !$omp ) collapse(2)
    do j = JS, JE
    do i = IS, IE
          
       if ( xy_SIceCon(i,j) >= IceMaskMin ) then
 
          SIceSfcAlbedo = xy_SIceAlbedoAI(i,j)
          SIceSfcTempK = degC2K( xy_SIceSfcTemp(i,j) )
          SIceCellAvgTempK =   (1d0 - xy_SIceCon(i,j))*xy_SeaSfcTemp(i,j) &
               &             + xy_SIceCon(i,j)*SIceSfcTempK
          
          if ( xy_SnowThick(i,j) <= EPS ) then
             emisiv = EmissivIce
          else
             emisiv = EmissivSnow
          end if
 
          if ( xy_SnowThick(i,j) <= EPS .and. xy_IceThick(i,j) < IceThickMax_PenSDRFlx) then
             xy_PenSWRFlxAI(i,j) = I0*xy_SfcHFlxAI0_sr(i,j)
          else
             xy_PenSWRFlxAI(i,j) = 0d0
          end if

          xy_SfcHFlxAI(i,j) = &
               &   xy_SfcHFlxAI0_ns(i,j) + xy_SfcHFlxAI0_sr(i,j) &
               & + xy_PenSWRFlxAI(i,j)
            
#ifdef DEBUG_SEAICE
          if(i==IS .and. j==j_dbg)then
             write(*,*) "SfcHFlxAI0", xy_SfcHFlxAI(i,j), "SIceSfcTemp:", xy_SIceSfcTemp(i,j), &
                  & "SfcAlbedo:", SIceSfcAlbedo
          end if
#endif              
           
          if (ThermBC_Surface == ThermBCTYPE_PresFlux_Han1984Method) then
             ! When we consider the change of latent and sensible heat flux due to
             ! surface temperature in the manner based on Han (1984), we assume that
             ! xy_SIceSfcTemp0 is set to a basic surface temperature (units: K)
             ! at each  grid point.  
          
             xy_SfcHFlxAI(i,j) = xy_SfcHFlxAI(i,j) + &
                  & xy_DSfcHFlxAIDTs(i,j)*(SIceSfcTempK - xy_SIceSfcTemp0(i,j))
          end if
                   
#ifdef DEBUG_SEAICE
          if(i==IS .and. j==j_dbg)then
             write(*,*) "SfcHFlxModAI:", &
                  & xy_DSfcHFlxAIDTs(i,j)*(SIceSfcTempK - xy_SIceSfcTemp0(i,j)), &
                  & xy_SIceSfcTemp(i,j), xy_SIceSfcTemp0(i,j)
          end if
#endif             
               
       else ! grid cell not covered by sea ice
          SIceSfcAlbedo       = UNDEFVAL
          xy_SfcHFlxAI(i,j)   = UNDEFVAL
       end if

       xy_SfcHFlxAO(i,j) = xy_SfcHFlxAO0(i,j) 

       if (ThermBC_Surface == ThermBCTYPE_PresFlux_Han1984Method) then          
          xy_SfcHFlxAO(i,j) = xy_SfcHFlxAO(i,j) + &
               &  xy_DSfcHFlxAODTs(i,j)*(xy_SeaSfcTemp(i,j) - xy_SeaSfcTemp0(i,j))
       end if
         
#ifdef DEBUG_SEAICE
       if(i==IS .and. j==j_dbg)then
          write(*,*) "SfcHFlxAO:", xy_SfcHFlxAO(i,j)
       end if
#endif       
 
    end do
    end do
  
  end subroutine DSIce_Boundary_common_CalcSfcHFlx_v2
  
  !-----------------------------------------

  subroutine DSIce_Boundary_common_CalcBtmHFlx( &
       & xy_BtmHFlxIO,                                  & ! (out)
       & xy_SIceCon, xy_UO, xy_VO, xy_SeaSfcTemp,       & ! (in)
       & xy_OcnFrzTemp, xy_OcnMixLyrDepth,              & ! (in)
       & dt                                             & ! (in)
       & )

    use SpmlUtil_mod
    
    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: xy_BtmHFlxIO(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)       
    real(DP), intent(in) :: xy_UO(IA,JA)
    real(DP), intent(in) :: xy_VO(IA,JA)
    real(DP), intent(in) :: xy_SeaSfcTemp(IA,JA)
    real(DP), intent(in) :: xy_OcnMixLyrDepth(IA,JA)
    real(DP), intent(in) :: xy_OcnFrzTemp(IA,JA)
    real(DP), intent(in) :: dt

    ! 局所変数
    ! Local variables
    !        
    integer :: i
    integer :: j
    real(DP) :: FricVel
    real(DP) :: FrzPot
    real(DP), parameter :: DragCoef = 0.00536d0   ! a constant ice-ocean drag coefficient
    
    ! 実行文; Executable statements
    !
  
    !$omp parallel do private(i, FricVel, FrzPot)
    do j = JS, JE
    do i = IS, IE

       if ( xy_SIceCon(i,j) >= IceMaskMin ) then
          FrzPot =    CpOcn*DensSeaWater*xy_OcnMixLyrDepth(i,j)/dt         &
               &    * (xy_OcnFrzTemp(i,j) - xy_SeaSfcTemp(i,j))
          if ( FrzPot <= 0d0 ) then
             FricVel = max(DragCoef*sqrt(xy_UO(i,j)**2 + xy_VO(i,j)**2), 5d-3)
             xy_BtmHFlxIO(i,j) = &
                  & - max( FrzPot ,                                                        &
                  &        BaseMeltHeatTransCoef*FricVel*FrzPot*dt/xy_OcnMixLyrDepth(i,j)  &
                  & )
          else
             xy_BtmHFlxIO(i,j) = - FrzPot
          end if
       else
          xy_BtmHFlxIO(i,j) = 0d0 !UNDEFVAL
       end if
  
#ifdef DEBUG_SEAICE
          if(i==IS .and. j==j_dbg)then
             write(*,*) "BtmHFlxIO:", xy_BtmHFlxIO(i,j), "FrzPot:", FrzPot
          end if
#endif
    end do
    end do
       
  end subroutine DSIce_Boundary_common_CalcBtmHFlx
  
end module DSIce_Boundary_common_mod
