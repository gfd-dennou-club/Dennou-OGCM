!-------------------------------------------------------------
! Copyright (c) 2015-2016 Yuta Kawai. All rights reserved.
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
       & DensSeaWater => RefDens

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
       & IceMaskMin
  
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
  
  use DSIce_Boundary_vars_mod, only: &
       & xy_SfcHFlxAI, xy_DSfcHFlxAIDTs, xy_SfcHFlxAO,        &
       & xy_SDwRFlx, xy_LDwRFlx, xy_LatHFlx, xy_SenHFlx,      &
       & xy_RainFall,                                         &
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

  public :: DSIce_Boundary_common_CalcSfcAlbedo
  
  ! 公開変数
  ! Public variable
  !
  
  ! 非公開変数
  ! Private variable
  !
  
  character(*), parameter:: module_name = 'DSIce_Boundary_common_mod' !< Module Name

  
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

    
  end subroutine DSIce_Boundary_common_Init

  !>
  !!
  !!
  subroutine DSIce_Boundary_common_Final()

    ! 実行文; Executable statements
    !

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
    
    where (xy_SIceCon(:,:) > IceMaskMin)
       xy_FreshWtFlxS = xy_RainFall / DensFreshWater
    elsewhere
       xy_FreshWtFlxS = 0d0
    end where

    !$omp end workshare
    !$omp end parallel

    !--------------------------------
    
    call DSIce_Boundary_common_CalcSfcHFlx( &
         & xy_SfcHFlxAI, xy_SfcHFlxAO, xy_PenSDRFlx,                       & ! (out)
         & xy_DSfcHFlxAIDTs, xy_SfcAlbedoAI,                               & ! (out)
         & xy_SDwRFlx, xy_LDwRFlx, xy_LatHFlx, xy_SenHFlx,                 & ! (in)
         & xy_DLatSenHFlxDTs,                                              & ! (in)
         & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp         & ! (in)
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

    xy_FreshWtFlxS(:,:) =   xy_FreshWtFlxS            &
         &                - xy_Wice / DensFreshWater

    !----------------------
    
    call DSIce_Boundary_common_CalcSfcAlbedo( &
       & xy_SfcAlbedoAI,                          & ! (out)
       & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp & ! (in)
       & )    
    
  end subroutine DSIce_Boundary_common_UpdateAfterTstep
  
  !-----------------------------------------  

  subroutine DSIce_Boundary_common_CalcSfcAlbedo( &
       & xy_SIceAlbedoAI,                         & ! (out)
       & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SIceAlbedoAI(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)

    ! 局所変数
    ! Local variables
    !        
    integer :: i
    integer :: j
    real(DP) :: SIceSfcTempK

    ! 実行文; Executable statements
    !
    
    !$omp parallel do private(SIceSfcTempK) collapse(2)
    do j = JS, JE
       do i = IS, IE
          if ( xy_SIceCon(i,j) > IceMaskMin ) then
             SIceSfcTempK = degC2K( xy_SIceSfcTemp(i,j) )             
             if ( xy_SnowThick(i,j) <= 0d0 ) then
                xy_SfcAlbedoAI(i,j) = AlbedoIce
             elseif( xy_SIceSfcTemp(i,j) >= FreezeTempWater ) then
                xy_SfcAlbedoAI(i,j) = AlbedoMeltSnow
             else
                xy_SfcAlbedoAI(i,j) = AlbedoSnow
             end if
          else ! grid cell not covered by sea ice
             xy_SfcAlbedoAI(i,j) = UNDEFVAL
          end if
       end do
    end do
    
  end subroutine DSIce_Boundary_common_CalcSfcAlbedo

  !-----------------------------------------  
  
  subroutine DSIce_Boundary_common_CalcSfcHFlx( &
       & xy_SfcHFlxAI, xy_SfcHFlxAO, xy_PenSWRFlxAI,               & ! (out)
       & xy_DSfcHFlxDTsAI, xy_SfcAlbedoAI,                         & ! (out)
       & xy_SDWRFlx, xy_LDWRFlx, xy_LatHFlx, xy_SenHFlx,           & ! (in)
       & xy_DLatSensHFlxDTs,                                       & ! (in)
       & xy_SIceCon, xy_SnowThick, xy_SIceSfcTemp, xy_SeaSfcTemp   & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xy_SfcHFlxAI(IA,JA)
    real(DP), intent(out) :: xy_SfcHFlxAO(IA,JA)
    real(DP), intent(out) :: xy_PenSWRFlxAI(IA,JA)
    real(DP), intent(out) :: xy_DSfcHFlxDTsAI(IA,JA)
    real(DP), intent(out) :: xy_SfcAlbedoAI(IA,JA)
    real(DP), intent(in) :: xy_SDWRFlx(IA,JA)
    real(DP), intent(in) :: xy_LDWRFlx(IA,JA)
    real(DP), intent(in) :: xy_LatHFlx(IA,JA)
    real(DP), intent(in) :: xy_SenHFlx(IA,JA)
    real(DP), intent(in) :: xy_DLatSensHFlxDTs(IA,JA)
    real(DP), intent(in) :: xy_SIceCon(IA,JA)
    real(DP), intent(in) :: xy_SnowThick(IA,JA)
    real(DP), intent(in) :: xy_SIceSfcTemp(IA,JA)  !< degC
    real(DP), intent(in) :: xy_SeaSfcTemp(IA,JA)   !< K

    ! 局所変数
    ! Local variables
    !    
    integer :: i
    integer :: j
    real(DP) :: emisiv
    real(DP) :: SIceSfcTempK

    ! 実行文; Executable statements
    !
    
    !$omp parallel do private(emisiv, SIceSfcTempK) collapse(2)
    do j = JS, JE
       do i = IS, IE
          
          if ( xy_SIceCon(i,j) > IceMaskMin ) then

             !--
             if ( xy_SnowThick(i,j) <= 0d0 ) then
                xy_SfcAlbedoAI(i,j) = AlbedoIce
                emisiv = EmissivIce
                xy_PenSWRFlxAI(i,j) = I0*(1d0 - AlbedoIce)*xy_SDWRFlx(i,j) ! Note: the sign is positive. 

             elseif( xy_SIceSfcTemp(i,j) >= FreezeTempWater ) then
                xy_SfcAlbedoAI(i,j) = AlbedoMeltSnow
                emisiv = EmissivSnow
                xy_PenSWRFlxAI(i,j) = 0d0

             else
                xy_SfcAlbedoAI(i,j) = AlbedoSnow
                emisiv = EmissivSnow
                xy_PenSWRFlxAI(i,j) = 0d0
             end if

             SIceSfcTempK = degC2K( xy_SIceSfcTemp(i,j) )
             
             xy_SfcHFlxAI(i,j) = &
                  &   xy_LatHFlx(i,j) + xy_SenHFlx(i,j)            &
                  & - emisiv*xy_LDWRFlx(i,j)                       &
                  & - (1d0 - xy_SfcAlbedoAI(i,j))*xy_SDWRFlx(i,j)  &
                  & + xy_PenSWRFlxAI(i,j)                          &
                  & + emisiv*(SBConst*SIceSfcTempK**4)

             xy_DSfcHFlxDTsAI(i,j) = &
                  &   xy_DLatSensHFlxDTs(i,j)                      &
                  & + 4d0*emisiv*(SBConst*SIceSfcTempK**3)
             
          else ! grid cell not covered by sea ice
             xy_SfcAlbedoAI(i,j) = UNDEFVAL
             xy_SfcHFlxAI(i,j)   = UNDEFVAL
          end if

          xy_SfcHFlxAO(i,j) = &
               &   xy_LatHFlx(i,j) + xy_SenHFlx(i,j)            &
               & - EmissivOcean*xy_LDWRFlx(i,j)                 &
               & - (1d0 - AlbedoOcean)*xy_SDWRFlx(i,j)          &
               & + EmissivOcean*(SBConst*xy_SeaSfcTemp(i,j)**4)
          
       end do
    end do
             
  end subroutine DSIce_Boundary_common_CalcSfcHFlx

  !-----------------------------------------

  subroutine DSIce_Boundary_common_CalcBtmHFlx( &
       & xy_BtmHFlxIO,                                  & ! (out)
       & xy_SIceCon, xy_UO, xy_VO, xy_SeaSfcTemp,       & ! (in)
       & xy_OcnFrzTemp, xy_OcnMixLyrDepth,              & ! (in)
       & dt                                             & ! (in)
       & )

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
             xy_BtmHFlxIO(i,j) = UNDEFVAL
          end if
          
       end do
    end do
    
  end subroutine DSIce_Boundary_common_CalcBtmHFlx
  
end module DSIce_Boundary_common_mod
