!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_TInt_common_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_GovernEq_mod, only: &
       & DynEqType, EOSType, &
       & OCNGOVERNEQ_DYN_HYDROBOUSSINESQ, &
       & OCNGOVERNEQ_DYN_NONDYN_MIXEDLYR

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM,       &
       & xyz_Z, xy_Topo,       &
       & xyz_Lat,              &
       & DOGCM_Admin_Grid_UpdateVCoord

  use DOGCM_Admin_Constants_mod, only: &
       & Grav, RefDens, Omega
  
  use EOSDriver_mod, only: &
       & EOSDriver_Init, EOSDriver_Final, &
       & EOSDriver_Eval
    
  use DOGCM_Admin_Variable_mod, only: &
       & xya_SSH, xyza_H, xyza_U, xyza_V, xyza_OMG, xyzaa_TRC, &
       & xyza_HydPres, xya_SfcPres,                               &
       & xyz_VViscCoef, xyz_VDiffCoef,                            &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Boundary_vars_mod, only:   &
       & xy_FreshWtFlx, xy_FreshWtFlxS
  
  use DOGCM_Boundary_driver_mod, only: &
       & DOGCM_Boundary_driver_ApplyBC
  
  use DOGCM_Dyn_driver_mod, only: &
       & DOGCM_Dyn_driver_SSHRHS,                 &
       & DOGCM_Dyn_driver_HTRCRHS,                &
       & DOGCM_Dyn_driver_MOMBarocRHS,            &
       & DOGCM_Dyn_driver_MOMBarotRHS,            &
       & DOGCM_Dyn_driver_BarotUpdate,            &
       & DOGCM_Dyn_driver_VorDivDiag,             &
       & DOGCM_Dyn_driver_OMGDiag,                &
       & DOGCM_Dyn_driver_HydPresDiag,            &
       & DOGCM_Dyn_driver_UVBarotDiag

  use DOGCM_Phys_driver_mod, only:    &
       & DOGCM_Phys_driver_Do,        &
       & DOGCM_Phys_driver_VImplUV,   &
       & DOGCM_Phys_driver_VImplTRC
       
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_TInt_common_Init, DOGCM_TInt_common_Final

  public :: DOGCM_TInt_common_advance_Dyn
  public :: DOGCM_TInt_common_advance_Phys
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_TInt_common_mod' !< Module Name


contains

  !>
  !!
  !!
  Subroutine DOGCM_TInt_common_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !


    !    call read_nmlData(configNmlName)

    
  end subroutine DOGCM_TInt_common_Init

  !>
  !!
  !!
  subroutine DOGCM_TInt_common_Final()

    ! 宣言文; Declaration statement
    !
    
    ! 実行文; Executable statements
    !


  end subroutine DOGCM_TInt_common_Final

  !-----------------------------------------------------------------------------------

  subroutine DOGCM_TInt_common_advance_Phys( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,         & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef,                           & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                  & ! (in)
       & xyz_Z, xy_Topo,                                         & ! (in)
       & dt                                                      & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP), intent(out) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP), intent(out) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(out) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(out) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_TOPO(IA,JA)
    real(DP), intent(in) :: dt

    ! 実行文; Executable statements
    !
    
    != Physical process ======================================================

    !$omp parallel
    !$omp workshare
    xyz_U_RHS_phy(:,:,:) = 0d0
    xyz_V_RHS_phy(:,:,:) = 0d0
    xyza_TRC_RHS_phy(:,:,:,:) = 0d0
    !$omp end workshare
    !$omp end parallel

    call DOGCM_Phys_driver_Do( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,       & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef,                         & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                & ! (in)
       & xyz_Z, xy_Topo,                                       & ! (in)
       & dt                                                    & ! (in)
       & )
    
    
  end subroutine DOGCM_TInt_common_advance_Phys

  !----------------------------------------------------------------------
  
  subroutine DOGCM_TInt_common_advance_Dyn(  & 
       & xyz_UA, xyz_VA, xyz_OMGA, xyz_HA,                    & ! (out)
       & xy_SSHA, xyza_TRCA, xy_SfcPresA, xyz_HydPresA,       & ! (out)
       & xyz_U, xyz_V, xyz_OMG, xyz_H,                        & ! (in)
       & xy_SSH, xyza_TRC, xy_SfcPres, xyz_HydPres,           & ! (in)
       & xyz_U0, xyz_V0, xyz_OMG0, xyz_H0,                    & ! (in)
       & xy_SSH0, xyza_TRC0, xy_SfcPres0, xyz_HydPres0,       & ! (in)
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,      & ! (in)
       & xyz_VViscCoef, xyz_VDiffCoef,                        & ! (in)
       & dt,                                                  & ! (in)
       & alpha, gamma, lambda                                 & ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !
    real(DP), intent(inout) :: xyz_UA(IA,JA,KA)
    real(DP), intent(inout) :: xyz_VA(IA,JA,KA)
    real(DP), intent(inout) :: xyz_OMGA(IA,JA,KA)
    real(DP), intent(inout) :: xyz_HA(IA,JA,KA)
    real(DP), intent(inout) :: xy_SSHA(IA,JA)
    real(DP), intent(inout) :: xyza_TRCA(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(inout) :: xy_SfcPresA(IA,JA)
    real(DP), intent(inout) :: xyz_HydPresA(IA,JA,KA)

    real(DP), intent(inout) :: xyz_U(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V(IA,JA,KA)
    real(DP), intent(inout) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(inout) :: xyz_H(IA,JA,KA)
    real(DP), intent(inout) :: xy_SSH(IA,JA)
    real(DP), intent(inout) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(inout) :: xy_SfcPres(IA,JA)
    real(DP), intent(inout) :: xyz_HydPres(IA,JA,KA)

    real(DP), intent(in) :: xyz_U0(IA,JA,KA)
    real(DP), intent(in) :: xyz_V0(IA,JA,KA)
    real(DP), intent(in) :: xyz_OMG0(IA,JA,KA)
    real(DP), intent(in) :: xyz_H0(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH0(IA,JA)
    real(DP), intent(in) :: xyza_TRC0(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xy_SfcPres0(IA,JA)
    real(DP), intent(in) :: xyz_HydPres0(IA,JA,KA)

    real(DP), intent(in) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP), intent(in) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP), intent(in) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(IA,JA,KA)

    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha
    real(DP), intent(in) :: gamma
    real(DP), intent(in) :: lambda
    

    
    ! 作業変数
    ! Work variables
    
    real(DP) :: xy_SSH_RHS(IA,JA)
    real(DP) :: xyz_H_RHS(IA,JA,KA)
    real(DP) :: xyza_HTRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP) :: xyz_U_RHS(IA,JA,KA)
    real(DP) :: xyz_V_RHS(IA,JA,KA)
    real(DP) :: xy_UBarot_RHS(IA,JA)
    real(DP) :: xy_VBarot_RHS(IA,JA)
    
    real(DP) :: xyz_Vor(IA,JA,KA)
    real(DP) :: xyz_Div(IA,JA,KA)
    real(DP) :: xyz_CoriU(IA,JA,KA)
    real(DP) :: xyz_CoriV(IA,JA,KA)
    real(DP) :: xyz_DensEdd(IA,JA,KA)
    real(DP) :: xyz_Pres(IA,JA,KA)
    real(DP) :: xyz_GeoPot(IA,JA,KA)

    real(DP) :: xy_UBarot(IA,JA)
    real(DP) :: xy_VBarot(IA,JA)
    real(DP) :: xy_UBarot0(IA,JA)
    real(DP) :: xy_VBarot0(IA,JA)
    real(DP) :: xy_UBarotA(IA,JA)
    real(DP) :: xy_VBarotA(IA,JA)
    real(DP) :: xy_UBarocForce(IA,JA)
    real(DP) :: xy_VBarocForce(IA,JA)
    real(DP) :: xy_SfcPresMomEq(IA,JA)
    real(DP) :: xy_CoriUBarot(IA,JA)
    real(DP) :: xy_CoriVBarot(IA,JA)

    real(DP) :: xy_Cori(IA,JA)    
    
    integer :: k
    integer :: n

    real(DP) :: dtSSH
    real(DP) :: PresTAvgCoefA

    real(DP) :: xy_U_RHSTmp(IA,JA)
    real(DP) :: xy_CoriImplFac(IA,JA)
    
    ! 実行文; Executable statements
    !

    != Preparation ===========================================================
    
    xyz_GeoPot(:,:,:) = Grav*xyz_Z(:,:,:)
    xy_Cori(:,:) = 2d0*Omega*sin(xyz_Lat(:,:,KS))
    
    != Dynamical process ====================================================

    !** Update sea surface height --------------------------------------------

    call DOGCM_Dyn_driver_SSHRHS( xy_SSH_RHS,              & ! (out)
         & xy_SSH, xy_Topo, xyz_U, xyz_V, xy_FreshWtFlx    & ! (in)
         & )
    xy_SSHA = xy_SSH0 + dt * xy_SSH_RHS

    !** Update vertical coordinate -------------------------------------------

    call DOGCM_Admin_Grid_UpdateVCoord( xyz_HA,   & ! (out)
         & xy_SSHA                       & ! (in)
         & )

    !** Diagnose horizontal divergence and relative vorcity ------------------

    call DOGCM_Dyn_driver_VorDivDiag( xyz_Vor, xyz_Div, &  ! (out)
         & xyz_U, xyz_V                                 &  ! (in)
         & )
    
    !** Diagnose vertical velocity -------------------------------------------

    call DOGCM_Dyn_driver_OMGDiag( xyz_OMG,    & ! (out)
         & xyz_Div, xyz_H0, xyz_HA, dt         & ! (in)
         & )
    
    !** Update tracers -------------------------------------------------------

    call DOGCM_Dyn_driver_HTRCRHS( xyza_HTRC_RHS,                              & ! (out)
         & xyza_TRC, xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H, xyza_TRC_RHS_phy   & ! (in)
         & )

    if (lambda > 0d0) then
       call DOGCM_Phys_driver_VImplTRC( xyza_TRCA,               & ! (out)
            & xyza_TRC0, xyza_HTRC_RHS,                          & ! (in)
            & xyz_HA, xyz_H0, xyz_VDiffCoef, dt, lambda          & ! (in)
            & )

    else
       do n = 1, TRC_TOT_NUM
          !$omp parallel
          !$omp workshare
          xyza_TRCA(:,:,:,n) = (xyz_H0*xyza_TRC0(:,:,:,n) + dt*xyza_HTRC_RHS(:,:,:,n))/xyz_HA
          !$omp end workshare
          !$omp end parallel

       end do
    end if
    
    !** Update momentum -----------------------------------------------------
    
    ! - Calculate density and Diagnose hydrostatic pressure 

    call EOSDriver_Eval( rhoEdd = xyz_DensEdd,                       & ! (out)
         & theta = xyza_TRC(:,:,:,TRCID_PTEMP),                  & ! (in)
         & s     = xyza_TRC(:,:,:,TRCID_SALT),                   & ! (in)
         & p     = -RefDens*xyz_GeoPot                               & ! (in)
         & )                 

    call DOGCM_Dyn_driver_HydPresDiag( xyz_HydPres,                  & ! (out)
         & xyz_DensEdd, xyz_H                                        & ! (in)
         & )

    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarot, xy_VBarot,         & ! (out)
         & xyz_U, xyz_V, xyz_H, xy_SSH, xy_TOPO                      & ! (in)
         & )

    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarot0, xy_VBarot0,       & ! (out)
         & xyz_U0, xyz_V0, xyz_H0, xy_SSH0, xy_TOPO                  & ! (in)
         & )
    
    !- Baroclinic mode

    !$omp parallel
    !$omp workshare
    xy_CoriUBarot(:,:) = xy_Cori*(gamma*xy_UBarot + (1d0 - gamma)*xy_UBarot0)
    xy_CoriVBarot(:,:) = xy_Cori*(gamma*xy_VBarot + (1d0 - gamma)*xy_VBarot0)
!    xy_SfcPresMomEq(:,:) = xy_SfcPres !gamma*xy_SfcPres + (1d0 - gamma)*xy_SfcPres0
    !$omp end workshare

    !$omp do
    do k = KS, KE
       ! Note: The surface pressure and colioris terms are not added
       ! if the splitting of baroclinic and barotropic modes is used.
       xyz_CoriU(:,:,k) =   xy_Cori*(gamma*xyz_U(:,:,k) + (1d0 - gamma)*xyz_U0(:,:,k))  &
            &             - xy_CoriUBarot
       xyz_CoriV(:,:,k) =   xy_Cori*(gamma*xyz_V(:,:,k) + (1d0 - gamma)*xyz_V0(:,:,k))  &
            &             - xy_CoriVBarot
       xyz_Pres(:,:,k) = xyz_HydPres(:,:,k) 
    end do
    !$omp end parallel
    
    call DOGCM_Dyn_driver_MOMBarocRHS( &
         & xyz_U_RHS, xyz_V_RHS,                      & ! (out)
         & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,   & ! (in)
         & xyz_H, xyz_Pres, xyz_DensEdd, xyz_GeoPot,  & ! (in)
         & xyz_CoriU, xyz_CoriV,                      & ! (in)
         & xyz_U_RHS_phy, xyz_V_RHS_phy               & ! (in)
         & )

    if (lambda > 0d0) then
       call DOGCM_Phys_driver_VImplUV( xyz_UA, xyz_VA,          & ! (out)
            & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,             & ! (in)
            & xyz_H, xyz_VViscCoef,                             & ! (in)
            & dt, lambda                                        & ! (in)
            & )
       !$omp parallel
       !$omp workshare
       xyz_U_RHS(:,:,:) = (xyz_UA - xyz_U0)/dt
       xyz_V_RHS(:,:,:) = (xyz_VA - xyz_V0)/dt
       !$omp end workshare
       !$omp end parallel
    end if
    
    if ( alpha > 0d0 ) then
       call DOGCM_Dyn_driver_UVBarotDiag( &
            & xy_UBarocForce, xy_VBarocForce,                      & ! (out)
            & xyz_U_RHS, xyz_V_RHS, xyz_H, xya_SSH, xy_Topo        & ! (in)
            & )
       
       xy_CoriImplFac(:,:) = alpha*dt*xy_Cori
       !$omp parallel do private(xy_U_RHSTmp)
       do k = KS, KE
          xy_U_RHSTmp(:,:) = xyz_U_RHS(:,:,k) - xy_UBarocForce
          xyz_V_RHS(:,:,k) = xyz_V_RHS(:,:,k) - xy_VBarocForce

          xyz_U_RHS(:,:,k) = (xy_U_RHSTmp + xy_CoriImplFac*xyz_V_RHS(:,:,k))/(1d0 + xy_CoriImplFac**2)
          xyz_V_RHS(:,:,k) = (xyz_V_RHS(:,:,k) - xy_CoriImplFac*xy_U_RHSTmp)/(1d0 + xy_CoriImplFac**2)

          xyz_UA(:,:,k) = xyz_U0(:,:,k) + dt*xyz_U_RHS(:,:,k)
          xyz_VA(:,:,k) = xyz_V0(:,:,k) + dt*xyz_V_RHS(:,:,k)          
       end do
    else
       !$omp parallel do
       do k = KS, KE
          xyz_UA(:,:,k) = xyz_U0(:,:,k) + dt*xyz_U_RHS(:,:,k)
          xyz_VA(:,:,k) = xyz_V0(:,:,k) + dt*xyz_V_RHS(:,:,k)
       end do

       call DOGCM_Dyn_driver_UVBarotDiag( &
            & xy_UBarocForce, xy_VBarocForce,                      & ! (out)
            & xyz_U_RHS, xyz_V_RHS, xyz_H, xya_SSH, xy_Topo        & ! (in)
            & )

    end if

    
    !- Barotropic mode
    
    call DOGCM_Dyn_driver_MOMBarotRHS( &
       & xy_UBarot_RHS, xy_VBarot_RHS,                      & ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,     & ! (in)
       & xy_UBarocForce, xy_VBarocForce                     & ! (in)
       & )

    if (alpha > 0d0) then
       !$omp parallel
       !$omp workshare
       xy_CoriImplFac(:,:) = alpha*dt*xy_Cori
       xy_UBarotA(:,:) = xy_UBarot0 + &
            & dt*(xy_UBarot_RHS + xy_CoriImplFac*xy_VBarot_RHS)/(1d0 + xy_CoriImplFac**2)       
       xy_VBarotA(:,:) = xy_VBarot0 + &
            & dt*(xy_VBarot_RHS - xy_CoriImplFac*xy_UBarot_RHS)/(1d0 + xy_CoriImplFac**2)       
       xy_SfcPresA(:,:) = xy_SfcPres !MomEq
       !$omp end workshare
       !$omp end parallel
    else
       !$omp parallel
       !$omp workshare
       xy_UBarotA(:,:) = xy_UBarot0 + dt*xy_UBarot_RHS
       xy_VBarotA(:,:) = xy_VBarot0 + dt*xy_VBarot_RHS
       xy_SfcPresA(:,:) = xy_SfcPres !MomEq
       !$omp end workshare
       !$omp end parallel
    end if
    
    PresTAvgCoefA = 1d0 !alpha
    dtSSH = dt
    call DOGCM_Dyn_driver_BarotUpdate( &
         & xy_UBarotA, xy_VBarotA,                   & ! (inout)
         & xy_SfcPresA, xy_SSHA,                     & ! (inout)
         & xy_Cori, dt, dtSSH, PresTAvgCoefA         & ! (in)
         & )


    !- Update full velocity
    
    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarot, xy_VBarot,    & ! (out)
         & xyz_UA, xyz_VA, xyz_HA, xy_SSHA, xy_TOPO             & ! (in)
         & )
    
    !$omp parallel
    !$omp do
    do k = KS, KE
       xyz_UA(:,:,k) =   xyz_UA(:,:,k) - xy_UBarot + xy_UBarotA  
       xyz_VA(:,:,k) =   xyz_VA(:,:,k) - xy_VBarot + xy_VBarotA 
    end do
    !$omp end parallel
    
  end subroutine DOGCM_TInt_common_advance_Dyn
  
end module DOGCM_TInt_common_mod
