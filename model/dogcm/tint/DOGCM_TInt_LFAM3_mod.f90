!-------------------------------------------------------------
! Copyright (c) 2017-2017 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_TInt_LFAM3_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_GovernEq_mod, only: &
       & SolverType, DynEqType, EOSType,  &
       & OCNGOVERNEQ_SOLVER_HSPM_VSPM,    &
       & OCNGOVERNEQ_SOLVER_HSPM_VFVM,    &
       & OCNGOVERNEQ_DYN_HYDROBOUSSINESQ, &
       & OCNGOVERNEQ_DYN_NONDYN_MIXEDLYR

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,               &
       & JA, JS, JE, JM, JBLOCK,       &
       & KA, KS, KE, KM,               &
       & xyz_Z, xy_Topo,               &
       & xyz_Lat,                      &
       & z_CDK,                        &
       & DOGCM_Admin_Grid_UpdateVCoord
#include "../admin/DOGCM_Admin_GaussSpmGridIndexDef.h"

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
       & DOGCM_Dyn_driver_OMGDiag2,               &
       & DOGCM_Dyn_driver_HydPresDiag,            &
       & DOGCM_Dyn_driver_UVBarotDiag

  use DOGCM_Phys_driver_mod, only:    &
       & DOGCM_Phys_driver_Do,        &
       & DOGCM_Phys_driver_ImplUV,   &
       & DOGCM_Phys_driver_ImplTRC
       
  use DOGCM_IO_History_mod, only: &
       & DOGCM_IO_History_RegistVar, &
       & DOGCM_IO_History_HistPut
  
  use SpmlUtil_mod    
  use VFvmUtil_mod

  use ProfUtil_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DOGCM_TInt_LFAM3_Init, DOGCM_TInt_LFAM3_Final

  public :: DOGCM_TInt_LFAM3_advance_Dyn
  public :: DOGCM_TInt_LFAM3_advance_Phys
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_TInt_LFAM3_mod' !< Module Name


contains

  !>
  !!
  !!
  Subroutine DOGCM_TInt_LFAM3_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !


    !    call read_nmlData(configNmlName)

!!$    call DOGCM_IO_History_RegistVar( 'PTempRHS_dyn', 'T', 'global mean of PTemp', 'K/s' )
!!$    call DOGCM_IO_History_RegistVar( 'PTempRHS_phy', 'T', 'global mean of PTemp', 'K/s' )
!!$    call DOGCM_IO_History_RegistVar( 'PTempRHS_impl', 'T', 'global mean of PTemp', 'K/s' )
!!$    call DOGCM_IO_History_RegistVar( 'SaltRHS_dyn', 'T', 'global mean of Salt', 'psu/s' )
!!$    call DOGCM_IO_History_RegistVar( 'SaltRHS_phy', 'T', 'global mean of Salt', 'psu/s' )
!!$    call DOGCM_IO_History_RegistVar( 'SaltRHS_impl', 'T', 'global mean of Salt', 'psu/s' )    
!!$    call DOGCM_IO_History_RegistVar( 'Div', 'IJKT', 'horizontal divergence', 's-1' )    
!!$    call DOGCM_IO_History_RegistVar( 'PTemp_t_dyn', 'IJKT', 'PTemp tendency of dynamics', 'K/s' )

    
  end subroutine DOGCM_TInt_LFAM3_Init

  !>
  !!
  !!
  subroutine DOGCM_TInt_LFAM3_Final()

    ! 宣言文; Declaration statement
    !
    
    ! 実行文; Executable statements
    !


  end subroutine DOGCM_TInt_LFAM3_Final

  !-----------------------------------------------------------------------------------


  subroutine DOGCM_TInt_LFAM3_advance_Phys( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,         & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,          & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                  & ! (in)
       & xyz_Z, xy_Topo,                                         & ! (in)
       & dt,                                                     & ! (in)
       & lhst_tend                                               & !(in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP), intent(out) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP), intent(out) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(out) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(out) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(out) :: xy_BtmFrictCoef(IA,JA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_TOPO(IA,JA)
    real(DP), intent(in) :: dt
    logical, intent(in), optional :: lhst_tend

    ! 作業変数
    ! Work variables
    !
    integer :: k
    integer :: n
    real(DP) :: PTempRHS_phy_tend
    real(DP) :: SaltRHS_phy_tend
    
    ! 実行文; Executable statements
    !
    
    != Physical process ======================================================

    call ProfUtil_RapStart('OcnPhys', 2)

    !$omp parallel
    !$omp do
    do k=KS, KE
       xyz_U_RHS_phy(:,:,k) = 0d0
       xyz_V_RHS_phy(:,:,k) = 0d0
    end do
    !$omp do collapse(2)
    do n=1, TRC_TOT_NUM
    do k=KS, KE
       xyza_TRC_RHS_phy(:,:,k,n) = 0d0
    end do
    end do
    !$omp end parallel

    call DOGCM_Phys_driver_Do( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,       & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,        & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                & ! (in)
       & xyz_Z, xy_Topo,                                       & ! (in)
       & dt, lhst_tend                                         & ! (in)
       & )
    
!!$    if (present(lhst_tend)) then
!!$       if (lhst_tend) then
!!$
!!$          select case(SolverType)
!!$          case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)
!!$             PTempRHS_phy_tend = AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyza_TRC_RHS_phy(IS:IE,JS:JE,KS:KE,TRCID_PTEMP)))
!!$             SaltRHS_phy_tend = AvrLonLat_xy(xy_IntSig_BtmToTop_xyz(xyza_TRC_RHS_phy(IS:IE,JS:JE,KS:KE,TRCID_Salt)))
!!$          case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
!!$             PTempRHS_phy_tend = AvrLonLat_xy(VFvm_Int_BtmToTop(xyza_TRC_RHS_phy(IS:IE,JS:JE,:,TRCID_PTEMP), xyz_H))
!!$             SaltRHS_phy_tend = AvrLonLat_xy(VFvm_Int_BtmToTop(xyza_TRC_RHS_phy(IS:IE,JS:JE,:,TRCID_Salt), xyz_H))
!!$          end select
!!$          call DOGCM_IO_History_HistPut( 'PTempRHS_phy', PTempRHS_phy_tend)
!!$          call DOGCM_IO_History_HistPut( 'SaltRHS_phy', SaltRHS_phy_tend )
!!$          
!!$       end if
!!$    end if

    call ProfUtil_RapEnd('OcnPhys', 2)
    
  end subroutine DOGCM_TInt_LFAM3_advance_Phys

  !----------------------------------------------------------------------
  
  subroutine DOGCM_TInt_LFAM3_advance_Dyn(  & 
       & xyz_UA, xyz_VA, xyz_OMGA, xyz_HA,                    & ! (out)
       & xy_SSHA, xyza_TRCA, xy_SfcPresA, xyz_HydPresA,       & ! (out)
       & xyz_UN, xyz_VN, xyz_OMGN, xyz_HN,                    & ! (in)
       & xy_SSHN, xyza_TRCN, xy_SfcPresN, xyz_HydPresN,       & ! (in)
       & xyz_UB, xyz_VB, xyz_OMGB, xyz_HB,                    & ! (in)
       & xy_SSHB, xyza_TRCB, xy_SfcPresB, xyz_HydPresB,       & ! (in)
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,      & ! (in)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,       & ! (in)
       & dt,                                                  & ! (in)
       & alpha, gamma, lambda,                                & ! (in)
       & lhst_tend                                            & ! (in)
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

    real(DP), intent(inout) :: xyz_UN(IA,JA,KA)
    real(DP), intent(inout) :: xyz_VN(IA,JA,KA)
    real(DP), intent(inout) :: xyz_OMGN(IA,JA,KA)
    real(DP), intent(inout) :: xyz_HN(IA,JA,KA)
    real(DP), intent(inout) :: xy_SSHN(IA,JA)
    real(DP), intent(inout) :: xyza_TRCN(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(inout) :: xy_SfcPresN(IA,JA)
    real(DP), intent(inout) :: xyz_HydPresN(IA,JA,KA)

    real(DP), intent(in) :: xyz_UB(IA,JA,KA)
    real(DP), intent(in) :: xyz_VB(IA,JA,KA)
    real(DP), intent(in) :: xyz_OMGB(IA,JA,KA)
    real(DP), intent(in) :: xyz_HB(IA,JA,KA)
    real(DP), intent(in) :: xy_SSHB(IA,JA)
    real(DP), intent(in) :: xyza_TRCB(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xy_SfcPresB(IA,JA)
    real(DP), intent(in) :: xyz_HydPresB(IA,JA,KA)

    real(DP), intent(in) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP), intent(in) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP), intent(in) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
    
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(in) :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: xy_BtmFrictCoef(IA,JA)
    
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha
    real(DP), intent(in) :: gamma
    real(DP), intent(in) :: lambda
    
    logical, intent(in)  :: lhst_tend
    
    ! 作業変数
    ! Work variables
    !
    
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
    real(DP) :: xy_UBarotB(IA,JA)
    real(DP) :: xy_VBarotB(IA,JA)
    real(DP) :: xy_UBarotN(IA,JA)
    real(DP) :: xy_VBarotN(IA,JA)
    real(DP) :: xy_UBarotA(IA,JA)
    real(DP) :: xy_VBarotA(IA,JA)
    real(DP) :: xy_UBarocForce(IA,JA)
    real(DP) :: xy_VBarocForce(IA,JA)
    real(DP) :: xy_SfcPresMomEq(IA,JA)
    real(DP) :: xy_CoriUBarot(IA,JA)
    real(DP) :: xy_CoriVBarot(IA,JA)
    real(DP) :: xy_UVolFlxBarot(IA,JA)
    real(DP) :: xy_VVolFlxBarot(IA,JA)    
    
    real(DP) :: xy_Cori(IA,JA)    

    integer :: i
    integer :: j
    integer :: jj
    integer :: k
    integer :: n

    real(DP) :: dtTRC
    real(DP) :: dtMom
    real(DP) :: dtSSH
    real(DP) :: PresTAvgCoefA

    real(DP) :: U_RHSTmp
    real(DP) :: CoriImplFac
    real(DP) :: BtmFrictImplFac
    
    real(DP) :: PTempRHS_dyn_tend
    real(DP) :: SaltRHS_dyn_tend
    real(DP) :: PTempRHS_impl_tend
    real(DP) :: SaltRHS_impl_tend

    real(DP) :: xyz_Tmp(IA,JA,KA)
    real(DP) :: xyz_Zero(IA,JA,KA)

    real(DP) :: xyz_U(IA,JA,KA)
    real(DP) :: xyz_V(IA,JA,KA)
    real(DP) :: xyz_OMG(IA,JA,KA)
    real(DP) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP) :: xyz_H(IA,JA,KA)
    real(DP) :: xyz_HydPres(IA,JA,KA)
    real(DP) :: xy_SSH(IA,JA)
    real(DP) :: xy_SfcPres(IA,JA)

    real(DP) :: alpha_    
    real(DP) :: gamma_
    real(DP) :: CoriCoefN, CoriCoefB
    
    ! 実行文; Executable statements
    !

    call ProfUtil_RapStart('OcnDyn', 2)
    
    != Preparation ===========================================================
    
    xy_Cori(:,:) = 2d0*Omega*sin(xyz_Lat(:,:,KS))
    xyz_Zero(:,:,:) = 0d0
    
    dtTRC = dt
    dtMom = dt 
    dtSSH = dt 

    != Predictor ============================================================

    
    !-
    
    call DOGCM_Dyn_driver_VorDivDiag( xyz_Vor, xyz_Div, &  ! (out)
         & xyz_UN, xyz_VN                               &  ! (in)
         & )
    
    !- vertical velocity
    !    & pseudo-compressible algorithm

    xyz_HA = xyz_HB
    call DOGCM_Dyn_driver_OMGDiag( xyz_OMGN,      & ! (out)
         & xyz_Div, xyz_HB, xyz_HA, dtSSH         & ! (in)
         & )
    xyz_H  = xyz_HA
    
    !- Updates Baroclinic momentum

    xyz_GeoPot(:,:,:) = Grav*xyz_Z(:,:,:)
    
    !$omp parallel do
    do k = KS, KE
       call EOSDriver_Eval( rhoEdd = xyz_DensEdd(:,:,k),   & ! (out)
            & theta =xyza_TRCN(:,:,k,TRCID_PTEMP),         & ! (in)
            & S     =xyza_TRCN(:,:,k,TRCID_SALT),          & ! (in)
            & p     = -RefDens*xyz_GeoPot(:,:,k)           & ! (in)
            & )
    end do
    
    call DOGCM_Dyn_driver_HydPresDiag( xyz_HydPresN,                 & ! (out)
         & xyz_DensEdd, xyz_HN                                       & ! (in)
         & )
    
    !$omp parallel private(i,j,k)
    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
       ! Note: The surface pressure and colioris terms are not added
       ! if the splitting of baroclinic and barotropic modes is used.
       xyz_CoriU(:,j,k) = xy_Cori(:,j)*(gamma*xyz_UN(:,j,k) + (1d0 - gamma)*xyz_UB(:,j,k))
       xyz_CoriV(:,j,k) = xy_Cori(:,j)*(gamma*xyz_VN(:,j,k) + (1d0 - gamma)*xyz_VB(:,j,k))
       xyz_Pres(:,j,k)  = xyz_HydPresN(:,j,k) + xy_SfcPresN(:,j)
    end do
    end do

    !$omp end parallel
    
    call DOGCM_Dyn_driver_MOMBarocRHS( &
         & xyz_U_RHS, xyz_V_RHS,                         & ! (out)
         & xyz_UN, xyz_VN, xyz_OMGN, xyz_Vor, xyz_Div,   & ! (in)
         & xyz_HN, xyz_Pres, xyz_DensEdd, xyz_GeoPot,    & ! (in)
         & xyz_CoriU, xyz_CoriV,                         & ! (in)
         & xyz_Zero, xyz_Zero                            & ! (in)
         & )
        
    if ( alpha > 0d0 ) then
       !$omp parallel do private(i,j,k,CoriImplFac,U_RHSTmp) collapse(2)
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          CoriImplFac      = 0.5d0*dtMom*xy_Cori(i,j)
          U_RHSTmp         = xyz_U_RHS(i,j,k)

          xyz_U_RHS(i,j,k) = (U_RHSTmp + CoriImplFac*xyz_V_RHS(i,j,k))/(1d0 + CoriImplFac**2)
          xyz_V_RHS(i,j,k) = (xyz_V_RHS(i,j,k) - CoriImplFac*U_RHSTmp)/(1d0 + CoriImplFac**2)

          xyz_U(i,j,k) = 5d0/6d0*dtMom*xyz_U_RHS(i,j,k) + (2d0*xyz_UN(i,j,k) + xyz_UB(i,j,k))/3d0
          xyz_V(i,j,k) = 5d0/6d0*dtMom*xyz_V_RHS(i,j,k) + (2d0*xyz_VN(i,j,k) + xyz_VB(i,j,k))/3d0
       end do
       end do
       end do
    else
       !$omp parallel do
       do k = KS, KE
          xyz_U(:,:,k) = 5d0/6d0*dtMom*xyz_U_RHS(:,:,k) + (2d0*xyz_UN(:,:,k) + xyz_UB(:,:,k))/3d0
          xyz_V(:,:,k) = 5d0/6d0*dtMom*xyz_V_RHS(:,:,k) + (2d0*xyz_VN(:,:,k) + xyz_VB(:,:,k))/3d0
       end do
    end if

    !** Update tracers -------------------------------------------------------
    
    call DOGCM_Dyn_driver_HTRCRHS( xyza_HTRC_RHS,                                                & ! (out)
         & xyza_TRCN, xyz_UN, xyz_VN, xyz_Div, xyz_OMGN, xyz_HN, spread(xyz_Zero,4,TRC_TOT_NUM) )  ! (in)
    
    do n = 1, TRC_TOT_NUM
       !$omp parallel
       !$omp workshare
       xyza_TRC(:,:,:,n) =  &
            & ( (xyz_HN + 0d0)*(2d0*xyza_TRCN(:,:,:,n) + xyza_TRCB(:,:,:,n))/3d0 &
            & + 5d0/6d0*dtTRC*xyza_HTRC_RHS(:,:,:,n) )/(xyz_HN + 0d0)
       !$omp end workshare
       !$omp end parallel
    end do
    
    != Corrector ===============================================================================

    !-

    call DOGCM_Dyn_driver_VorDivDiag( xyz_Vor, xyz_Div, &  ! (out)
         & xyz_U, xyz_V                                 &  ! (in)
         & )
    
    !- vertical velocity
    !    & pseudo-compressible algorithm

    xyz_HA = xyz_HB
    call DOGCM_Dyn_driver_OMGDiag( xyz_OMG,       & ! (out)
         & xyz_Div, xyz_HB, xyz_HA, dtSSH         & ! (in)
         & )
    xyz_H = xyz_HA
    
    !--    
    
    xyz_GeoPot(:,:,:) = Grav*xyz_Z(:,:,:)
    
    !$omp parallel do
    do k = KS, KE
       call EOSDriver_Eval( rhoEdd = xyz_DensEdd(:,:,k),   & ! (out)
            & theta = xyza_TRC(:,:,k,TRCID_PTEMP),         & ! (in)
            & S     = xyza_TRC(:,:,k,TRCID_SALT),          & ! (in)
            & p     = -RefDens*xyz_GeoPot(:,:,k)           & ! (in)
            & )
    end do
    
    call DOGCM_Dyn_driver_HydPresDiag( xyz_HydPres,         & ! (out)
         & xyz_DensEdd, xyz_H                               & ! (in)
         & )


    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarotN, xy_VBarotN,         & ! (out)
         & xyz_UN, xyz_VN, xyz_HN, xy_SSHN, xy_TOPO                    & ! (in)
         & )

    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarotB, xy_VBarotB,       & ! (out)
         & xyz_UB, xyz_VB, xyz_HB, xy_SSHB, xy_TOPO                  & ! (in)
         & )
    
    alpha_ = 0.5d0
    gamma_ = 0.5d0
    CoriCoefN = gamma_ + alpha_
    CoriCoefB = 1d0 - alpha_ - gamma_
    
    !$omp parallel private(i,j,k)
    !$omp do
    do j = JS, JE
       xy_CoriUBarot(:,j)   =   xy_Cori(:,j)*(CoriCoefN*xy_UBarotN(:,j) + CoriCoefB*xy_UBarotB(:,j))
       xy_CoriVBarot(:,j)   =   xy_Cori(:,j)*(CoriCoefN*xy_VBarotN(:,j) + CoriCoefB*xy_VBarotB(:,j))
       
       xy_SfcPresMomEq(:,j) =    CoriCoefN*xy_SfcPresN(:,j) + CoriCoefB*xy_SfcPresB(:,j) &
            &                  - xy_SfcPresN(:,j)
    end do
    
    !$omp do collapse(2)
    do k=KS, KE
    do j=JS, JE
       ! Note: The surface pressure and colioris terms are not added
       ! if the splitting of baroclinic and barotropic modes is used.
       xyz_CoriU(:,j,k) = xy_Cori(:,j)*(CoriCoefN*xyz_UN(:,j,k) + CoriCoefB*xyz_UB(:,j,k))
       xyz_CoriV(:,j,k) = xy_Cori(:,j)*(CoriCoefN*xyz_VN(:,j,k) + CoriCoefB*xyz_VB(:,j,k))
       xyz_Pres(:,j,k)  = xyz_HydPres(:,j,k) + xy_SfcPresN(:,j)
    end do
    end do
    !$omp end parallel
    
    call DOGCM_Dyn_driver_MOMBarocRHS( &
         & xyz_U_RHS, xyz_V_RHS,                         & ! (out)
         & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,      & ! (in)
         & xyz_H, xyz_Pres, xyz_DensEdd, xyz_GeoPot,     & ! (in)
         & xyz_CoriU, xyz_CoriV,                         & ! (in)
         & xyz_U_RHS_phy, xyz_V_RHS_phy                  & ! (in)
         & )

    !-------

    call ProfUtil_RapStart('-OcnImpl_Mom', 3)
    call DOGCM_Phys_driver_ImplUV( xyz_UA, xyz_VA,           & ! (out)
         & xyz_UN, xyz_VN, xyz_U_RHS, xyz_V_RHS,             & ! (in)
         & xyz_H, xyz_VViscCoef, xy_BtmFrictCoef,            & ! (in)
         & dtMom, lambda                                     & ! (in)
         & )

    !$omp parallel do
    do k=KS,KE
       xyz_U_RHS(:,:,k) = (xyz_UA(:,:,k) - xyz_UN(:,:,k))/dtMom
       xyz_V_RHS(:,:,k) = (xyz_VA(:,:,k) - xyz_VN(:,:,k))/dtMom
    end do
    call ProfUtil_RapEnd('-OcnImpl_Mom', 3)       

    !--------    
    

    call DOGCM_Dyn_driver_UVBarotDiag( &
         & xy_UBarocForce, xy_VBarocForce,                      & ! (out)
         & xyz_U_RHS, xyz_V_RHS, xyz_H, xya_SSH, xy_Topo        & ! (in)
         & )
    !$omp parallel do
    do j=JS, JE
       xy_UBarocForce(:,j) = xy_UBarocForce(:,j) - xy_CoriVBarot(:,j)
       xy_VBarocForce(:,j) = xy_VBarocForce(:,j) + xy_CoriUBarot(:,j)
    end do

    if ( alpha_ > 0d0 ) then
       !$omp parallel do private(i,j,k,CoriImplFac,U_RHSTmp) collapse(2)
       do k = KS, KE
       do j = JS, JE
       do i = IS, IE
          CoriImplFac      = alpha_*dtMom*xy_Cori(i,j)
          U_RHSTmp         = xyz_U_RHS(i,j,k)

          xyz_U_RHS(i,j,k) = (U_RHSTmp + CoriImplFac*xyz_V_RHS(i,j,k))/(1d0 + CoriImplFac**2)
          xyz_V_RHS(i,j,k) = (xyz_V_RHS(i,j,k) - CoriImplFac*U_RHSTmp)/(1d0 + CoriImplFac**2)

          xyz_UA(i,j,k) = xyz_UN(i,j,k) + dtMom*xyz_U_RHS(i,j,k)
          xyz_VA(i,j,k) = xyz_VN(i,j,k) + dtMom*xyz_V_RHS(i,j,k)          
       end do
       end do
       end do
    else
       !$omp parallel do
       do k = KS, KE
          xyz_UA(:,:,k) = xyz_UN(:,:,k) + dtMom*xyz_U_RHS(:,:,k)
          xyz_VA(:,:,k) = xyz_VN(:,:,k) + dtMom*xyz_V_RHS(:,:,k)
       end do
    end if

    !- Barotropic mode ----------------------------------------------------------------------------

    call DOGCM_Dyn_driver_MOMBarotRHS( &
       & xy_UBarot_RHS, xy_VBarot_RHS,                      & ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPresMomEq,     & ! (in)
       & xy_UBarocForce, xy_VBarocForce                     & ! (in)
       & )

    if (alpha_ > 0d0) then
       !$omp parallel do private(CoriImplFac) 
       do j = JS, JE
       do i = IS, IE
          CoriImplFac     = alpha_*dtSSH*xy_Cori(i,j)
          xy_UBarotA(i,j) = xy_UBarotN(i,j) + &
               & dtSSH*(xy_UBarot_RHS(i,j) + CoriImplFac*xy_VBarot_RHS(i,j)) / (1d0 + CoriImplFac**2)       

          xy_VBarotA(i,j) = xy_VBarotN(i,j) + &
               & dtSSH*(xy_VBarot_RHS(i,j) - CoriImplFac*xy_UBarot_RHS(i,j)) / (1d0 + CoriImplFac**2)

          xy_SfcPresA(i,j) = xy_SfcPresN(i,j)
       end do
       end do
       
    else
       !$omp parallel
       !$omp workshare
       xy_UBarotA(:,:) = xy_UBarotN + dtSSH*xy_UBarot_RHS
       xy_VBarotA(:,:) = xy_VBarotN + dtSSH*xy_VBarot_RHS
       xy_SfcPresA(:,:) = xy_SfcPresN(:,:)
       !$omp end workshare
       !$omp end parallel
    end if

    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarotA, xy_VBarotA,           & ! (out)
         & xyz_UA, xyz_VA, xyz_H, xy_SSH, xy_TOPO                        & ! (in)
         & )
    
    PresTAvgCoefA = alpha_
    call DOGCM_Dyn_driver_BarotUpdate( &
         & xy_UBarotA, xy_VBarotA,                      & ! (inout)
         & xy_SfcPresA, xy_SSHA,                        & ! (inout)
         & xy_Cori, dtMom, dtSSH, PresTAvgCoefA         & ! (in)
         & )
    
    !----------------------------------------------------------------------
    

    !** Update vertical coordinate -------------------------------------------

    call DOGCM_Admin_Grid_UpdateVCoord( xyz_HA,   & ! (out)
         & xy_SSHA                                & ! (in)
         & )

    
    !- Update full velocity    

    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarot, xy_VBarot,           & ! (out)
         & xyz_UA, xyz_VA, xyz_HA, xy_SSHA, xy_TOPO                    & ! (in)
         & )

    !$omp parallel do collapse(2) private(j,k)
    do k = KS, KE
    do j = JS, JE
       xyz_UA(:,j,k) =   xyz_UA(:,j,k) - xy_UBarot(:,j) + xy_UBarotA(:,j)  
       xyz_VA(:,j,k) =   xyz_VA(:,j,k) - xy_VBarot(:,j) + xy_VBarotA(:,j)

       xyz_U(:,j,k) = (5d0*xyz_UA(:,j,k) + 8d0*xyz_UN(:,j,k) - xyz_UB(:,j,k))/12d0
       xyz_V(:,j,k) = (5d0*xyz_VA(:,j,k) + 8d0*xyz_VN(:,j,k) - xyz_VB(:,j,k))/12d0       
       xyz_H(:,j,k) = (5d0*xyz_HA(:,j,k) + 8d0*xyz_HN(:,j,k) - xyz_HB(:,j,k))/12d0       
    end do
    end do
    do j = JS, JE
       xy_SSH(:,j) = (5d0*xy_SSHA(:,j) + 8d0*xy_SSHN(:,j) - xy_SSHB(:,j))/12d0
    end do

    != Update tracers ****************************************************
    !

    call DOGCM_Dyn_driver_UVBarotDiag( xy_UBarot, xy_VBarot,           & ! (out)
         & xyz_U, xyz_V, xyz_H, xy_SSH, xy_TOPO                        & ! (in)
         & )
    
    !** Modify horizontal momentum such that the vertical integrated mass fluxes is
    !   consistent with surface height equation. 

!!$    !$omp parallel do collapse(2) private(j,k)
!!$    do k = KS, KE
!!$    do j = JS, JE
!!$       xyz_U(:,j,k) = xyz_U(:,j,k) + (xy_UVolFlxBarot(:,j)/(xy_SSH(:,j) + xy_Topo(:,j))) - xy_UBarot(:,j)
!!$       xyz_V(:,j,k) = xyz_V(:,j,k) + (xy_VVolFlxBarot(:,j)/(xy_SSH(:,j) + xy_Topo(:,j))) - xy_VBarot(:,j)
!!$    end do
!!$    end do
 
    !- Diagnose horizontal divergence and relative vorcity ------------------

    call DOGCM_Dyn_driver_VorDivDiag( xyz_Vor, xyz_Div, &  ! (out)
         & xyz_U, xyz_V                                 &  ! (in)
         & )
    
    !- Diagnose vertical velocity -------------------------------------------

    call DOGCM_Dyn_driver_OMGDiag( xyz_OMG,       & ! (out)
         & xyz_Div, xyz_HN, xyz_HA, dtSSH         & ! (in)
         & )
    
    
    call DOGCM_Dyn_driver_HTRCRHS( xyza_HTRC_RHS,                              & ! (out)
         & xyza_TRC, xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H, xyza_TRC_RHS_phy   & ! (in)
         & )

    if (lambda > 0d0) then
       call DOGCM_Phys_driver_ImplTRC( xyza_TRCA,                    & ! (out)
            & xyza_TRCN, xyza_HTRC_RHS,                              & ! (in)
            & xyz_HA, xyz_HN, xyz_VDiffCoef, dtTRC, lambda,          & ! (in)
            & lhst_tend )                                              ! (in)
    else
       do n = 1, TRC_TOT_NUM
          !$omp parallel
          !$omp workshare
          xyza_TRCA(:,:,:,n) = (xyz_HN*xyza_TRCN(:,:,:,n) + dtTRC*xyza_HTRC_RHS(:,:,:,n))/xyz_HA
          !$omp end workshare
          !$omp end parallel
       end do
    end if

    call ProfUtil_RapEnd('OcnDyn', 2)

!!$    write(*,*) "(C) U:", xyz_UA(IS,JS:JE,KS)
!!$    write(*,*) "(C) PTemp:", xyza_TRCA(IS,JS:JE,KS,TRCID_PTEMP)
!!$    if (isNan(xyz_UA(IS,JS,KS)))stop
    
  end subroutine DOGCM_TInt_LFAM3_advance_Dyn
  
end module DOGCM_TInt_LFAM3_mod
