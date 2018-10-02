module OcnDiag_TendencyAnalysis_mod

  ! モジュール引用; Use statement
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use gtool_history

  !* Dennou-OGCM / SeaIce

  use DOGCM_Admin_Constants_mod, only: &
       & PI, RPlanet,               &
       & hViscCoef, hHyperViscCoef, &
       & hDiffCoef, hHyperDiffCoef
         
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM,        &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  use DOGCM_Admin_GovernEq_mod, only: &
       & isPhysicsCompActivated, &
       & OCNGOVERNEQ_LPHYS_MIXMOM_NAME,  &
       & OCNGOVERNEQ_LPHYS_MIXTRC_NAME,  &
       & OCNGOVERNEQ_LPHYS_REDIGM_NAME, &
       & OCNGOVERNEQ_VPHYS_MIXMOM_NAME,  &
       & OCNGOVERNEQ_VPHYS_MIXTRC_NAME,  &
       & OCNGOVERNEQ_VPHYS_CONVEC_NAME

  !-- Module for the parametrizations of lateral oceanic physics

  use LPhys_DIFF_spm_mod, only: &
       & LPhys_DIFF_spm_Init, LPhys_DIFF_spm_Final, &
       & LPhys_DIFF_spm_LMixMOMRHS,                       &
       & LPhys_DIFF_spm_LMixMOMRHSImpl,                   &
       & LPhys_DIFF_spm_LMixTRCRHS,                       &
       & LPhys_DIFF_spm_LMixTRCRHSImpl
       
  use LPhys_RediGM_hspm_vfvm_mod, only: &
       & LPhys_RediGM_hspm_vfvm_Init, LPhys_RediGM_hspm_vfvm_Final, &
       & LPhys_RediGM_hspm_vfvm_AddMixingTerm

  !-- Module for the parametrizations of vertical oceanic physics
  
  use DOGCM_VPhys_ConvAdjust_mod, only: &
       & DOGCM_VPhys_ConvAdjust_Init, DOGCM_VPhys_ConvAdjust_Final, &
       & DOGCM_VPhys_ConvAdjust_AddMixingTerm

  
  
  ! 宣言文; Declareration statements
  !  
  implicit none

  public :: OcnDiag_TendencyAnalysis_Init
  public :: OcnDiag_TendencyAnalysis_Final

  public :: OcnDiag_TendencyAnalysis_perform
  
contains

  subroutine OcnDiag_TendencyAnalysis_Init()

  end subroutine OcnDiag_TendencyAnalysis_Init

  subroutine OcnDiag_TendencyAnalysis_Final()

  end subroutine OcnDiag_TendencyAnalysis_Final

  subroutine OcnDiag_TendencyAnalysis_prepair_Output()

    call DOGCM_IO_History_RegistVar( 'PTemp_t_lmix')
    
  end subroutine OcnDiag_TendencyAnalysis_prepair_Output
  
  subroutine OcnDiag_TendencyAnalysis_perform( &
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                & ! (in)
       & xyz_Z, xy_Topo                                        & ! (in)
       & )

    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)    
    real(DP), intent(in) :: xy_SSH(IA,JA)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_TOPO(IA,JA)

    real(DP) :: xyza_TRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    
!!$    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_MIXTRC_NAME )  ) then
!!$       xyza_TRC_RHS(:,:,:,:) = 0d0       
!!$       call LPhys_DIFF_spm_LMixTRCRHSImpl( xyza_TRC_RHS,           & ! (inout)
!!$            & xyza_TRC, xyz_H, hDiffCoef, hHyperDiffCoef, dt       &  ! (in)
!!$            & )
!!$       
!!$    end if
!!$
!!$    !-- Lateral mixing of tracers by eddy induced velocity -------------------------
!!$    
!!$    if ( isPhysicsCompActivated( OCNGOVERNEQ_LPHYS_REDIGM_NAME )  ) then
!!$       call LPhys_RediGM_hspm_vfvm_AddMixingTerm( &
!!$            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP), xyza_TRC_RHS_phy(:,:,:,TRCID_SALT),  & ! (inout)
!!$            & xyz_VDiffCoef,                                                            & ! (inout)
!!$            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),                  & ! (in)
!!$            & xyz_H, xyz_Z, xy_Topo                                                     & ! (in)
!!$            & )
!!$    end if
!!$
!!$    !-- Vertical mixing of tracers by non-penetrative convection -------------------
!!$    
!!$    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_CONVEC_NAME )  ) then
!!$       call DOGCM_VPhys_ConvAdjust_AddMixingTerm( &
!!$            & xyza_TRC_RHS_phy(:,:,:,TRCID_PTEMP),                         & ! (inout)
!!$            & xyza_TRC_RHS_phy(:,:,:,TRCID_SALT),                          & ! (inout)
!!$            & xyz_ConvIndex,                                               & ! (inout)
!!$            & xyza_TRC(:,:,:,TRCID_PTEMP), xyza_TRC(:,:,:,TRCID_SALT),     & ! (in)
!!$            & xyz_Z, z_KAXIS_Weight, dt                                    & ! (in)
!!$            & )
!!$    end if
!!$
!!$
!!$    if ( isPhysicsCompActivated( OCNGOVERNEQ_VPHYS_MIXTRC_NAME )  ) then
!!$       call DOGCM_Phys_hspm_vfvm_VMixTRCRHS( xyza_TRC_RHS_phy,   & ! (inout)
!!$            & xyza_TRC, xyz_H, xyz_VDiffCoef                     & ! (in)
!!$            & )
!!$    end if
!!$    
  end subroutine OcnDiag_TendencyAnalysis_perform
  
end module OcnDiag_TendencyAnalysis_mod
