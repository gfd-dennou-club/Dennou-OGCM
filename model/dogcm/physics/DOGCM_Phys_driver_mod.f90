!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Phys_driver_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM,       &       
       & xyz_Z, xy_Topo

  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM,             &
       & TRCID_PTEMP, TRCID_SALT
  
  use DOGCM_VPhys_driver_mod, only: &
       & DOGCM_VPhys_driver_Init, DOGCM_VPhys_driver_Final, &
       & DOGCM_VPhys_driver_UpdateVViscDiffCoef

  use DOGCM_Phys_spm_mod, only: &
       & DOGCM_Phys_spm_Init, DOGCM_Phys_spm_Final, &
       & DOGCM_Phys_spm_Do,           &
       & DOGCM_Phys_spm_VImplUV,      &
       & DOGCM_Phys_spm_VImplTRC
       
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_Phys_driver_mod' !< Module Name


  public :: DOGCM_Phys_driver_Init, DOGCM_Phys_driver_Final
  public :: DOGCM_Phys_driver_Do

  public :: DOGCM_Phys_driver_VImplUV
  public :: DOGCM_Phys_driver_VImplTRC

contains

  !>
  !!
  !!
  Subroutine DOGCM_Phys_driver_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)

    ! Initialize some modules for physical process
    !

    call DOGCM_Phys_spm_Init( configNmlName )
    
    call DOGCM_VPhys_driver_Init( configNmlName )

    
  end subroutine DOGCM_Phys_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_Phys_driver_Final()

    ! 実行文; Executable statements
    !

    call DOGCM_VPhys_driver_Final()

    call DOGCM_Phys_spm_Final()
    
  end subroutine DOGCM_Phys_driver_Final

  

  !-------------------------------------

  subroutine DOGCM_Phys_driver_Do( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,       & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef,                         & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                & ! (in)
       & xyz_Z, xy_Topo,                                       & ! (in)
       & dt                                                    & ! (in)
       & )

    ! 宣言文; Declareration statements
    !    
    real(DP), intent(inout) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP), intent(inout) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
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
    
    call DOGCM_VPhys_driver_UpdateVViscDiffCoef(   &
         & xyz_VViscCoef, xyz_VDiffCoef,               & ! (out)
         & xyz_U, xyz_V, xyza_TRC, xyz_Z               & ! (in)
         & )
    
    call DOGCM_Phys_spm_Do(     &
         & xyz_U_RHS_phy(IS:IE,JS:JE,KS:KE),           & ! (out)
         & xyz_V_RHS_phy(IS:IE,JS:JE,KS:KE),           & ! (out)
         & xyza_TRC_RHS_phy(IS:IE,JS:JE,KS:KE,:),      & ! (out)
         & xyz_VViscCoef(IS:IE,JS:JE,KS:KE),           & ! (in)
         & xyz_VDiffCoef(IS:IE,JS:JE,KS:KE),           & ! (in)
         & xyz_U(IS:IE,JS:JE,KS:KE),                   & ! (in)
         & xyz_V(IS:IE,JS:JE,KS:KE),                   & ! (in)
         & xyz_H(IS:IE,JS:JE,KS:KE),                   & ! (in)
         & xy_SSH(IS:IE,JS:JE),                        & ! (in)
         & xyza_TRC(IS:IE,JS:JE,KS:KE,:),              & ! (in)
         & xyz_Z(IS:IE,JS:JE,KS:KE),                   & ! (in)
         & xy_Topo(IS:IE,JS:JE),                       & ! (in)
         & dt                                          & ! (in)
         & )    
    
  end subroutine DOGCM_Phys_driver_Do

  subroutine DOGCM_Phys_driver_VImplTRC( xyza_TRCA,    & ! (out)
       & xyza_TRC0, xyza_TRC_RHS,                      & ! (in)
       & xyz_HA, xyz_H0, xyz_VDiffCoef, dt, alpha      & ! (in)
       & )

    ! 宣言文; Declareration statements
    !    
    real(DP), intent(out) :: xyza_TRCA(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC0(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_HA(IA,JA,KA)    
    real(DP), intent(in) :: xyz_H0(IA,JA,KA)    
    real(DP), intent(in)  :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 実行文; Executable statements
    !
    
    call DOGCM_Phys_spm_VImplTRC( xyza_TRCA(IS:IE,JS:JE,KS:KE,:),                & ! (out)
         & xyza_TRC0(IS:IE,JS:JE,KS:KE,:), xyza_TRC_RHS(IS:IE,JS:JE,KS:KE,:),    & ! (in)
         & xyz_HA(IS:IE,JS:JE,KS:KE), xyz_H0(IS:IE,JS:JE,KS:KE),                 & ! (in)
         & xyz_VDiffCoef(IS:IE,JS:JE,KS:KE),                                     & ! (in)
         & dt,  alpha                                                            & ! (in)
         & )

  end subroutine DOGCM_Phys_driver_VImplTRC

  subroutine DOGCM_Phys_driver_VImplUV( xyz_UA, xyz_VA,    & ! (out)
       & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,             & ! (in)
       & xyz_H, xyz_VViscCoef,                             & ! (in)
       & dt,  alpha                                        & ! (in)
       & )
    
    ! 宣言文; Declaration statement
    !          
    real(DP), intent(out) :: xyz_UA(IA,JA,KA)
    real(DP), intent(out) :: xyz_VA(IA,JA,KA)
    real(DP), intent(in) :: xyz_U0(IA,JA,KA)
    real(DP), intent(in) :: xyz_V0(IA,JA,KA)
    real(DP), intent(in) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)    
    real(DP), intent(in) :: xyz_VViscCoef(IA,JA,KA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 実行文; Executable statements
    !

    call DOGCM_Phys_spm_VImplUV( &
         & xyz_UA(IS:IE,JS:JE,KS:KE), xyz_VA(IS:IE,JS:JE,KS:KE),       & ! (out)
         & xyz_U0(IS:IE,JS:JE,KS:KE), xyz_V0(IS:IE,JS:JE,KS:KE),       & ! (in)
         & xyz_U_RHS(IS:IE,JS:JE,KS:KE), xyz_V_RHS(IS:IE,JS:JE,KS:KE), & ! (in)
         & xyz_H(IS:IE,JS:JE,KS:KE), xyz_VViscCoef(IS:IE,JS:JE,KS:KE), & ! (in)
         & dt, alpha )

  end subroutine DOGCM_Phys_driver_VImplUV

  
end module DOGCM_Phys_driver_mod
