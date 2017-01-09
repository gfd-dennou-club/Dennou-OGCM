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

  use ProfUtil_mod
  
  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM,       &
       & JA, JS, JE, JM,       &
       & KA, KS, KE, KM,       &       
       & xyz_Z, xy_Topo

  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM,             &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Admin_GovernEq_mod, only: &
       & SolverType,                   &
       & OCNGOVERNEQ_SOLVER_HSPM_VSPM, &
       & OCNGOVERNEQ_SOLVER_HSPM_VFVM       
  
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

  public :: DOGCM_Phys_driver_ImplUV
  public :: DOGCM_Phys_driver_ImplTRC

contains

  !>
  !!
  !!
  Subroutine DOGCM_Phys_driver_Init(configNmlName)

    use DOGCM_VPhys_driver_mod, only: &
         & DOGCM_VPhys_driver_Init

    use DOGCM_Phys_spm_mod, only: &
         & DOGCM_Phys_spm_Init

    use DOGCM_Phys_hspm_vfvm_mod, only: &
         & DOGCM_Phys_hspm_vfvm_Init
    
    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)

    ! Initialize some modules for physical process
    !

    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)
       call DOGCM_Phys_spm_Init( configNmlName )
    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
       call DOGCM_Phys_hspm_vfvm_Init( configNmlName )
    end select
    
    call DOGCM_VPhys_driver_Init( configNmlName )

    
  end subroutine DOGCM_Phys_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_Phys_driver_Final()

    use DOGCM_VPhys_driver_mod, only: &
       & DOGCM_VPhys_driver_Final

    use DOGCM_Phys_spm_mod, only: &
         & DOGCM_Phys_spm_Final

    use DOGCM_Phys_hspm_vfvm_mod, only: &
         & DOGCM_Phys_hspm_vfvm_Final

    ! 実行文; Executable statements
    !

    call DOGCM_VPhys_driver_Final()

    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)
       call DOGCM_Phys_spm_Final()
    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
       call DOGCM_Phys_hspm_vfvm_Final()
    end select
    
    
  end subroutine DOGCM_Phys_driver_Final

  

  !--------------------------------------------------

  subroutine DOGCM_Phys_driver_Do( &
       & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,       & ! (out)
       & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,        & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,                & ! (in)
       & xyz_Z, xy_Topo,                                       & ! (in)
       & dt                                                    & ! (in)
       & )

    use DOGCM_VPhys_driver_mod, only: &
         & DOGCM_VPhys_driver_UpdateVViscDiffCoef

    use DOGCM_Phys_spm_mod, only: &
         & DOGCM_Phys_spm_Do

    use DOGCM_Phys_hspm_vfvm_mod, only: &
         & DOGCM_Phys_hspm_vfvm_Do
    
    ! 宣言文; Declareration statements
    !    
    real(DP), intent(inout) :: xyz_U_RHS_phy(IA,JA,KA)
    real(DP), intent(inout) :: xyz_V_RHS_phy(IA,JA,KA)
    real(DP), intent(inout) :: xyza_TRC_RHS_phy(IA,JA,KA,TRC_TOT_NUM)
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

    ! 実行文; Executable statements
    !

    call ProfUtil_RapStart('OcnPhys_Do', 3)
    
    call DOGCM_VPhys_driver_UpdateVViscDiffCoef(   &
         & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,  & ! (out)
         & xyz_U, xyz_V, xyza_TRC,                         & ! (in)
         & xyz_H, xyz_Z, xy_TOPO                           & ! (in)
         & )

    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)
       call DOGCM_Phys_spm_Do(     &
            & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,  & ! (inout)
            & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,   & ! (out)
            & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,           & ! (in)
            & xyz_Z, xy_Topo,                                  & ! (in)
            & dt                                               & ! (in)
            & )    
    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
       call DOGCM_Phys_hspm_vfvm_Do(     &
            & xyz_U_RHS_phy, xyz_V_RHS_phy, xyza_TRC_RHS_phy,  & ! (inout)
            & xyz_VViscCoef, xyz_VDiffCoef, xy_BtmFrictCoef,   & ! (out)
            & xyz_U, xyz_V, xyz_H, xy_SSH, xyza_TRC,           & ! (in)
            & xyz_Z, xy_Topo,                                  & ! (in)
            & dt                                               & ! (in)
            & )    
    end select
    
    call ProfUtil_RapEnd('OcnPhys_Do', 3)
    
  end subroutine DOGCM_Phys_driver_Do

  !------------------------------------------------------------
  
  subroutine DOGCM_Phys_driver_ImplTRC( xyza_TRCA,    & ! (out)
       & xyza_TRC0, xyza_HTRC_RHS,                     & ! (in)
       & xyz_HA, xyz_H0, xyz_VDiffCoef, dt, alpha      & ! (in)
       & )

    use DOGCM_Phys_spm_mod, only: &
         & DOGCM_Phys_spm_ImplTRC

    use DOGCM_Phys_hspm_vfvm_mod, only: &
         & DOGCM_Phys_hspm_vfvm_ImplTRC


    ! 宣言文; Declareration statements
    !    
    real(DP), intent(out) :: xyza_TRCA(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC0(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_HTRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_HA(IA,JA,KA)    
    real(DP), intent(in) :: xyz_H0(IA,JA,KA)    
    real(DP), intent(in)  :: xyz_VDiffCoef(IA,JA,KA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 実行文; Executable statements
    !

    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)
       call DOGCM_Phys_spm_ImplTRC( xyza_TRCA,                                 & ! (out)
            & xyza_TRC0, xyza_HTRC_RHS, xyz_HA, xyz_H0, xyz_VDiffCoef,          & ! (in)
            & dt,  alpha                                                        & ! (in)
            & )
    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
       call DOGCM_Phys_hspm_vfvm_ImplTRC( xyza_TRCA,                           & ! (out)
            & xyza_TRC0, xyza_HTRC_RHS, xyz_HA, xyz_H0, xyz_VDiffCoef,          & ! (in)
            & dt,  alpha                                                        & ! (in)
            & )
    end select
    

  end subroutine DOGCM_Phys_driver_ImplTRC

  !--------------------------------------------------
  
  subroutine DOGCM_Phys_driver_ImplUV( xyz_UA, xyz_VA,    & ! (out)
       & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,             & ! (in)
       & xyz_H, xyz_VViscCoef, xy_BtmFrictCoef,            & ! (in)
       & dt,  alpha                                        & ! (in)
       & )

    use DOGCM_Phys_spm_mod, only: &
         & DOGCM_Phys_spm_ImplUV

    use DOGCM_Phys_hspm_vfvm_mod, only: &
         & DOGCM_Phys_hspm_vfvm_ImplUV
    
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
    real(DP), intent(in) :: xy_BtmFrictCoef(IA,JA)
    real(DP), intent(in) :: dt
    real(DP), intent(in) :: alpha

    ! 実行文; Executable statements
    !


    select case(SolverType)
    case(OCNGOVERNEQ_SOLVER_HSPM_VSPM)
       call DOGCM_Phys_spm_ImplUV( &
            & xyz_UA, xyz_VA,                                             & ! (out)
            & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,                       & ! (in)
            & xyz_H, xyz_VViscCoef,  xy_BtmFrictCoef,                     & ! (in)            
            & dt, alpha )

    case(OCNGOVERNEQ_SOLVER_HSPM_VFVM)
       call DOGCM_Phys_hspm_vfvm_ImplUV( &
            & xyz_UA, xyz_VA,                                             & ! (out)
            & xyz_U0, xyz_V0, xyz_U_RHS, xyz_V_RHS,                       & ! (in)
            & xyz_H, xyz_VViscCoef,  xy_BtmFrictCoef,                     & ! (in)            
            & dt, alpha )
    end select
    
    
  end subroutine DOGCM_Phys_driver_ImplUV

  
end module DOGCM_Phys_driver_mod
