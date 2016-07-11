!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Dyn_driver_mod

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM

  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Dyn_spm_mod
  
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
  character(*), parameter:: module_name = 'DOGCM_Dyn_driver_mod' !< Module Name


  public :: DOGCM_Dyn_driver_Init, DOGCM_Dyn_driver_Final

  public :: DOGCM_Dyn_driver_SSHRHS
  public :: DOGCM_Dyn_driver_HTRCRHS
  public :: DOGCM_Dyn_driver_MOMBarocRHS
  public :: DOGCM_Dyn_driver_MOMBarotRHS
  public :: DOGCM_Dyn_driver_BarotUpdate
  
  public :: DOGCM_Dyn_driver_OMGDiag
  public :: DOGCM_Dyn_driver_HydPresDiag
  public :: DOGCM_Dyn_driver_VorDivDiag
  public :: DOGCM_Dyn_driver_UVBarotDiag
  
contains

  !>
  !!
  !!
  Subroutine DOGCM_Dyn_driver_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)

    call DOGCM_Dyn_spm_Init( configNmlName )
    
  end subroutine DOGCM_Dyn_driver_Init

  !>
  !!
  !!
  subroutine DOGCM_Dyn_driver_Final()

    ! 実行文; Executable statements
    !

    call DOGCM_Dyn_spm_Final()
    
  end subroutine DOGCM_Dyn_driver_Final

  

  !-------------------------------------

  subroutine DOGCM_Dyn_driver_SSHRHS( xy_SSH_RHS,                  &  ! (out)
       & xy_SSH, xy_TotDepBasic, xyz_U, xyz_V, xy_FreshWtFlx   )      ! (in)

    real(DP), intent(out) :: xy_SSH_RHS(IA,JA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    real(DP), intent(in) :: xy_TotDepBasic(IA,JA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xy_FreshWtFlx(IA,JA)

!    call MessageNotify('M', module_name, "SSHRHS..")

    call DOGCM_Dyn_spm_SSHRHS( xy_SSH_RHS(IS:IE,JS:JE),        &  ! (out)
         & xy_SSH(IS:IE,JS:JE), xy_TotDepBasic(IS:IE,JS:JE),   &  ! (in)
         & xyz_U(IS:IE,JS:JE,KS:KE), xyz_V(IS:IE,JS:JE,KS:KE), &  ! (in)
         & xy_FreshWtFlx(IS:IE,JS:JE)   )                         ! (in)
    
  end subroutine DOGCM_Dyn_driver_SSHRHS

  !-------------------------------------

  subroutine DOGCM_Dyn_driver_HTRCRHS( xyza_HTRC_RHS,                 & ! (out)
       & xyza_TRC, xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H,             & ! (in)
       & xyza_HTRC_RHS_phys                                          )  ! (in)

    real(DP), intent(out) :: xyza_HTRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyza_HTRC_RHS_phys(IA,JA,KA,TRC_TOT_NUM)

!    call MessageNotify('M', module_name, "HTRCRHS..")

    call DOGCM_Dyn_spm_HTRCRHS( xyza_HTRC_RHS(IS:IE,JS:JE,KS:KE,:),      &  ! (out)
         & xyza_TRC(IS:IE,JS:JE,KS:KE,:),                                                  &  ! (in)
         & xyz_U(IS:IE,JS:JE,KS:KE), xyz_V(IS:IE,JS:JE,KS:KE), xyz_Div(IS:IE,JS:JE,KS:KE), &  ! (in)
         & xyz_OMG(IS:IE,JS:JE,KS:KE), xyz_H(IS:IE,JS:JE,KS:KE),                           &  ! (in)
         & xyza_HTRC_RHS_phys(IS:IE,JS:JE,KS:KE,:)                                         &  ! (in)
         & )
    
  end subroutine DOGCM_Dyn_driver_HTRCRHS

  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_driver_MOMBarocRHS( xyz_U_RHS, xyz_V_RHS,          &  ! (out)
       & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,                         &  ! (in)
       & xyz_H, xyz_Pres, xyz_DensEdd,                                    &  ! (in)
       & xyz_GeoPot, xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys   &  ! (in)
       & )

    real(DP), intent(out) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(out) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_Vor(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Pres(IA,JA,KA)
    real(DP), intent(in) :: xyz_DensEdd(IA,JA,KA)
    real(DP), intent(in) :: xyz_GeoPot(IA,JA,KA)
    real(DP), intent(in) :: xyz_CoriU(IA,JA,KA)
    real(DP), intent(in) :: xyz_CoriV(IA,JA,KA)
    real(DP), intent(in) :: xyz_URHS_phys(IA,JA,KA)
    real(DP), intent(in) :: xyz_VRHS_phys(IA,JA,KA)
    
!    call MessageNotify('M', module_name, "MOMBarocRHS..")

    call DOGCM_Dyn_spm_MOMBarocRHS( &
         & xyz_U_RHS(IS:IE,JS:JE,KS:KE), xyz_V_RHS(IS:IE,JS:JE,KS:KE),                     &  ! (out)
         & xyz_U(IS:IE,JS:JE,KS:KE), xyz_V(IS:IE,JS:JE,KS:KE), xyz_OMG(IS:IE,JS:JE,KS:KE), &  ! (in)
         & xyz_Vor(IS:IE,JS:JE,KS:KE), xyz_Div(IS:IE,JS:JE,KS:KE),                         &  ! (in)
         & xyz_H(IS:IE,JS:JE,KS:KE),                                                       &  ! (in)
         & xyz_Pres(IS:IE,JS:JE,KS:KE), xyz_DensEdd(IS:IE,JS:JE,KS:KE),                    &  ! (in)
         & xyz_GeoPot(IS:IE,JS:JE,KS:KE),                                                  &  ! (in)
         & xyz_CoriU(IS:IE,JS:JE,KS:KE), xyz_CoriV(IS:IE,JS:JE,KS:KE),                     &  ! (in)
         & xyz_URHS_phys(IS:IE,JS:JE,KS:KE), xyz_VRHS_phys(IS:IE,JS:JE,KS:KE)              &  ! (in)
         & )
    
  end subroutine DOGCM_Dyn_driver_MOMBarocRHS

  !--------------------------------------------------------------

  subroutine DOGCM_Dyn_driver_MOMBarotRHS( xy_UBarot_RHS, xy_VBarot_RHS,  &  ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,                        &  ! (in)
       & xy_UBarocForce, xy_VBarocForce )                                    ! (in)

    real(DP), intent(out) :: xy_UBarot_RHS(IA,JA)
    real(DP), intent(out) :: xy_VBarot_RHS(IA,JA)
    real(DP), intent(in) :: xy_CoriUBarot(IA,JA)
    real(DP), intent(in) :: xy_CoriVBarot(IA,JA)
    real(DP), intent(in) :: xy_SfcPres(IA,JA)    
    real(DP), intent(in) :: xy_UBarocForce(IA,JA)
    real(DP), intent(in) :: xy_VBarocForce(IA,JA)

 !   call MessageNotify('M', module_name, "MOMBarotRHS..")
    
    call DOGCM_Dyn_spm_MOMBarotRHS( &
         & xy_UBarot_RHS(IS:IE,JS:JE), xy_VBarot_RHS(IS:IE,JS:JE),                             &  ! (out)
         & xy_CoriUBarot(IS:IE,JS:JE), xy_CoriVBarot(IS:IE,JS:JE), xy_SfcPres(IS:IE,JS:JE),    &  ! (in)
         & xy_UBarocForce(IS:IE,JS:JE), xy_VBarocForce(IS:IE,JS:JE)                            &  ! (in)
         & )
    
  end subroutine DOGCM_Dyn_Driver_MOMBarotRHS

  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_driver_BarotUpdate( &
       & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,                     & ! (out)
       & xy_Cori, DelTime, DelTimeSSH, PresTAvgCoefA                       & ! (in)
       & )
    
    real(DP), intent(inout) :: xy_UBarotA(IA,JA)
    real(DP), intent(inout) :: xy_VBarotA(IA,JA)
    real(DP), intent(inout) :: xy_SfcPresA(IA,JA)
    real(DP), intent(inout) :: xy_SSHA(IA,JA)
    real(DP), intent(in) :: xy_Cori(IA,JA)
    real(DP), intent(in) :: DelTime
    real(DP), intent(in) :: DelTimeSSH
    real(DP), intent(in) :: PresTAvgCoefA

    call DOGCM_Dyn_spm_BarotUpdate( &
         & xy_UBarotA(IS:IE,JS:JE), xy_VBarotA(IS:IE,JS:JE),               & ! (out)
         & xy_SfcPresA(IS:IE,JS:JE), xy_SSHA(IS:IE,JS:JE),                 & ! (out)
         & xy_Cori(IS:IE,JS:JE), DelTime, DelTimeSSH, PresTAvgCoefA        & ! (in)
         & )
    
  end subroutine DOGCM_Dyn_driver_BarotUpdate
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_driver_OMGDiag( xyz_OMG,      & ! (out)
       & xyz_Div, xyz_H, xyz_HA, DelTime )             ! (in)

    real(DP), intent(out) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_HA(IA,JA,KA)
    real(DP), intent(in) :: DelTime

    call DOGCM_Dyn_spm_OMGDiag( xyz_OMG(IS:IE,JS:JE,KS:KE),         & ! (out)
         & xyz_Div(IS:IE,JS:JE,KS:KE), xyz_H(IS:IE,JS:JE,KS:KE),    & ! (in)
         & xyz_HA(IS:IE,JS:JE,KS:KE), DelTime                       & ! (in)
         & )

  end subroutine DOGCM_Dyn_driver_OMGDiag

  !--------------------------------------------------------------

  subroutine DOGCM_Dyn_driver_VorDivDiag( xyz_Vor, xyz_Div,    & ! (out)
       & xyz_U, xyz_V )                                          ! (in)

    real(DP), intent(out) :: xyz_Vor(IA,JA,KA)
    real(DP), intent(out) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)

    call DOGCM_Dyn_spm_VorDivDiag( &
         & xyz_Vor(IS:IE,JS:JE,KS:KE), xyz_Div(IS:IE,JS:JE,KS:KE),    & ! (out)
         & xyz_U(IS:IE,JS:JE,KS:KE), xyz_V(IS:IE,JS:JE,KS:KE)         & ! (in)
         & )

  end subroutine DOGCM_Dyn_driver_VorDivDiag
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_driver_HydPresDiag( xyz_HydPres, & ! (out)
       & xyz_DensEdd, xyz_H )                             ! (in)

    real(DP), intent(out) :: xyz_HydPres(IA,JA,KA)
    real(DP), intent(in) :: xyz_DensEdd(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)

    call DOGCM_Dyn_spm_HydPresDiag( &
         & xyz_HydPres(IS:IE,JS:JE,KS:KE),                          & ! (out)
         & xyz_DensEdd(IS:IE,JS:JE,KS:KE), xyz_H(IS:IE,JS:JE,KS:KE) & ! (in)
         & )
    
  end subroutine DOGCM_Dyn_driver_HydPresDiag


  subroutine DOGCM_Dyn_driver_UVBarotDiag( xy_UBarot, xy_VBarot,    & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xy_Topo )                       ! (in)

    real(DP), intent(out) :: xy_UBarot(IA,JA)
    real(DP), intent(out) :: xy_VBarot(IA,JA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    real(DP), intent(in) :: xy_Topo(IA,JA)

    call DOGCM_Dyn_spm_UVBarotDiag( &
         & xy_UBarot(IS:IE,JS:JE), xy_VBarot(IS:IE,JS:JE),                               & ! (out)
         & xyz_U(IS:IE,JS:JE,KS:KE), xyz_V(IS:IE,JS:JE,KS:KE), xyz_H(IS:IE,JS:JE,KS:KE), & ! (in)
         & xy_SSH(IS:IE,JS:JE), xy_Topo(IS:IE,JS:JE)                                     & ! (in)
         & )
    
  end subroutine DOGCM_Dyn_driver_UVBarotDiag
  
end module DOGCM_Dyn_driver_mod

  
