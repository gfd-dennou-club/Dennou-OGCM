!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DOGCM_Dyn_hspm_vfvm_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM
  use DOGCM_Admin_Grid_mod, only: &
       & IA, IS, IE, IM, &
       & JA, JS, JE, JM, &
       & KA, KS, KE, KM, &
       & iMax, jMax, kMax, lMax

  use SpmlUtil_mod
  
  use DOGCM_Admin_Variable_mod, only: &
       & TRC_TOT_NUM, &
       & TRCID_PTEMP, TRCID_SALT

  use DOGCM_Admin_BC_mod, only: &
       & KinBC_Surface,         &
       & KinBCTYPE_RigidLid,    &
       & KinBCTYPE_LinFreeSurf
  
  use HBEBaroc_hspm_vfvm_mod  
  use HBEBarot_hspm_vfvm_mod
  use HBEDiagnose_hspm_vfvm_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public DOGCM_Dyn_hspm_vfvm_Init, DOGCM_Dyn_hspm_vfvm_Final

  public :: DOGCM_Dyn_hspm_vfvm_SSHRHS
  public :: DOGCM_Dyn_hspm_vfvm_HTRCRHS
  public :: DOGCM_Dyn_hspm_vfvm_MOMBarocRHS
  public :: DOGCM_Dyn_hspm_vfvm_MOMBarotRHS
  public :: DOGCM_Dyn_hspm_vfvm_BarotUpdate
  
  public :: DOGCM_Dyn_hspm_vfvm_OMGDiag
  public :: DOGCM_Dyn_hspm_vfvm_OMGDiag2

  public :: DOGCM_Dyn_hspm_vfvm_HydPresDiag
  public :: DOGCM_Dyn_hspm_vfvm_VorDivDiag
  public :: DOGCM_Dyn_hspm_vfvm_UVBarotDiag

  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DOGCM_Dyn_hspm_vfvm_mod' !< Module Name




contains

  !>
  !!
  !!
  Subroutine DOGCM_Dyn_hspm_vfvm_Init(configNmlName)

    ! 宣言文; Declaration statement
    !
    character(*), intent(in) :: configNmlName

    ! 実行文; Executable statements
    !

!    call read_nmlData(configNmlName)

    call HBEBaroc_Init()    
    call HBEBarot_Init()
    
  end subroutine DOGCM_Dyn_hspm_vfvm_Init

  !>
  !!
  !!
  subroutine DOGCM_Dyn_hspm_vfvm_Final()

    ! 実行文; Executable statements
    !

    call HBEBaroc_Final()    
    call HBEBarot_Final()
    
  end subroutine DOGCM_Dyn_hspm_vfvm_Final

  !-------------------------------------

  !-------------------------------------

  subroutine DOGCM_Dyn_hspm_vfvm_SSHRHS( xy_SSH_RHS,                   &  ! (out)
       & xy_SSH, xy_TotDepBasic, xyz_U, xyz_V, xy_FreshWtFlx   )    ! (in)
    
    real(DP), intent(out) :: xy_SSH_RHS(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_TotDepBasic(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xy_FreshWtFlx(0:iMax-1,jMax)    

    select case( KinBC_Surface )
    case (KinBCTYPE_RigidLid)
       xy_SSH_RHS(:,:) = 0d0
    case (KinBCTYPE_LinFreeSurf)       
       call HBEBarot_SSHRHS_LinFreeSfc( xy_SSH_RHS,                   &    ! (out)
            & xy_SSH, xy_TotDepBasic, xyz_U, xyz_V, xy_FreshWtFlx   )      ! (in)
    end select
    
  end subroutine DOGCM_Dyn_hspm_vfvm_SSHRHS

  !-------------------------------------

  subroutine DOGCM_Dyn_hspm_vfvm_HTRCRHS( xyza_HTRC_RHS,                             &  ! (out)
       & xyza_TRC, xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H, xyza_HTRC_RHS_phys   )  ! (in)

    real(DP), intent(out) :: xyza_HTRC_RHS(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyza_HTRC_RHS_phys(0:iMax-1,jMax,KA,TRC_TOT_NUM)

    real(DP) :: xyz_HTRC_RHS(0:iMax-1,jMax,KA)
    real(DP) :: wz_HTRC_RHS(lMax,KA)
    integer :: n
    
!    call MessageNotify('M', module_name, "HTRCRHS..")

    do n = 1, TRC_TOT_NUM 
!!$       call HBEBaroc_HTRCRHS( wz_HTRC_RHS,                               &  ! (out)
       call HBEBaroc_HTRCRHS( xyz_HTRC_RHS,                               &  ! (out)
            & xyza_TRC(:,:,:,n), xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H,  &  ! (in)
            & xyza_HTRC_RHS_phys(:,:,:,n)   )                               ! (in)

       xyza_HTRC_RHS(:,:,:,n) = xyz_HTRC_RHS
!!$       xyza_HTRC_RHS(:,:,:,n) = xya_wa(wz_HTRC_RHS)
       
    end do
    
  end subroutine DOGCM_Dyn_hspm_vfvm_HTRCRHS

  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_MOMBarocRHS( xyz_U_RHS, xyz_V_RHS,               &  ! (out)
       & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,                           &  ! (in)
       & xyz_H, xyz_Pres, xyz_DensEdd,                                      &  ! (in)
       & xyz_GeoPot, xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys     &  ! (in)
       & )

    real(DP), intent(out) :: xyz_U_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_V_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Vor(0:iMax-1,jMax,KA)    
    real(DP), intent(in) :: xyz_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Pres(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_GeoPot(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_CoriU(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_CoriV(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_URHS_phys(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_VRHS_phys(0:iMax-1,jMax,KA)

    real(DP) :: wz_Vor_RHS(lMax,KA)
    real(DP) :: wz_Div_RHS(lMax,KA)

    call HBEBaroc_MOMRHS_VorDivForm( wz_Vor_RHS, wz_Div_RHS,              & ! (out)
         & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,                       & ! (in)
         & xyz_H, xyz_Pres, xyz_DensEdd, xyz_GeoPot,                      & ! (in)
         & xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys             & ! (in)
         & )

    call calc_VorDiv2UV( wz_Vor_RHS, wz_Div_RHS, & ! (in)
         & xyz_U_RHS, xyz_V_RHS )                  ! (out)
         
  end subroutine DOGCM_Dyn_hspm_vfvm_MOMBarocRHS

  !--------------------------------------------------------------

  subroutine DOGCM_Dyn_hspm_vfvm_MOMBarotRHS( xy_UBarot_RHS, xy_VBarot_RHS,  &  ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,                     &  ! (in)
       & xy_UBarocForce, xy_VBarocForce )                                    ! (in)


    real(DP), intent(out) :: xy_UBarot_RHS(0:iMax-1,jMax)
    real(DP), intent(out) :: xy_VBarot_RHS(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_CoriUBarot(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_CoriVBarot(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_SfcPres(0:iMax-1,jMax)    
    real(DP), intent(in) :: xy_UBarocForce(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_VBarocForce(0:iMax-1,jMax)

    call HBEBarot_MOMRHS( xy_UBarot_RHS, xy_VBarot_RHS,  &  ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,       &  ! (in)
       & xy_UBarocForce, xy_VBarocForce )                   ! (in)
    
  end subroutine DOGCM_Dyn_hspm_vfvm_MOMBarotRHS

  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_BarotUpdate( &
       & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,                     & ! (out)
       & xy_Cori, DelTime, DelTimeSSH, PresTAvgCoefA                       & ! (in)
       & )
    
    real(DP), intent(inout) :: xy_UBarotA(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_VBarotA(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_SfcPresA(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_SSHA(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_Cori(0:iMax-1,jMax)
    real(DP), intent(in) :: DelTime
    real(DP), intent(in) :: DelTimeSSH
    real(DP), intent(in) :: PresTAvgCoefA

    select case( KinBC_Surface )
    case (KinBCTYPE_RigidLid)
       call HBEBarot_Update_LinFreeSfc( &
         & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,        & ! (out)
         & xy_Cori, DelTime, DelTimeSSH, PresTAvgCoefA          & ! (in)
         & )
    case (KinBCTYPE_LinFreeSurf)       
       call HBEBarot_Update_LinFreeSfc( &
         & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,        & ! (out)
         & xy_Cori, DelTime, DelTimeSSH, PresTAvgCoefA          & ! (in)
         & )
    end select
    
    
  end subroutine DOGCM_Dyn_hspm_vfvm_BarotUpdate
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_OMGDiag( xyz_OMG,         & ! (out)
       & xyz_Div, xyz_H, xyz_HA, DelTime )             ! (in)

    real(DP), intent(out) :: xyz_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: DelTime

    call HBEDiagnose_OMG( xyz_OMG,          & ! (out)
       & xyz_Div, xyz_H, xyz_HA, DelTime )    ! (in)    

  end subroutine DOGCM_Dyn_hspm_vfvm_OMGDiag

  subroutine DOGCM_Dyn_hspm_vfvm_OMGDiag2( xyz_OMG,         & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_HA, DelTime )             ! (in)

    real(DP), intent(out) :: xyz_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: DelTime

    call HBEDiagnose_OMG2( xyz_OMG,                & ! (out)
         & xyz_U, xyz_V, xyz_H, xyz_HA, DelTime )    ! (in)    

  end subroutine DOGCM_Dyn_hspm_vfvm_OMGDiag2
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_VorDivDiag( xyz_Vor, xyz_Div,    & ! (out)
       & xyz_U, xyz_V )                                       ! (in)

    real(DP), intent(out) :: xyz_Vor(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)

    call HBEDiagnose_VorDiv( xyz_Vor,  xyz_Div,   & ! (out)
         & xyz_U, xyz_V                           & ! (in)
         & )

  end subroutine DOGCM_Dyn_hspm_vfvm_VorDivDiag
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_HydPresDiag( xyz_HydPres,    & ! (out)
       & xyz_DensEdd, xyz_H )                             ! (in)

    real(DP), intent(out) :: xyz_HydPres(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)

    call HBEDiagnose_HydPres( xyz_HydPres,     & ! (out)
       & xyz_DensEdd, xyz_H )                    ! (in)
    
  end subroutine DOGCM_Dyn_hspm_vfvm_HydPresDiag

  !-------------------------------------

  subroutine DOGCM_Dyn_hspm_vfvm_UVBarotDiag( xy_UBarot, xy_VBarot,       & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xy_Topo )                       ! (in)

    real(DP), intent(out) :: xy_UBarot(0:iMax-1,jMax)
    real(DP), intent(out) :: xy_VBarot(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_Topo(0:iMax-1,jMax)

    call HBEDiagnose_UVBarot( xy_UBarot, xy_VBarot,       & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xy_Topo )             ! (in)

  end subroutine DOGCM_Dyn_hspm_vfvm_UVBarotDiag
  
end module DOGCM_Dyn_hspm_vfvm_mod
