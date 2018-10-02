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
       & KinBCTYPE_LinFreeSurf, &
       & KinBCTYPE_FreeSurf
  
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

  !-------------------------------------
  
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

  subroutine DOGCM_Dyn_hspm_vfvm_SSHRHS( xy_SSH_RHS,                   &  ! (out)
       & xy_SSH, xy_TotDepBasic, xy_UBarot, xy_VBarot, xy_FreshWtFlx   )  ! (in)
    
    real(DP), intent(out) :: xy_SSH_RHS(IA,JA)
    real(DP), intent(in) :: xy_SSH(IA,JA)
    real(DP), intent(in) :: xy_TotDepBasic(IA,JA)
    real(DP), intent(in) :: xy_UBarot(IA,JA)
    real(DP), intent(in) :: xy_VBarot(IA,JA)
    real(DP), intent(in) :: xy_FreshWtFlx(IA,JA)    

    select case( KinBC_Surface )
    case (KinBCTYPE_RigidLid)
       xy_SSH_RHS(:,:) = 0d0
    case (KinBCTYPE_LinFreeSurf)       
       call HBEBarot_SSHRHS_LinFreeSfc( xy_SSH_RHS,                   &    ! (out)
            & xy_SSH, xy_TotDepBasic, xy_UBarot, xy_VBarot, xy_FreshWtFlx   )      ! (in)
    case (KinBCTYPE_FreeSurf)       
       call HBEBarot_SSHRHS_NonLinFreeSfc( xy_SSH_RHS,                &    ! (out)
            & xy_SSH, xy_TotDepBasic, xy_UBarot, xy_VBarot, xy_FreshWtFlx   )      ! (in)
    end select
    
  end subroutine DOGCM_Dyn_hspm_vfvm_SSHRHS

  !-------------------------------------

!!$  subroutine DOGCM_Dyn_hspm_vfvm_HTRCRHS( xyza_HTRC_RHS,                             &  ! (out)
!!$       & xyza_TRC, xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H, xyza_HTRC_RHS_phys   )  ! (in)
!!$
!!$    use ProfUtil_mod
!!$    
!!$    real(DP), intent(out) :: xyza_HTRC_RHS(0:iMax-1,jMax,KA,TRC_TOT_NUM)
!!$    real(DP), intent(in) :: xyza_TRC(0:iMax-1,jMax,KA,TRC_TOT_NUM)
!!$    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
!!$    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
!!$    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,KA)
!!$    real(DP), intent(in) :: xyz_OMG(0:iMax-1,jMax,KA)
!!$    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
!!$    real(DP), intent(in) :: xyza_HTRC_RHS_phys(0:iMax-1,jMax,KA,TRC_TOT_NUM)
!!$
!!$    real(DP) :: xyz_HTRC_RHS(0:iMax-1,jMax,KA)
!!$    integer :: n
!!$    
!!$!    call MessageNotify('M', module_name, "HTRCRHS..")
!!$
!!$    call ProfUtil_RapStart('HTRCRHS_driver', 3)
!!$    
!!$    do n = 1, TRC_TOT_NUM 
!!$       call HBEBaroc_HTRCRHS( xyz_HTRC_RHS,                               &  ! (out)
!!$            & xyza_TRC(:,:,:,n), xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H,   &  ! (in)
!!$            & xyza_HTRC_RHS_phys(:,:,:,n)   )                               ! (in)
!!$
!!$       xyza_HTRC_RHS(:,:,:,n) = xyz_HTRC_RHS
!!$    end do
!!$
!!$    call ProfUtil_RapEnd('HTRCRHS_driver', 3)
!!$    
!!$  end subroutine DOGCM_Dyn_hspm_vfvm_HTRCRHS
!!$
  subroutine DOGCM_Dyn_hspm_vfvm_HTRCRHS( xyza_HTRC_RHS,                             &  ! (out)
       & xyza_TRC, xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H, xyza_HTRC_RHS_phys   )  ! (in)

    use ProfUtil_mod
    
    real(DP), intent(out) :: xyza_HTRC_RHS(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyza_TRC(IA,JA,KA,TRC_TOT_NUM)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyza_HTRC_RHS_phys(IA,JA,KA,TRC_TOT_NUM)

    integer :: n
    
!    call MessageNotify('M', module_name, "HTRCRHS..")

    do n = 1, TRC_TOT_NUM 
       call HBEBaroc_HTRCRHS( xyza_HTRC_RHS(:,:,:,n),                     &  ! (out)
            & xyza_TRC(:,:,:,n), xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H,   &  ! (in)
            & xyza_HTRC_RHS_phys(:,:,:,n), KinBC_Surface )                   ! (in)
    end do
    
  end subroutine DOGCM_Dyn_hspm_vfvm_HTRCRHS
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_MOMBarocRHS( xyz_U_RHS, xyz_V_RHS,               &  ! (out)
       & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,                           &  ! (in)
       & xyz_H, xyz_Pres, xyz_DensEdd,                                      &  ! (in)
       & xyz_GeoPot, xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys     &  ! (in)
       & )

    real(DP), intent(out) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(out) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_Vor(IA,JA,KA)    
    real(DP), intent(in) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Pres(IA,JA,KA)
    real(DP), intent(in) :: xyz_DensEdd(IA,JA,KA)
    real(DP), intent(in) :: xyz_GeoPot(IA,JA,KA)
    real(DP), intent(in) :: xyz_CoriU(IA,JA,KA)
    real(DP), intent(in) :: xyz_CoriV(IA,JA,KA)
    real(DP), intent(in) :: xyz_URHS_phys(IA,JA,KA)
    real(DP), intent(in) :: xyz_VRHS_phys(IA,JA,KA)

    call HBEBaroc_MOMRHS_VorDivForm( xyz_U_RHS, xyz_V_RHS,                & ! (out)
         & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,                       & ! (in)
         & xyz_H, xyz_Pres, xyz_DensEdd, xyz_GeoPot,                      & ! (in)
         & xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys             & ! (in)
         & )
         
  end subroutine DOGCM_Dyn_hspm_vfvm_MOMBarocRHS

  !--------------------------------------------------------------

  subroutine DOGCM_Dyn_hspm_vfvm_MOMBarotRHS( xy_UBarot_RHS, xy_VBarot_RHS,  &  ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,                           &  ! (in)
       & xy_UBarocForce, xy_VBarocForce )                                       ! (in)


    real(DP), intent(out) :: xy_UBarot_RHS(IA,JA)
    real(DP), intent(out) :: xy_VBarot_RHS(IA,JA)
    real(DP), intent(in) :: xy_CoriUBarot(IA,JA)
    real(DP), intent(in) :: xy_CoriVBarot(IA,JA)
    real(DP), intent(in) :: xy_SfcPres(IA,JA)    
    real(DP), intent(in) :: xy_UBarocForce(IA,JA)
    real(DP), intent(in) :: xy_VBarocForce(IA,JA)

    call HBEBarot_MOMRHS( xy_UBarot_RHS, xy_VBarot_RHS,  &  ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,       &  ! (in)
       & xy_UBarocForce, xy_VBarocForce )                   ! (in)
    
  end subroutine DOGCM_Dyn_hspm_vfvm_MOMBarotRHS

  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_BarotUpdate( &
       & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,                     & ! (out)
       & xy_Cori, xy_TotDepthBasic, DelTime, DelTimeSSH, PresTAvgCoefA,    & ! (in)
       & xy_FreshWtFlx                                                     & ! (in)
       & )
    
    real(DP), intent(inout) :: xy_UBarotA(IA,JA)
    real(DP), intent(inout) :: xy_VBarotA(IA,JA)
    real(DP), intent(inout) :: xy_SfcPresA(IA,JA)
    real(DP), intent(inout) :: xy_SSHA(IA,JA)
    real(DP), intent(in) :: xy_Cori(IA,JA)
    real(DP), intent(in) :: xy_TotDepthBasic(IA,JA)
    real(DP), intent(in) :: DelTime
    real(DP), intent(in) :: DelTimeSSH
    real(DP), intent(in) :: PresTAvgCoefA
    real(DP), intent(in) :: xy_FreshWtFlx(IA,JA)

    select case( KinBC_Surface )
    case (KinBCTYPE_RigidLid)
       call HBEBarot_Update_LinFreeSfc( &
         & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,                   & ! (out)
         & xy_Cori, xy_TotDepthBasic, DelTime, DelTimeSSH, PresTAvgCoefA,  & ! (in)
         & 0d0, xy_FreshWtFlx                                              & ! (in)
         & )
    case (KinBCTYPE_LinFreeSurf, KinBCTYPE_FreeSurf)       
       call HBEBarot_Update_LinFreeSfc( &
         & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,                   & ! (out)
         & xy_Cori, xy_TotDepthBasic, DelTime, DelTimeSSH, PresTAvgCoefA,  & ! (in)
         & 1d0, xy_FreshWtFlx                                              & ! (in)
         & )
    end select
    
    
  end subroutine DOGCM_Dyn_hspm_vfvm_BarotUpdate
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_OMGDiag( xyz_OMG,         & ! (out)
       & xyz_Div, xyz_H, xyz_HA,                           & ! (in)
       & xyz_Z, xy_Topo, DelTime )                           ! (in)

    real(DP), intent(out) :: xyz_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_HA(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)
    real(DP), intent(in) :: xy_Topo(IA,JA)
    real(DP), intent(in) :: DelTime

    call HBEDiagnose_OMG( xyz_OMG,          & ! (out)
         & xyz_Div, xyz_H, xyz_HA,          & ! (in)
         & xyz_Z, xy_Topo, DelTime )          ! (in)    

  end subroutine DOGCM_Dyn_hspm_vfvm_OMGDiag

  subroutine DOGCM_Dyn_hspm_vfvm_OMGDiag2( xyz_OMG,         & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_HA, DelTime )               ! (in)

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
       & xyz_U, xyz_V )                                             ! (in)

    real(DP), intent(out) :: xyz_Vor(IA,JA,KA)
    real(DP), intent(out) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)

    call HBEDiagnose_VorDiv( xyz_Vor,  xyz_Div,   & ! (out)
         & xyz_U, xyz_V                           & ! (in)
         & )

  end subroutine DOGCM_Dyn_hspm_vfvm_VorDivDiag
  
  !--------------------------------------------------------------
  
  subroutine DOGCM_Dyn_hspm_vfvm_HydPresDiag( xyz_HydPres,    & ! (out)
       & xyz_DensEdd, xyz_H )                                   ! (in)

    real(DP), intent(out) :: xyz_HydPres(IA,JA,KA)
    real(DP), intent(in) :: xyz_DensEdd(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)

    call HBEDiagnose_HydPres( xyz_HydPres,     & ! (out)
       & xyz_DensEdd, xyz_H )                    ! (in)
    
  end subroutine DOGCM_Dyn_hspm_vfvm_HydPresDiag

  !-------------------------------------

  subroutine DOGCM_Dyn_hspm_vfvm_UVBarotDiag( xy_UBarot, xy_VBarot,       & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_Topo )                                     ! (in)

    real(DP), intent(out) :: xy_UBarot(IA,JA)
    real(DP), intent(out) :: xy_VBarot(IA,JA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xy_Topo(IA,JA)

    call HBEDiagnose_UVBarot( xy_UBarot, xy_VBarot,       & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_Topo )                     ! (in)

  end subroutine DOGCM_Dyn_hspm_vfvm_UVBarotDiag
  
end module DOGCM_Dyn_hspm_vfvm_mod
