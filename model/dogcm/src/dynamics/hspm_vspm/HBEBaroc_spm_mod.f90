!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module HBEBaroc_spm_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Constants_mod
  
  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax, kMax, lMax

  use SpmlUtil_mod, only: &
       & w_xy, xy_w,                            &
       & xy_AlphaOptr_w, w_AlphaOptr_xy,        &
       & xya_AlphaOptr_wa, wa_AlphaOptr_xya,    &
       & xy_GradLon_w, xy_GradLat_w, w_Lapla_w, &
       & xyz_DSig_xyz,                          &
       & xy_CosLat
  
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
  character(*), parameter:: module_name = 'HBEBaroc_spm_mod' !< Module Name


  public :: HBEBaroc_Init, HBEBaroc_Final
  public :: HBEBaroc_HTRCRHS
  public :: HBEBaroc_MOMRHS_VorDivForm

  
contains

  !>
  !!
  !!
  Subroutine HBEBaroc_Init()

    ! 宣言文; Declaration statement
    !

    ! 実行文; Executable statements
    !

  end subroutine HBEBaroc_Init

  !>
  !!
  !!
  subroutine HBEBaroc_Final()

    ! 実行文; Executable statements
    !

  end subroutine HBEBaroc_Final


!- RHS of tracer equation ------------------------------------

!!$  subroutine HBEBaroc_HTRCRHS( wz_HTRC_RHS,              &  ! (out)
  subroutine HBEBaroc_HTRCRHS( xyz_HTRC_RHS,                &  ! (out)
       & xyz_TRC, xyz_U, xyz_V, xyz_Div, xyz_OMG, xyz_H,    &  ! (in)
       & xyz_HTRCRHS_phys )                                    ! (in)
use SpmlUtil_mod
!!$    real(DP), intent(out) :: wz_HTRC_RHS(lMax,0:kMax)
    real(DP), intent(out) :: xyz_HTRC_RHS(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_TRC(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_OMG(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_HTRCRHS_phys(0:iMax-1,jMax,0:kMax)

    integer :: k

    real(DP) :: xyz_DSigTRC(0:iMax-1,jMax,0:kMax)
!!$    real(DP) :: xyz_DSigOMGTRC(0:iMax-1,jMax,0:kMax)
    real(DP) :: xy_HTRCCosLat(0:iMax-1,jMax)
!!$    real(DP) :: xyz_TRCFlxSig(0:iMax-1,jMax,0:kMax)
!!$    real(DP) :: xy_Tmp(0:iMax-1,jMax)

    xyz_DSigTRC(:,:,:) = xyz_DSig_xyz(xyz_TRC)

!!$    xyz_TRCFlxSig = xyz_TRC*xyz_OMG
!!$    xyz_TRCFlxSig(:,:,0) = 0d0    
!!$    xyz_TRCFlxSig(:,:,kMax) = 0d0
!!$    xyz_DSigOMGTRC(:,:,:) = xyz_DSig_xyz(xyz_TRCFlxSig)

    !$omp parallel do private(xy_HTRCCosLat)
    do k=0, kMax
       xy_HTRCCosLat(:,:) = xyz_H(:,:,k)*xy_CosLat(:,:)*xyz_TRC(:,:,k)
!!$       wz_HTRC_RHS(:,k) = &
!!$            & - w_AlphaOptr_xy( xyz_U(:,:,k)*xy_HTRCCosLat, xyz_V(:,:,k)*xy_HTRCCosLat )    &
!!$            & + w_xy(                                                                       &
!!$            &    - xyz_OMG(:,:,k)*xyz_DSigTRC(:,:,k)                                        &
!!$            &    + xyz_H(:,:,k)*( xyz_Div(:,:,k)*xyz_TRC(:,:,k) + xyz_HTRCRHS_phys(:,:,k))  &
!!$            & )
       xyz_HTRC_RHS(:,:,k) = &
            & + xy_w( &
            &         - w_AlphaOptr_xy( xyz_U(:,:,k)*xy_HTRCCosLat, xyz_V(:,:,k)*xy_HTRCCosLat )          &
            &         + w_xy( - xyz_OMG(:,:,k)*xyz_DSigTRC(:,:,k)                                         &
            &                 + xyz_H(:,:,k)*(xyz_Div(:,:,k)*xyz_TRC(:,:,k) + xyz_HTRCRHS_phys(:,:,k)) )  &
            &   )
!!$       xyz_HTRC_RHS(:,:,k) = &
!!$            & + xy_w( &
!!$            &         - w_AlphaOptr_xy( xyz_U(:,:,k)*xy_HTRCCosLat, xyz_V(:,:,k)*xy_HTRCCosLat )          &
!!$            &         + w_xy( - xyz_DSigOMGTRC(:,:,k)                                                     &
!!$            &                 + xyz_H(:,:,k)*( + xyz_HTRCRHS_phys(:,:,k)) )  &
!!$            &   )
!!$       xyz_HTRC_RHS(:,:,k) = xy_w( w_xy( xyz_H(:,:,k)*(                 &
!!$            & - xyz_U(:,:,k)*xy_GradLon_w(w_xy(xyz_TRC(:,:,k)))/RPlanet &
!!$            & - xyz_V(:,:,k)*xy_GradLat_w(w_xy(xyz_TRC(:,:,k)))/RPlanet &            
!!$            & - xyz_OMG(:,:,k)/xyz_H(:,:,k)*xyz_DSigTRC(:,:,k)          &
!!$            & + xyz_HTRCRHS_phys(:,:,k)                                 &
!!$            & ) ) )


    end do

  end subroutine HBEBaroc_HTRCRHS

  
!- RHS of momentum equation ------------------------------------

  subroutine HBEBaroc_MOMRHS_VorDivForm( wz_Vor_RHS, wz_Div_RHS,        & ! (out)
       & xyz_U, xyz_V, xyz_OMG, xyz_Vor, xyz_Div,                       & ! (in)
       & xyz_H, xyz_Pres, xyz_DensEdd, xyz_GeoPot,                      & ! (in)
       & xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys             & ! (in)
       & )

    real(DP), intent(out) :: wz_Vor_RHS(lMax,0:kMax)
    real(DP), intent(out) :: wz_Div_RHS(lMax,0:kMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Vor(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_OMG(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Pres(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_GeoPot(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_CoriU(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_CoriV(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_URHS_phys(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_VRHS_phys(0:iMax-1,jMax,0:kMax)
    
    integer :: k
    real(DP) :: xy_A(0:iMax-1,jMax)
    real(DP) :: xy_B(0:iMax-1,jMax)
    real(DP) :: w_GeoPot(lMax)
    real(DP) :: xyz_DSigU(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_DSigV(0:iMax-1,jMax,0:kMax)
    
    xyz_DSigU(:,:,:) = xyz_DSig_xyz(xyz_U)
    xyz_DSigV(:,:,:) = xyz_DSig_xyz(xyz_V)
    
    !$omp parallel do private(w_GeoPot, xy_A, xy_B)
    do k=0, kMax
       w_GeoPot(:) = w_xy(xyz_GeoPot(:,:,k))

       xy_A(:,:) = xy_CosLat*( &
            &   xyz_Vor(:,:,k)*xyz_U(:,:,k) + xyz_CoriU(:,:,k)               &
            & + xyz_OMG(:,:,k)/xyz_H(:,:,k) * xyz_DSigV(:,:,k)               &
            & + xyz_DensEdd(:,:,k)*xy_GradLat_w(w_GeoPot)/(RefDens*RPlanet)  &
            & - xyz_VRHS_phys(:,:,k)                                         &
            & )

       xy_B(:,:) = xy_CosLat*( &
            &   xyz_Vor(:,:,k)*xyz_V(:,:,k) + xyz_CoriV(:,:,k)               &
            & - xyz_OMG(:,:,k)/xyz_H(:,:,k) * xyz_DSigU(:,:,k)               &
            & + xyz_DensEdd(:,:,k)*xy_GradLon_w(w_GeoPot)/(RefDens*RPlanet)  &
            & + xyz_URHS_phys(:,:,k)                                         &
            & )

       wz_Vor_RHS(:,k) = - w_AlphaOptr_xy(xy_A, xy_B)

       wz_Div_RHS(:,k) = &
            &   w_AlphaOptr_xy(xy_B, - xy_A)                                 &
            & - w_Lapla_w(w_xy(                                              &
            &       0.5d0*(xyz_U(:,:,k)**2 + xyz_V(:,:,k)**2)                &
            &     + xyz_Pres(:,:,k)/RefDens                                  &
            &   ))/RPlanet**2
    end do

  end subroutine HBEBaroc_MOMRHS_VorDivForm

end module HBEBaroc_spm_mod

