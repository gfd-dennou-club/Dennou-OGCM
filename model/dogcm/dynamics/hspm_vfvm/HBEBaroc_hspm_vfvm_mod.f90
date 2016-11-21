!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief The module calculates some terms in hydrostatic boussinesq baroclinic equations.  
!! 
!! @author Yuta Kawai
!!
!!
module HBEBaroc_hspm_vfvm_mod

  ! モジュール引用; Use statements
  !

  !* gtool5
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM

  use DOGCM_Admin_Constants_mod, only: &
       & RPlanet, RefDens
  
  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax, lMax,             &
       & KS, KE, KA,                   &
       & z_CDK, z_RCDK, z_FDK, z_RFDK

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
  public :: HBEBaroc_Init, HBEBaroc_Final
  public :: HBEBaroc_HTRCRHS
  public :: HBEBaroc_MOMRHS_VorDivForm
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HBEBaroc_hspm_vfvm_mod' !< Module Name

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
       & xyz_TRC, xyz_U, xyz_V, xyz_Div, xyr_OMG, xyz_H,    &  ! (in)
       & xyz_HTRCRHS_phys )                                    ! (in)

    
    ! 宣言文; Declaration statement
    !
!!$    real(DP), intent(out) :: wz_HTRC_RHS(lMax,KA)
    real(DP), intent(out) :: xyz_HTRC_RHS(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_TRC(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_HTRCRHS_phys(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyz_DSigTRC(0:iMax-1,jMax,KA)
    real(DP) :: xy_HTRCCosLat(0:iMax-1,jMax)
    real(DP) :: xyr_ADVFlxZ(0:iMax-1,jMax,KA)

    integer :: k
    
    ! 実行文; Executable statements
    !
    
!!$    call calc_ADVFlxZ_Cen2( xyr_ADVFlxZ,    & ! (out)
!!$         & xyr_OMG, xyz_TRC            )      ! (in)   

    call calc_ADVFlxZ_QUICK( xyr_ADVFlxZ,   & ! (out)
         & xyr_OMG, xyz_TRC            )      ! (in)   

    !$omp parallel do private(xy_HTRCCosLat)
    do k=KS, KE
       xy_HTRCCosLat(:,:) = xyz_H(:,:,k)*xy_CosLat(:,:)*xyz_TRC(:,:,k)
       xyz_HTRC_RHS(:,:,k) = &
            & + xy_w( &
            &         - w_AlphaOptr_xy( xyz_U(:,:,k)*xy_HTRCCosLat, xyz_V(:,:,k)*xy_HTRCCosLat )        &
            &         + w_xy( - (xyr_ADVFlxZ(:,:,k-1) - xyr_ADVFlxZ(:,:,k))*z_RCDK(k)                   &
            &                 +  xyz_H(:,:,k)*( + xyz_HTRCRHS_phys(:,:,k))                     )        &
            &   ) 
    end do
    
  end subroutine HBEBaroc_HTRCRHS

  subroutine calc_ADVFlxZ_Cen2( xyr_ADVFlxZ, & ! (out)
       & xyr_OMG, xyz_TRC                    & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyr_ADVFlxZ(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_TRC(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !    
    integer :: k

    ! 実行文; Executable statements
    !
    
    xyr_ADVFlxZ(:,:,KS-1) = 0d0
    !$omp parallel do 
    do k=KS, KE-1
       xyr_ADVFlxZ(:,:,k) = xyr_OMG(:,:,k)*( &
            & z_CDK(k+1)*xyz_TRC(:,:,k) + z_CDK(k)*xyz_TRC(:,:,k+1) &
            & )/(z_CDK(k) + z_CDK(k+1))
    end do
    xyr_ADVFlxZ(:,:,KE) = 0d0    
    
  end subroutine calc_ADVFlxZ_Cen2

  subroutine calc_ADVFlxZ_QUICK( xyr_ADVFlxZ, & ! (out)
       & xyr_OMG, xyz_TRC                     & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyr_ADVFlxZ(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_TRC(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: xyz_Diff(0:iMax-1,jMax,KA)
    real(DP) :: m1(KA)
    real(DP) :: m2(KA)
    real(DP) :: dfac(KA)
    integer :: k

    ! 実行文; Executable statements
    !
    
    do k=KS, KE-1
       m1(k) = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       m2(k) = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))
       dfac(k) = 0.125d0*z_CDK(k+1)*z_CDK(k)
    end do
    
    !$omp parallel
    !$omp do
    do k=KS+1, KE-1
       xyz_Diff(:,:,k) =  (   (xyz_TRC(:,:,k-1) - xyz_TRC(:,:,k  ))*z_RFDK(k-1)     &
            &               - (xyz_TRC(:,:,k  ) - xyz_TRC(:,:,k+1))*z_RFDK(k  )     &
            &             )*0.25d0*z_CDK(k)*(z_RFDK(k-1) + z_RFDK(k))
    end do
    
    !$omp do
    do k=KS+1, KE-2
       xyr_ADVFlxZ(:,:,k) = &
            &  xyr_OMG(:,:,k)*(                                                      &
            &      (m1(k)*xyz_TRC(:,:,k) + m2(k)*xyz_TRC(:,:,k+1))                   &
            &    - dfac(k)*(xyz_Diff(:,:,k) + xyz_Diff(:,:,k+1))                     &
            &  )                                                                     &
            &  + abs(xyr_OMG(:,:,k))*dfac(k)*(xyz_Diff(:,:,k) - xyz_Diff(:,:,k+1))
    end do
    !$omp workshare
    xyr_ADVFlxZ(:,:,KS-1) = 0d0
    xyr_ADVFlxZ(:,:,KE  ) = 0d0    
    xyr_ADVFlxZ(:,:,KS  ) = xyr_OMG(:,:,KS  )*(m1(KS)*xyz_TRC(:,:,KS) + m2(KS)*xyz_TRC(:,:,KS+1))
    xyr_ADVFlxZ(:,:,KE-1) = xyr_OMG(:,:,KE-1)*(m1(KE-1)*xyz_TRC(:,:,KE-1) + m2(KE-1)*xyz_TRC(:,:,KE))
    !$omp end workshare
    !$omp end parallel
    
  end subroutine calc_ADVFlxZ_QUICK
  
  !- RHS of momentum equation ------------------------------------

  subroutine HBEBaroc_MOMRHS_VorDivForm( wz_Vor_RHS, wz_Div_RHS,        & ! (out)
       & xyz_U, xyz_V, xyr_OMG, xyz_Vor, xyz_Div,                       & ! (in)
       & xyz_H, xyz_Pres, xyz_DensEdd, xyz_GeoPot,                      & ! (in)
       & xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys             & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: wz_Vor_RHS(lMax,KA)
    real(DP), intent(out) :: wz_Div_RHS(lMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Vor(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyr_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Pres(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_GeoPot(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_CoriU(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_CoriV(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_URHS_phys(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_VRHS_phys(0:iMax-1,jMax,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_A(0:iMax-1,jMax)
    real(DP) :: xy_B(0:iMax-1,jMax)
    real(DP) :: w_GeoPot(lMax)
    real(DP) :: xyr_VAdvU(0:iMax-1,jMax,KA)
    real(DP) :: xyr_VadvV(0:iMax-1,jMax,KA)

    integer :: k

    ! 実行文; Executable statements
    !

    xyr_VAdvU(:,:,KS-1) = 0d0
    xyr_VAdvV(:,:,KS-1) = 0d0
    !$omp parallel do
    do k = KS, KE-1
       xyr_VAdvU(:,:,k) = xyr_OMG(:,:,k)*(xyz_U(:,:,k) - xyz_U(:,:,k+1))*z_RFDK(k)
       xyr_VAdvV(:,:,k) = xyr_OMG(:,:,K)*(xyz_V(:,:,k) - xyz_V(:,:,k+1))*z_RFDK(k)
    end do
    xyr_VAdvU(:,:,KE)   = 0d0
    xyr_VAdvV(:,:,KE)   = 0d0

    !$omp parallel do private(w_GeoPot, xy_A, xy_B)
    do k=KS, KE
       w_GeoPot(:) = w_xy(xyz_GeoPot(:,:,k))

       xy_A(:,:) = xy_CosLat*( &
            &   xyz_Vor(:,:,k)*xyz_U(:,:,k) + xyz_CoriU(:,:,k)               &
            & + 0.5d0*(xyr_VAdvV(:,:,k) + xyr_VAdvV(:,:,k-1))/xyz_H(:,:,k)   &
            & + xyz_DensEdd(:,:,k)*xy_GradLat_w(w_GeoPot)/(RefDens*RPlanet)  &
            & - xyz_VRHS_phys(:,:,k)                                         &
            & )

       xy_B(:,:) = xy_CosLat*( &
            &   xyz_Vor(:,:,k)*xyz_V(:,:,k) + xyz_CoriV(:,:,k)               &
            & - 0.5d0*(xyr_VAdvU(:,:,k) + xyr_VAdvU(:,:,k-1))/xyz_H(:,:,k)   &
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

end module HBEBaroc_hspm_vfvm_mod

