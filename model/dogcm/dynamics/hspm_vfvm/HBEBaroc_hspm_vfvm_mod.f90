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

  use ProfUtil_mod
  
  use DOGCM_Admin_Constants_mod, only: &
       & RPlanet, RefDens
  
  use DOGCM_Admin_Grid_mod, only: &
       & iMax, jMax, lMax,                      &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA,    &
       & z_CDK, z_RCDK, z_FDK, z_RFDK,          &
       & JBLOCK
#include "../../admin/DOGCM_Admin_GaussSpmGridIndexDef.h"

  use SpmlUtil_mod, only: &
       & w_xy, xy_w,                            &
       & xy_AlphaOptr_w, w_AlphaOptr_xy,        &
       & xya_AlphaOptr_wa, wa_AlphaOptr_xya,    &
       & xy_GradLon_w, xy_GradLat_w, w_Lapla_w, &
       & calc_VorDiv2UV,                        &
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
    real(DP), intent(out) :: xyz_HTRC_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_TRC(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyr_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_HTRCRHS_phys(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyr_ADVFlxZ(IA,JA,KA)
    real(DP) :: xy_Tmp(0:iMax-1,jMax)
    real(DP) :: wz_HTRC_RHS(lMax,KA)
    integer :: k
    
    ! 実行文; Executable statements
    !
    
!!$    call calc_ADVFlxZ_Cen2( xyr_ADVFlxZ,    & ! (out)
!!$         & xyr_OMG, xyz_TRC            )      ! (in)   

    call calc_ADVFlxZ_QUICK( xyr_ADVFlxZ,   & ! (out)
         & xyr_OMG, xyz_TRC            )      ! (in)   

    !$omp parallel private(xy_Tmp)
    !$omp do 
    do k=KS, KE
       xy_Tmp(:,:) = xyz_H(IS:IE,JS:JE,k)*xy_CosLat(:,:)*xyz_TRC(IS:IE,JS:JE,k)
       wz_HTRC_RHS(:,k) = - w_AlphaOptr_xy( xyz_U(IS:IE,JS:JE,k)*xy_Tmp, xyz_V(IS:IE,JS:JE,k)*xy_Tmp )
    end do
    !$omp do 
    do k=KS, KE
       wz_HTRC_RHS(:,k) = wz_HTRC_RHS(:,k) +  w_xy( &
            &        - (xyr_ADVFlxZ(IS:IE,JS:JE,k-1) - xyr_ADVFlxZ(IS:IE,JS:JE,k))*z_RCDK(k)   &
            &        +  xyz_H(IS:IE,JS:JE,k)*( + xyz_HTRCRHS_phys(IS:IE,JS:JE,k) )             &
            & )

       xyz_HTRC_RHS(IS:IE,JS:JE,k) = xy_w(wz_HTRC_RHS(:,k))       
    end do
    !$omp end parallel

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
    real(DP), intent(out) :: xyr_ADVFlxZ(IA,JA,KA)
    real(DP), intent(in) :: xyr_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_TRC(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !    
    real(DP) :: xyz_Diff(IA,JA,KA)
    real(DP) :: m1(KA)
    real(DP) :: m2(KA)
    real(DP) :: dfac(KA)

    integer :: i
    integer :: j
    integer :: jj
    integer :: k

    ! 実行文; Executable statements
    !
    
    do k=KS, KE-1
       m1(k) = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       m2(k) = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))
       dfac(k) = 0.125d0*z_CDK(k+1)*z_CDK(k)
    end do
    
    !$omp parallel private(i,j,k)
    !$omp do collapse(2)
    do k=KS+1, KE-1
    do j=JS, JE
    do i=IS, IS+_IMAX_-1
       xyz_Diff(i,j,k) =  (   (xyz_TRC(i,j,k-1) - xyz_TRC(i,j,k  ))*z_RFDK(k-1)     &
            &               - (xyz_TRC(i,j,k  ) - xyz_TRC(i,j,k+1))*z_RFDK(k  )     &
            &             )*0.25d0*z_CDK(k)*(z_RFDK(k-1) + z_RFDK(k))
    end do
    end do
    end do
    
    !$omp do collapse(2)
    do k=KS+1, KE-2
    do j=JS, JE
    do i=IS, IS+_IMAX_-1
       xyr_ADVFlxZ(i,j,k) = &
            &  xyr_OMG(i,j,k)*(                                                      &
            &      (m1(k)*xyz_TRC(i,j,k) + m2(k)*xyz_TRC(i,j,k+1))                   &
            &    - dfac(k)*(xyz_Diff(i,j,k) + xyz_Diff(i,j,k+1))                     &
            &  )                                                                     &
            &  + abs(xyr_OMG(i,j,k))*dfac(k)*(xyz_Diff(i,j,k) - xyz_Diff(i,j,k+1))
    end do
    end do
    end do
    !$omp do
    do j=JS,JE
    do i=IS, IS+_IMAX_-1
       xyr_ADVFlxZ(i,j,KS-1) = 0d0
       xyr_ADVFlxZ(i,j,KE  ) = 0d0    
       xyr_ADVFlxZ(i,j,KS  ) = xyr_OMG(i,j,KS  )*(m1(KS)*xyz_TRC(i,j,KS) + m2(KS)*xyz_TRC(i,j,KS+1))
       xyr_ADVFlxZ(i,j,KE-1) = xyr_OMG(i,j,KE-1)*(m1(KE-1)*xyz_TRC(i,j,KE-1) + m2(KE-1)*xyz_TRC(i,j,KE))
    end do
    end do
     
    !$omp end parallel
    
  end subroutine calc_ADVFlxZ_QUICK
  
  !- RHS of momentum equation ------------------------------------

  subroutine HBEBaroc_MOMRHS_VorDivForm( xyz_U_RHS, xyz_V_RHS,          & ! (out)
       & xyz_U, xyz_V, xyr_OMG, xyz_Vor, xyz_Div,                       & ! (in)
       & xyz_H, xyz_Pres, xyz_DensEdd, xyz_GeoPot,                      & ! (in)
       & xyz_CoriU, xyz_CoriV, xyz_URHS_phys, xyz_VRHS_phys             & ! (in)
       & )

    ! 宣言文; Declaration statement
    !    
    real(DP), intent(out) :: xyz_U_RHS(IA,JA,KA)
    real(DP), intent(out) :: xyz_V_RHS(IA,JA,KA)
    real(DP), intent(in) :: xyz_U(IA,JA,KA)
    real(DP), intent(in) :: xyz_V(IA,JA,KA)
    real(DP), intent(in) :: xyz_Vor(IA,JA,KA)
    real(DP), intent(in) :: xyz_Div(IA,JA,KA)
    real(DP), intent(in) :: xyr_OMG(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Pres(IA,JA,KA)
    real(DP), intent(in) :: xyz_DensEdd(IA,JA,KA)
    real(DP), intent(in) :: xyz_GeoPot(IA,JA,KA)
    real(DP), intent(in) :: xyz_CoriU(IA,JA,KA)
    real(DP), intent(in) :: xyz_CoriV(IA,JA,KA)
    real(DP), intent(in) :: xyz_URHS_phys(IA,JA,KA)
    real(DP), intent(in) :: xyz_VRHS_phys(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xyr_VAdvU(IA,JA,KA)
    real(DP) :: xyr_VadvV(IA,JA,KA)

    real(DP) :: w_Vor_RHS(lMax)
    real(DP) :: w_Div_RHS(lMax)
    real(DP) :: w_GeoPot(lMax)
    
    real(DP) :: xy_A(IA,JA)
    real(DP) :: xy_B(IA,JA)
    real(DP) :: xy_C(IA,JA)
    
    integer :: i
    integer :: j
    integer :: jj
    integer :: k

    ! 実行文; Executable statements
    !

    !$omp parallel private(i,j,jj,k)
    !$omp do collapse(2)
    do k=KS, KE-1
    do j=JS, JE
    do i=IS, IS+_IMAX_-1
       xyr_VAdvU(i,j,k) = xyr_OMG(i,j,k)*(xyz_U(i,j,k) - xyz_U(i,j,k+1))*z_RFDK(k)
       xyr_VAdvV(i,j,k) = xyr_OMG(i,j,K)*(xyz_V(i,j,k) - xyz_V(i,j,k+1))*z_RFDK(k)
    end do
    end do
    end do
    !$omp do
    do j=JS,JE
    do i=IS, IS+_IMAX_-1
       xyr_VAdvU(i,j,KS-1) = 0d0
       xyr_VAdvU(i,j,KE)   = 0d0
       xyr_VAdvV(i,j,KS-1) = 0d0
       xyr_VAdvV(i,j,KE)   = 0d0
    end do
    end do
    !$omp end parallel
    
    !$omp parallel do private(w_Vor_RHS, w_Div_RHS, w_GeoPot, xy_A, xy_B, xy_C, i, j, k)
    do k=KS, KE

       w_GeoPot(:) = w_xy(xyz_GeoPot(IS:IE,JS:JE,k))

       !
       do j=JS,JE
       do i=IS,IS+_IMAX_-1
          xy_A(i,j) =   xyz_Vor(i,j,k)*xyz_U(i,j,k) + xyz_CoriU(i,j,k)               &
               &     +  0.5d0*(xyr_VAdvV(i,j,k) + xyr_VAdvV(i,j,k-1))/xyz_H(i,j,k)   &
               &     -  xyz_VRHS_phys(i,j,k)
       end do
       end do
       xy_A(IS:IE,JS:JE) = xy_CosLat(:,:)*( &
            &   xy_A(IS:IE,JS:JE) +                                                  &
            & + xyz_DensEdd(IS:IE,JS:JE,k)*xy_GradLat_w(w_GeoPot)/(RefDens*RPlanet)  &
            & )

       !
       do j=JS,JE
       do i=IS,IS+_IMAX_-1
          xy_B(i,j) =   xyz_Vor(i,j,k)*xyz_V(i,j,k) + xyz_CoriV(i,j,k)               &
               &     -  0.5d0*(xyr_VAdvU(i,j,k) + xyr_VAdvU(i,j,k-1))/xyz_H(i,j,k)   &
               &     +  xyz_URHS_phys(i,j,k)
       end do
       end do
       xy_B(IS:IE,JS:JE) = xy_CosLat(:,:)*( &
            &   xy_B(IS:IE,JS:JE) +                                                  &
            & + xyz_DensEdd(IS:IE,JS:JE,k)*xy_GradLon_w(w_GeoPot)/(RefDens*RPlanet)  &
            & )


       !
       do j=JS,JE
       do i=IS,IS+_IMAX_-1
          xy_C(i,j) = (  0.5d0*(xyz_U(i,j,k)**2 + xyz_V(i,j,k)**2)                   &
               &       + xyz_Pres(i,j,k)/RefDens )/RPlanet**2
       end do
       end do

       !--
       
       w_Vor_RHS(:) = - w_AlphaOptr_xy(xy_A(IS:IE,JS:JE), xy_B(IS:IE,JS:JE))

       w_Div_RHS(:) = &
            &   w_AlphaOptr_xy(xy_B(IS:IE,JS:JE), - xy_A(IS:IE,JS:JE))              &
            & - w_Lapla_w( w_xy(xy_C(IS:IE,JS:JE)) )

       !---
       
       call calc_VorDiv2UV( w_Vor_RHS(:), w_Div_RHS(:),             & ! (in)
            & xyz_U_RHS(IS:IE,JS:JE,k), xyz_V_RHS(IS:IE,JS:JE,k) )    ! (out)
       
    end do

  end subroutine HBEBaroc_MOMRHS_VorDivForm

end module HBEBaroc_hspm_vfvm_mod

