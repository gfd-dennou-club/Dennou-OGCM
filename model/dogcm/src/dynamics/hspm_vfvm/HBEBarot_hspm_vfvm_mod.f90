!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module HBEBarot_hspm_vfvm_mod

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
       & iMax, jMax, lMax,                      &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA,    &
       & z_CDK, z_RCDK, z_FDK, z_RFDK,          &
       & JBLOCK
#include "../../admin/DOGCM_Admin_GaussSpmGridIndexDef.h"

  use SpmlUtil_mod, only: &
       & nm_l, w_xy, xy_w,                      &
       & calc_VorDiv2UV, calc_UVCosLat2VorDiv,  &
       & xy_AlphaOptr_w, w_AlphaOptr_xy,        &
       & xya_AlphaOptr_wa, wa_AlphaOptr_xya,    &
       & xy_GradLon_w, xy_GradLat_w, w_Lapla_w, &
       & w_Lapla_w, w_LaplaInv_w,               &
       & xy_IntSig_BtmToTop_xyz,                &
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
  character(*), parameter:: module_name = 'HBEBarot_hspm_vfvm_mod' !< Module Name


  public :: HBEBarot_Init, HBEBarot_Final
  
  public :: HBEBarot_SSHRHS_LinFreeSfc, HBEBarot_SSHRHS_NonLinFreeSfc
  public :: HBEBarot_MOMRHS
  public :: HBEBarot_Update_LinFreeSfc
  
contains

  !>
  !!
  !!
  Subroutine HBEBarot_Init()
    
    ! 宣言文; Declaration statement
    !

    ! 実行文; Executable statements
    !
    
  end subroutine HBEBarot_Init

  
  !>
  !!
  !!
  subroutine HBEBarot_Final()

    ! 実行文; Executable statements
    !

  end subroutine HBEBarot_Final

  !- RHS of sea surface height equation ------------------------------------

  subroutine HBEBarot_SSHRHS_LinFreeSfc( xy_SSH_RHS,                  &  ! (out)
       & xy_SSH, xy_TotDepBasic, xy_UBarot, xy_VBarot, xy_FreshWtFlx )   ! (in)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(out) :: xy_SSH_RHS(IA,JA)
    real(DP), intent(in)  :: xy_SSH(IA,JA)
    real(DP), intent(in)  :: xy_TotDepBasic(IA,JA)    
    real(DP), intent(in)  :: xy_UBarot(IA,JA)
    real(DP), intent(in)  :: xy_VBarot(IA,JA)
    real(DP), intent(in)  :: xy_FreshWtFlx(IA,JA)

    ! 作業変数
    ! Work variables
    !
    real(DP) :: xy_Tmp(0:iMax-1,jMax)
    integer :: j
    integer :: j_

    
    ! 実行文; Executable statements
    !

    !$omp parallel do private(j_)
    do j=1, jMax
       j_ = JS + j -1
       xy_Tmp(:,j) = xy_TotDepBasic(IS:IE,j_)*xy_CosLat(:,j)
    end do
    
    xy_SSH_RHS(IS:IE,JS:JE) = xy_w( &
         & - w_AlphaOptr_xy( &
         &      xy_Tmp*xy_UBarot(IS:IE,JS:JE),   &
         &      xy_Tmp*xy_VBarot(IS:IE,JS:JE)  ) &
         & + w_xy( xy_FreshWtFlx(IS:IE,JS:JE)                          ) &
         & )

  end subroutine HBEBarot_SSHRHS_LinFreeSfc

  !-------------------------------------
  
  subroutine HBEBarot_SSHRHS_NonLinFreeSfc( xy_SSH_RHS,               &  ! (out)
       & xy_SSH, xy_TotDepBasic, xy_UBarot, xy_VBarot, xy_FreshWtFlx )   ! (in)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(out) :: xy_SSH_RHS(IA,JA)
    real(DP), intent(in)  :: xy_SSH(IA,JA)
    real(DP), intent(in)  :: xy_TotDepBasic(IA,JA)    
    real(DP), intent(in)  :: xy_UBarot(IA,JA)
    real(DP), intent(in)  :: xy_VBarot(IA,JA)
    real(DP), intent(in)  :: xy_FreshWtFlx(IA,JA)

    ! 作業変数
    ! Work variables
    
    real(DP) :: xy_Tmp(0:iMax-1,jMax)
    integer :: j
    integer :: j_
    
    ! 実行文; Executable statements
    !
    
    !$omp parallel do private(j_)
    do j=1, jMax
       j_ = JS + j -1
       xy_Tmp(:,j) = (xy_TotDepBasic(IS:IE,j_) + xy_SSH(IS:IE,j_))*xy_CosLat(:,j)
    end do
    
    xy_SSH_RHS(IS:IE,JS:JE) = xy_w( &
         & - w_AlphaOptr_xy( xy_Tmp*xy_UBarot(IS:IE,JS:JE),     &
         &                   xy_Tmp*xy_VBarot(IS:IE,JS:JE) )    &
         & + w_xy( xy_FreshWtFlx(IS:IE,JS:JE))                  &
         & )
    
  end subroutine HBEBarot_SSHRHS_NonLinFreeSfc

  !-------------------------------------
  
  subroutine HBEBarot_MOMRHS( xy_UBarot_RHS, xy_VBarot_RHS,  &  ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,           &  ! (in)
       & xy_UBarocForce, xy_VBarocForce )                       ! (in)

    ! 宣言文; Declaration statement
    !    

    real(DP), intent(out) :: xy_UBarot_RHS(IA,JA)
    real(DP), intent(out) :: xy_VBarot_RHS(IA,JA)
    real(DP), intent(in) :: xy_CoriUBarot(IA,JA)
    real(DP), intent(in) :: xy_CoriVBarot(IA,JA)
    real(DP), intent(in) :: xy_SfcPres(IA,JA)    
    real(DP), intent(in) :: xy_UBarocForce(IA,JA)
    real(DP), intent(in) :: xy_VBarocForce(IA,JA)

    ! 作業変数
    ! Work variables

    integer :: j
    integer :: j_
    real(DP) :: xy_A(0:iMax-1,jMax)
    real(DP) :: xy_B(0:iMax-1,jMax)
    real(DP) :: w_Vor_RHS(lMax)
    real(DP) :: w_Div_RHS(lMax)

    ! 実行文; Executable statements
    !

    !$omp parallel do private(j_)
    do j=1, jMax
       j_ = JS + j -1
       xy_A(:,j) = xy_CosLat(:,j)*(xy_CoriUBarot(IS:IE,j_) - xy_VBarocForce(IS:IE,j_))
       xy_B(:,j) = xy_CosLat(:,j)*(xy_CoriVBarot(IS:IE,j_) + xy_UBarocForce(IS:IE,j_))
    end do
    
    w_Vor_RHS(:) = - w_AlphaOptr_xy(xy_A, xy_B)
    
    w_Div_RHS(:) =   w_AlphaOptr_xy(xy_B, - xy_A)        &
         &         - w_Lapla_w(w_xy(xy_SfcPres(IS:IE,JS:JE)))/(RefDens*RPlanet**2)

    call calc_VorDiv2UV( w_Vor_RHS, w_Div_RHS,                        & ! (in)
         & xy_UBarot_RHS(IS:IE,JS:JE), xy_VBarot_RHS(IS:IE,JS:JE) )     ! (out)
    
  end subroutine HBEBarot_MOMRHS

  !-------------------------------------
  
  subroutine HBEBarot_Update_LinFreeSfc( &
       & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,   & ! (out)
       & xy_Cori, xy_TotDepthBasic,                      & ! (in)
       & DelTime, DelTimeSSH, PresTAvgCoefA,             & ! (in)
       & vareps1, xy_FreshWtFlx                          & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(inout) :: xy_UBarotA(IA,JA)
    real(DP), intent(inout) :: xy_VBarotA(IA,JA)
    real(DP), intent(inout) :: xy_SfcPresA(IA,JA)
    real(DP), intent(inout) :: xy_SSHA(IA,JA)
    real(DP), intent(in) :: xy_Cori(IA,JA)
    real(DP), intent(in) :: xy_TotDepthBasic(IA,JA)
    real(DP), intent(in) :: DelTime
    real(DP), intent(in) :: DelTimeSSH
    real(DP), intent(in) :: PresTAvgCoefA
    real(DP), intent(in) :: vareps1    ! 0: ridgid lid, 1: linear free surface
    real(DP), intent(in) :: xy_FreshWtFlx(IA,JA)

    ! 作業変数
    ! Work variables
    
    real(DP) :: w_DDiv(lMax)
    real(DP) :: w_DVor(lMax)
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    
    real(DP) :: xy_DUBarot(IA,JA)
    real(DP) :: xy_DVBarot(IA,JA)
    real(DP) :: xy_DSfcPres(IA,JA)
    real(DP) :: w_Coef(lMax)
    real(DP) :: w_Tmp(lMax)
    real(DP) :: w_DSfcPres(lMax)
    
    integer :: itr

    integer :: i
    integer :: j
    integer :: l

    integer :: nm(2)
    real(DP) :: diag
    real(DP) :: TotDepth
    
    ! 実行文; Executable statements
    !


    call calc_UVCosLat2VorDiv( &
         & xy_UBarotA(IS:IE,JS:JE)*xy_CosLat, xy_VBarotA(IS:IE,JS:JE)*xy_CosLat, & ! (in)
         & w_Vor, w_Div )                                                          ! (out)

    !------------------------------------------------------------------------    
!!$    do itr = 1, 2
!!$       w_DVor(:) = - w_xy( &
!!$            & 1d0*DelTime*( &
!!$            & xy_Cori*xy_w(w_DDiv) + 2d0*Omega*xy_DVBarot*xy_CosLat/RPlanet &
!!$            & ) )
!!$       
!!$       call calc_VorDiv2UV(w_DVor, w_DDiv,  & ! (in)
!!$            & xy_DUBarot, xy_DVBarot )        ! (out)
!!$    end do
    !------------------------------------------------------------------------    

    TotDepth = xy_TotDepthBasic(IS,JS)
    diag = vareps1/(Grav * PresTAvgCoefA * DelTimeSSH**2 * TotDepth)
    w_Tmp = - w_Div
    
    if (vareps1 > 0d0) then
       w_Coef(1) = 1d0/diag
       w_Tmp(:) = w_Tmp + w_xy(xy_FreshWtFlx(IS:IE,JS:JE)/TotDepth)
    else
       w_Coef(1) = 0d0
    end if    
    do l=2, lMax
       nm = nm_l(l)
       w_Coef(l) =  1d0/(nm(1)*(nm(1) + 1d0)/RPlanet**2 + diag)
    end do

    do l=1, lMax
       nm = nm_l(l)       
       w_DSfcPres(l) = RefDens * w_Coef(l) * w_Tmp(l) / (PresTAvgCoefA*DelTimeSSH)
       w_DDiv(l) = nm(1)*(nm(1) + 1d0)/RPlanet**2 * w_Coef(l) * w_Tmp(l)
    end do
    xy_DSfcPres(IS:IE,JS:JE) = xy_w(w_DSfcPres)

    ! Calculate the correction of barotropic velocity
    
    w_DVor(:) = 0d0
    call calc_VorDiv2UV( w_DVor, w_DDiv,                                        &  ! (in)
         & xy_DUBarot(IS:IE,JS:JE), xy_DVBarot(IS:IE,JS:JE) )                      ! (out)

    !------------------------------------------------------------------------
    
    !$omp parallel do
    do j=JS, JE
       xy_SfcPresA(:,j) = xy_SfcPresA(:,j) + xy_DSfcPres(:,j)
       xy_SSHA(:,j)     = vareps1*xy_SfcPresA(:,j)/(RefDens*Grav)
       xy_UBarotA(:,j)  = xy_UBarotA(:,j) + xy_DUBarot(:,j)
       xy_VBarotA(:,j)  = xy_VBarotA(:,j) + xy_DVBarot(:,j)
    end do

  end subroutine HBEBarot_Update_LinFreeSfc

  subroutine HBEBarot_Update_FreeSfc( &
       & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,   & ! (out)
       & xy_Cori, xy_TotDepthBasic,                      & ! (in)
       & DelTime, DelTimeSSH, PresTAvgCoefA,             & ! (in)
       & vareps1, xy_FreshWtFlx                          & ! (in)
       & )

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(inout) :: xy_UBarotA(IA,JA)
    real(DP), intent(inout) :: xy_VBarotA(IA,JA)
    real(DP), intent(inout) :: xy_SfcPresA(IA,JA)
    real(DP), intent(inout) :: xy_SSHA(IA,JA)
    real(DP), intent(in) :: xy_Cori(IA,JA)
    real(DP), intent(in) :: xy_TotDepthBasic(IA,JA)
    real(DP), intent(in) :: DelTime
    real(DP), intent(in) :: DelTimeSSH
    real(DP), intent(in) :: PresTAvgCoefA
    real(DP), intent(in) :: vareps1    ! 0: ridgid lid, 1: linear free surface
    real(DP), intent(in) :: xy_FreshWtFlx(IA,JA)

    ! 作業変数
    ! Work variables
    
    real(DP) :: w_DDiv(lMax)
    real(DP) :: w_DVor(lMax)
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    
    real(DP) :: xy_DUBarot(IA,JA)
    real(DP) :: xy_DVBarot(IA,JA)
    real(DP) :: xy_DSfcPres(IA,JA)
    real(DP) :: w_Coef(lMax)
    real(DP) :: w_Tmp(lMax)
    real(DP) :: w_DSfcPres(lMax)
    
    integer :: itr

    integer :: i
    integer :: j
    integer :: l

    integer :: nm(2)
    real(DP) :: diag
    real(DP) :: TotDepth
    
    ! 実行文; Executable statements
    !


    call calc_UVCosLat2VorDiv( &
         & xy_UBarotA(IS:IE,JS:JE)*xy_CosLat, xy_VBarotA(IS:IE,JS:JE)*xy_CosLat, & ! (in)
         & w_Vor, w_Div )                                                          ! (out)

    !------------------------------------------------------------------------    
!!$    do itr = 1, 2
!!$       w_DVor(:) = - w_xy( &
!!$            & 1d0*DelTime*( &
!!$            & xy_Cori*xy_w(w_DDiv) + 2d0*Omega*xy_DVBarot*xy_CosLat/RPlanet &
!!$            & ) )
!!$       
!!$       call calc_VorDiv2UV(w_DVor, w_DDiv,  & ! (in)
!!$            & xy_DUBarot, xy_DVBarot )        ! (out)
!!$    end do
    !------------------------------------------------------------------------    

    TotDepth = xy_TotDepthBasic(IS,JS)
    diag = vareps1/(Grav * PresTAvgCoefA * DelTimeSSH**2 * TotDepth)
    w_Tmp = - w_Div
    
    if (vareps1 > 0d0) then
       w_Coef(1) = 1d0/diag
       w_Tmp(:) = w_Tmp + w_xy(xy_FreshWtFlx(IS:IE,JS:JE)/TotDepth)
    else
       w_Coef(1) = 0d0
    end if    
    do l=2, lMax
       nm = nm_l(l)
       w_Coef(l) =  1d0/(nm(1)*(nm(1) + 1d0)/RPlanet**2 + diag)
    end do

    do l=1, lMax
       nm = nm_l(l)       
       w_DSfcPres(l) = RefDens * w_Coef(l) * w_Tmp(l) / (PresTAvgCoefA*DelTimeSSH)
       w_DDiv(l) = nm(1)*(nm(1) + 1d0)/RPlanet**2 * w_Coef(l) * w_Tmp(l)
    end do
    xy_DSfcPres(IS:IE,JS:JE) = xy_w(w_DSfcPres)

    ! Calculate the correction of barotropic velocity
    
    w_DVor(:) = 0d0
    call calc_VorDiv2UV( w_DVor, w_DDiv,                                        &  ! (in)
         & xy_DUBarot(IS:IE,JS:JE), xy_DVBarot(IS:IE,JS:JE) )                      ! (out)

    !------------------------------------------------------------------------
    
    !$omp parallel do
    do j=JS, JE
       xy_SfcPresA(:,j) = xy_SfcPresA(:,j) + xy_DSfcPres(:,j)
       xy_SSHA(:,j)     = vareps1*xy_SfcPresA(:,j)/(RefDens*Grav)
       xy_UBarotA(:,j)  = xy_UBarotA(:,j) + xy_DUBarot(:,j)
       xy_VBarotA(:,j)  = xy_VBarotA(:,j) + xy_DVBarot(:,j)
    end do

  end subroutine HBEBarot_Update_FreeSfc
  
end module HBEBarot_hspm_vfvm_mod
