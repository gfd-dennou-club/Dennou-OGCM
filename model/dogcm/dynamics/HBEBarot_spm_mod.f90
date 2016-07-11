!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module HBEBarot_spm_mod

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
  character(*), parameter:: module_name = 'HBEBarot_spm_mod' !< Module Name


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

  subroutine HBEBarot_SSHRHS_LinFreeSfc( w_SSH_RHS,                              &  ! (out)
       & xy_SSH, xy_TotDepBasic, xy_UCosBarot, xy_VCosBarot, xy_FreshWtFlx )   ! (in)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(out) :: w_SSH_RHS(lMax)
    real(DP), intent(in)  :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in)  :: xy_TotDepBasic(0:iMax-1,jMax)    
    real(DP), intent(in)  :: xy_UCosBarot(0:iMax-1,jMax)
    real(DP), intent(in)  :: xy_VCosBarot(0:iMax-1,jMax)
    real(DP), intent(in)  :: xy_FreshWtFlx(0:iMax-1,jMax)

    ! 作業変数
    ! Work variables

    
    ! 実行文; Executable statements
    !

    w_SSH_RHS(:) = &
         & - w_AlphaOptr_xy( xy_TotDepBasic*xy_UCosBarot, xy_TotDepBasic*xy_VCosBarot ) &
         & + w_xy( xy_FreshWtFlx )

  end subroutine HBEBarot_SSHRHS_LinFreeSfc

  subroutine HBEBarot_SSHRHS_NonLinFreeSfc( w_SSH_RHS,                              &  ! (out)
       & xy_SSH, xy_TotDepBasic, xy_UCosBarot, xy_VCosBarot, xy_FreshWtFlx )   ! (in)

    ! 宣言文; Declaration statement
    !
    
    real(DP), intent(out) :: w_SSH_RHS(lMax)
    real(DP), intent(in)  :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in)  :: xy_TotDepBasic(0:iMax-1,jMax)    
    real(DP), intent(in)  :: xy_UCosBarot(0:iMax-1,jMax)
    real(DP), intent(in)  :: xy_VCosBarot(0:iMax-1,jMax)
    real(DP), intent(in)  :: xy_FreshWtFlx(0:iMax-1,jMax)

    ! 作業変数
    ! Work variables
    
    real(DP) :: xy_TotDep(0:iMax-1,jMax)
    
    ! 実行文; Executable statements
    !

    xy_TotDep(:,:) = xy_TotDepBasic + xy_SSH
    w_SSH_RHS(:) = &
         & - w_AlphaOptr_xy( xy_TotDep*xy_UCosBarot, xy_TotDep*xy_VCosBarot ) &
         & + w_xy( xy_FreshWtFlx )
    
  end subroutine HBEBarot_SSHRHS_NonLinFreeSfc
  
  subroutine HBEBarot_MOMRHS( xy_UBarot_RHS, xy_VBarot_RHS,  &  ! (out)
       & xy_CoriUBarot, xy_CoriVBarot, xy_SfcPres,           &  ! (in)
       & xy_UBarocForce, xy_VBarocForce )                       ! (in)


    real(DP), intent(out) :: xy_UBarot_RHS(0:iMax-1,jMax)
    real(DP), intent(out) :: xy_VBarot_RHS(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_CoriUBarot(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_CoriVBarot(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_SfcPres(0:iMax-1,jMax)    
    real(DP), intent(in) :: xy_UBarocForce(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_VBarocForce(0:iMax-1,jMax)

    real(DP) :: w_SfcPres(lMax)

    w_SfcPres(:) = w_xy(xy_SfcPres)

    xy_UBarot_RHS(:,:) = &
         & + xy_CoriVBarot                                 &
         & - xy_GradLon_w(w_SfcPres) / (RefDens*RPlanet)   &
         & + xy_UBarocForce

    xy_VBarot_RHS(:,:) = &
         & - xy_CoriUBarot                                 &
         & - xy_GradLat_w(w_SfcPres) / (RefDens*RPlanet)   &
         & + xy_VBarocForce

!!$    write(*,*) "Cori=", -xy_CoriUBarot
!!$    write(*,*) "Pres=", -xy_GradLat_w(w_SfcPres)/(RefDens*RPlanet)
!!$    write(*,*) "Fbaroc=", xy_VBarocForce
    
  end subroutine HBEBarot_MOMRHS

  subroutine HBEBarot_Update_LinFreeSfc( &
       & xy_UBarotA, xy_VBarotA, xy_SfcPresA, xy_SSHA,                       & ! (out)
       & xy_Cori, DelTime, DelTimeSSH, PresTAvgCoefA                         & ! (in)
       & )
    
    real(DP), intent(inout) :: xy_UBarotA(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_VBarotA(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_SfcPresA(0:iMax-1,jMax)
    real(DP), intent(inout) :: xy_SSHA(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_Cori(0:iMax-1,jMax)
    real(DP), intent(in) :: DelTime
    real(DP), intent(in) :: DelTimeSSH
    real(DP), intent(in) :: PresTAvgCoefA

    real(DP) :: w_DDiv(lMax)
    real(DP) :: w_DVor(lMax)
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    
    real(DP) :: xy_DUBarot(0:iMax-1,jMax)
    real(DP) :: xy_DVBarot(0:iMax-1,jMax)
    real(DP) :: xy_DSfcPres(0:iMax-1,jMax)

    integer :: itr
    
    call calc_UVCosLat2VorDiv( xy_UBarotA*xy_CosLat, xy_VBarotA*xy_CosLat, & ! (in)
         & w_Vor, w_Div )                                                    ! (out)

    w_DDiv(:) = - w_Div; w_DVor(:) = 0d0
    call calc_VorDiv2UV( w_DVor, w_DDiv,   &  ! (in)
         & xy_DUBarot, xy_DVBarot )           ! (out)

!!$    do itr = 1, 2
!!$       w_DVor(:) = - w_xy( &
!!$            & 1.0d0*DelTime*( &
!!$            & xy_Cori*xy_w(w_DDiv) + 2d0*Omega*xy_DVBarot*xy_CosLat/RPlanet &
!!$            & ) )
!!$       
!!$       call calc_VorDiv2UV(w_DVor, w_DDiv,  & ! (in)
!!$            & xy_DUBarot, xy_DVBarot )        ! (out)
!!$    end do
    
    xy_DSfcPres(:,:) = RefDens * xy_w( RPlanet**2 * w_LaplaInv_w(     &
         & w_Div/(PresTAvgCoefA*DelTimeSSH)                           &
         & ))

    !$omp parallel
    !$omp workshare
    xy_SfcPresA(:,:) = xy_SfcPresA + xy_DSfcPres
    xy_SSHA(:,:) = xy_SfcPresA/(RefDens*Grav)
    xy_UBarotA(:,:) = xy_UBarotA + xy_DUBarot
    xy_VBarotA(:,:) = xy_VBarotA + xy_DVBarot
    !$omp end workshare
    !$omp end parallel

  end subroutine HBEBarot_Update_LinFreeSfc
       
end module HBEBarot_spm_mod
