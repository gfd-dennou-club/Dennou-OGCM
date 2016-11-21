!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module HBEDiagnose_hspm_vfvm_mod

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
       & iMax, jMax, lMax,             &
       & KS, KE, KA,                   &
       & z_CDK, z_RCDK, z_FDK, z_RFDK
  
  use SpmlUtil_mod, only: &
       & w_xy, xy_w,                                      &
       & calc_UVCosLat2VorDiv, calc_VorDiv2UV,            &
       & xy_AlphaOptr_w, w_AlphaOptr_xy,                  &
       & xya_AlphaOptr_wa, wa_AlphaOptr_xya,              &
       & xy_GradLon_w, xy_GradLat_w, w_Lapla_w,           &
       & xy_IntSig_BtmToTop_xyz, xyz_IntSig_SigToTop_xyz, &
       & xyz_DSig_xyz,                                    &
       & xy_CosLat
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HBEDiagnose_Init, HBEDiagnose_Final
  public :: HBEDiagnose_OMG
  public :: HBEDiagnose_OMG2
  public :: HBEDiagnose_HydPres
  public :: HBEDiagnose_VorDiv
  public :: HBEDiagnose_UVBarot
  
  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HBEDiagnose_hspm_vfvm_mod' !< Module Name


contains

  !>
  !!
  !!
  Subroutine HBEDiagnose_Init()

    ! 宣言文; Declaration statement
    !

    ! 実行文; Executable statements
    !

  end subroutine HBEDiagnose_Init

  !>
  !!
  !!
  subroutine HBEDiagnose_Final()

    ! 実行文; Executable statements
    !

  end subroutine HBEDiagnose_Final

  !-----------

  subroutine HBEDiagnose_OMG( xyr_OMG,      & ! (out)
       & xyz_Div, xyz_H, xyz_HA, DelTime )    ! (in)


    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: xyr_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: DelTime

    ! 作業変数
    ! Work variables
    
    integer :: i
    integer :: j
    integer :: k

    real(DP) :: z_Div(KA)
    real(DP) :: z_H(KA)
    real(DP) :: r_OMG(KA)

    ! 実行文; Executable statements
    !
    
    !$omp parallel do private(k, r_OMG, z_Div, z_H) collapse(2)
    do j=1, jMax
       do i=0, iMax-1
          z_H(:) = xyz_H(i,j,:)
          z_Div(:) = xyz_Div(i,j,:)

          r_OMG(KS-1) = 0d0
          do k=KS, KE
             r_OMG(k) = r_OMG(k-1) + &
                  & z_H(k)*z_Div(k)*z_CDK(k)
          end do
          xyr_OMG(i,j,:) = r_OMG(:)
       end do
    end do

!!$    write(*,*) "check: OMG:", xyr_OMG(0,:,KE)
!!$    xyr_OMG(:,:,KE) = 0d0
    
  end subroutine HBEDiagnose_OMG

  !-----------

  subroutine HBEDiagnose_OMG2( xyr_OMG,           & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_HA, DelTime )    ! (in)

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: xyr_OMG(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: DelTime

    ! 作業変数
    ! Work variables
    
    integer :: i
    integer :: j
    integer :: k

    real(DP) :: z_U(KA)
    real(DP) :: z_V(KA)
    real(DP) :: z_H(KA)
    real(DP) :: r_UHInt(KA)
    real(DP) :: r_VHInt(KA)
    
    real(DP) :: xyr_UHInt(0:iMax-1,jMax,KA)
    real(DP) :: xyr_VHInt(0:iMax-1,jMax,KA)
    real(DP) :: int_w
    real(DP) :: xy_BtmBCFlx(0:iMax-1,jMax)

    ! 実行文; Executable statements
    !
    
    !$omp parallel do private(z_U, z_V, z_H, r_UHInt, r_VHInt, int_w, k) collapse(2)
    do j=1, jMax
       do i=0, iMax-1
          z_U(:) = xyz_U(i,j,:)
          z_V(:) = xyz_V(i,j,:)
          z_H(:) = xyz_H(i,j,:)

          r_UHInt(KS-1) = 0d0
          r_VHInt(KS-1) = 0d0          
          do k=KS, KE
             int_w = z_CDK(k)*z_H(k)
             r_UHInt(k) = r_UHInt(k-1) + int_w*z_U(k)
             r_VHInt(k) = r_VHInt(k-1) + int_w*z_V(k)
          end do

          xyr_UHInt(i,j,KS-1:KE) = r_UHInt(KS-1:KE)
          xyr_VHInt(i,j,KS-1:KE) = r_VHInt(KS-1:KE)
       end do
    end do

    xy_BtmBCFlx(:,:) = 0d0  ! If the topography at sea bottom exists, this flux must be nonzero. 

    !
    
    xyr_OMG(:,:,KS-1) = 0d0
    !$omp parallel do
    do k=KS, KE
       xyr_OMG(:,:,k) = xy_w( w_AlphaOptr_xy(     &
            & xyr_UHInt(:,:,k)*xy_CosLat,         &
            & xyr_VHInt(:,:,k)*xy_CosLat  ))      &
            & + xy_BtmBCFlx(:,:)
    end do

    write(*,*) "check: V:", xyr_VHInt(0,:,KE)
    write(*,*) "check: OMG:", xyr_OMG(0,:,KE)
    xyr_OMG(:,:,KE) = 0d0

    
  end subroutine HBEDiagnose_OMG2
  
  !--------------
  
  subroutine HBEDiagnose_UVBarot( xy_UBarot, xy_VBarot,       & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xy_Topo )                 ! (in)

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: xy_UBarot(0:iMax-1,jMax)
    real(DP), intent(out) :: xy_VBarot(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_Topo(0:iMax-1,jMax)

    ! 作業変数
    ! Work variables
    
    integer :: k
    real(DP) :: xy_TotDep(0:iMax-1,jMax)

    ! 実行文; Executable statements
    !
    
    xy_TotDep = xy_SSH + xy_Topo
    xy_UBarot(:,:) = 0d0
    xy_VBarot(:,:) = 0d0

    !$omp parallel do reduction(+: xy_UBarot, xy_VBarot)
    do k=KS, KE
       xy_UBarot(:,:) = xy_UBarot + xyz_H(:,:,k)*xyz_U(:,:,k)*z_CDK(k)
       xy_VBarot(:,:) = xy_VBarot + xyz_H(:,:,k)*xyz_V(:,:,k)*z_CDK(k)       
    end do
    xy_UBarot(:,:) = xy_UBarot/xy_TotDep
    xy_VBarot(:,:) = xy_VBarot/xy_TotDep
    
  end subroutine HBEDiagnose_UVBarot

  !---------------

  subroutine HBEDiagnose_VorDiv( xyz_Vor, xyz_Div,       & ! (out)
       & xyz_U, xyz_V )                                    ! (in)

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: xyz_Vor(0:iMax-1,jMax,KA)
    real(DP), intent(out) :: xyz_Div(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,KA)

    ! 作業変数
    ! Work variables
    
    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    integer :: k
    
    !$omp parallel do private(w_Vor, w_Div)
    do k=KS, KE
       call calc_UVCosLat2VorDiv( &
            & xyz_U(:,:,k)*xy_CosLat, xyz_V(:,:,k)*xy_CosLat, & ! (in)
            & w_Vor, w_Div                                    & ! (out)
            & )
       
       xyz_Vor(:,:,k) = xy_w(w_Vor)
       xyz_Div(:,:,k) = xy_w(w_Div)       
    end do

  end subroutine HBEDiagnose_VorDiv

  !---------------
  
  subroutine HBEDiagnose_HydPres( xyz_HydPres, & ! (out)
       & xyz_DensEdd, xyz_H )                    ! (in)

    ! 宣言文; Declaration statement
    !        
    real(DP), intent(out) :: xyz_HydPres(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,KA)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,KA)

    ! 作業変数
    ! Work variables
    
    integer :: i
    integer :: j
    integer :: k

    real(DP) :: r_HydPres(KA)
    real(DP) :: z_H(KA)
    real(DP) :: z_DensEdd(KA)
    
    ! 実行文; Executable statements
    !
    
    !$omp parallel do private(k, z_H, z_DensEdd, r_HydPres) collapse(2)
    do j=1, jMax
       do i=0, iMax-1
          z_H(:) = xyz_H(i,j,:)
          z_DensEdd(:) = xyz_DensEdd(i,j,:)

          r_HydPres(KS-1) = 0d0
          do k=KS, KE
             r_HydPres(k) = r_HydPres(k-1) + &
                  & Grav*z_H(k)*z_DensEdd(k)*z_CDK(k)
          end do
          
          xyz_HydPres(i,j,KS:KE) = 0.5d0*( r_HydPres(KS-1:KE-1) + r_HydPres(KS:KE) ) 
       end do
    end do

  end subroutine HBEDiagnose_HydPres
  
end module HBEDiagnose_hspm_vfvm_mod
