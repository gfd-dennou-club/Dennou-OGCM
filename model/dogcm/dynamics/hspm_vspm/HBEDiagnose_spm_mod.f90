!-------------------------------------------------------------
! Copyright (c) 2016-2016 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module HBEDiagnose_spm_mod

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
  character(*), parameter:: module_name = 'HBEDiagnose_spm_mod' !< Module Name


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

  subroutine HBEDiagnose_OMG( xyz_OMG,      & ! (out)
       & xyz_Div, xyz_H, xyz_HA, DelTime )    ! (in)


    use SpmlUtil_mod, only: &
         & DifMat, IntMat, RMatM1, RMatM2
    
    real(DP), intent(out) :: xyz_OMG(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: DelTime

    integer :: k
    real(DP) :: xyz_IntSig_kernel(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_DSigOMG(0:iMax-1,jMax,0:kMax)

    real(DP) :: AMat(0:kMax,0:kMax)
    real(DP) :: b(0:kMax)
    integer :: IPIV(0:kMax)
    integer :: info


    integer :: i
    integer :: j

    
    !$omp parallel
    !$omp workshare
    xyz_IntSig_kernel(:,:,:) = &
         &   xyz_H*xyz_Div &
         & + (xyz_HA - xyz_H)/DelTime
    !$omp end workshare
    !$omp end parallel
    xyz_OMG(:,:,:) = xyz_IntSig_SigToTop_xyz(xyz_IntSig_kernel)
    xyz_OMG(:,:,kMax) = 0d0

!!$    do j=1, jMax
!!$       do i=0, iMax-1
!!$          AMat(1:kMax,:) = matmul(RMatM1, DifMat)
!!$          AMat(0,:) =  0d0; AMat(0,0) = 1d0
!!$          AMat(kMax,:) =  0d0; AMat(kMax,0) = 1d0
!!$          b(1:kMax-1) = - matmul(RMatM2, xyz_H(i,j,:)*xyz_Div(i,j,:))
!!$          b(0) = 0d0; b(kMax) = 0d0
!!$          call DGESV(kMax+1, 1, AMat, kMax+1, IPIV, b, kMax+1, info)
!!$          xyz_OMG(i,j,:) = b(:)
!!$          xyz_OMG(i,j,kMax) = 0d0
!!$       end do
!!$    end do
!!$    
!!$    xyz_DSigOMG = xyz_DSig_xyz(xyz_OMG)/xyz_H
!!$    do k=0, kMax
!!$       write(*,'(a,i2,3(E18.5e2))') "k=", k, xyz_Div(0,jMax/2,k), xyz_DSigOMG(0,jMax/2,k), xyz_OMG(0,jMax/2,k)
!!$    end do
!!$    write(*,*) "----------"
!!$    stop
  end subroutine HBEDiagnose_OMG

  !-----------

  subroutine HBEDiagnose_OMG2( xyz_OMG,           & ! (out)
       & xyz_U, xyz_V, xyz_H, xyz_HA, DelTime )    ! (in)

    real(DP), intent(out) :: xyz_OMG(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_HA(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: DelTime

    integer :: k

    real(DP) :: xyz_IntSig_U(0:iMax-1,jMax,0:kMax)
    real(DP) :: xyz_IntSig_V(0:iMax-1,jMax,0:kMax)

    xyz_IntSig_U(:,:,:) = xyz_IntSig_SigToTop_xyz(xyz_U)
    xyz_IntSig_V(:,:,:) = xyz_IntSig_SigToTop_xyz(xyz_V)

    xyz_OMG(:,:,0) = 0d0
    !$omp parallel do
    do k=1, kMax
       xyz_OMG(:,:,k) = xy_w(w_AlphaOptr_xy( &
            & xyz_H(:,:,k)*xyz_IntSig_U(:,:,k)*xy_CosLat, &
            & xyz_H(:,:,k)*xyz_IntSig_V(:,:,k)*xy_CosLat ))
    end do
    
!!$    write(*,*) "check: hDIV:", xyz_OMG(0,jMax/2,:)    

  end subroutine HBEDiagnose_OMG2
  
  !--------------
  
  subroutine HBEDiagnose_UVBarot( xy_UBarot, xy_VBarot,       & ! (out)
       & xyz_U, xyz_V, xyz_H, xy_SSH, xy_Topo )                 ! (in)

    real(DP), intent(out) :: xy_UBarot(0:iMax-1,jMax)
    real(DP), intent(out) :: xy_VBarot(0:iMax-1,jMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xy_SSH(0:iMax-1,jMax)
    real(DP), intent(in) :: xy_Topo(0:iMax-1,jMax)

    integer :: k
    real(DP) :: xy_TotDep(0:iMax-1,jMax)

    xy_TotDep = xy_SSH + xy_Topo
    xy_UBarot(:,:) = xy_IntSig_BtmToTop_xyz(xyz_H*xyz_U)/xy_TotDep
    xy_VBarot(:,:) = xy_IntSig_BtmToTop_xyz(xyz_H*xyz_V)/xy_TotDep

  end subroutine HBEDiagnose_UVBarot

  !---------------

  subroutine HBEDiagnose_VorDiv( xyz_Vor, xyz_Div,       & ! (out)
       & xyz_U, xyz_V )                                    ! (in)

    real(DP), intent(out) :: xyz_Vor(0:iMax-1,jMax,0:kMax)
    real(DP), intent(out) :: xyz_Div(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_U(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_V(0:iMax-1,jMax,0:kMax)

    real(DP) :: w_Vor(lMax)
    real(DP) :: w_Div(lMax)
    integer :: k
    
    !$omp parallel do private(w_Vor, w_Div)
    do k=0, kMax
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

    real(DP), intent(out) :: xyz_HydPres(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_DensEdd(0:iMax-1,jMax,0:kMax)
    real(DP), intent(in) :: xyz_H(0:iMax-1,jMax,0:kMax)

    real(DP) :: xyz_IntSig_kernel(0:iMax-1,jMax,0:kMax)

    !$omp parallel
    !$omp workshare
    xyz_IntSig_kernel(:,:,:) = Grav*xyz_H*xyz_DensEdd
    !$omp end workshare
    !$omp end parallel
    xyz_HydPres(:,:,:) = xyz_IntSig_SigToTop_xyz(xyz_IntSig_kernel)
    
  end subroutine HBEDiagnose_HydPres
  
end module HBEDiagnose_spm_mod
