!-------------------------------------------------------------
! Copyright (c) 2017-2017 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Kawai Yuta
!!
!!
module OcnDiag_DiagUtil_mod

  ! モジュール引用; Use statement
  !

  !* gtool5
  
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  !* Dennou-OGCM / SeaIce

  use GridIndex_mod, only: &
       & IS, IE, IA, JS, JE, JA, KS, KE, KA

  use DOGCM_Admin_Grid_mod, only: &
       & z_CDK

  use DOGCM_Admin_Constants_mod, only: &
       & RefDens, Grav
  
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !  
  public :: OcnDiag_DiagUtil_Init, OcnDiag_DiagUtil_Final
  public :: DiagUtil_GetBVFreq
  public :: DiagUtil_GetDensPot

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  
contains

  subroutine OcnDiag_DiagUtil_Init()
  end subroutine OcnDiag_DiagUtil_Init

  subroutine OcnDiag_DiagUtil_Final()
  end subroutine OcnDiag_DiagUtil_Final
  
  subroutine DiagUtil_GetBVFreq( xyz_BVFreq, &
       & xyz_PTemp, xyz_Salt, xyz_H, xyz_Z )
    
    use EOSDriver_mod, only: &
         & EOSDriver_alpha_beta
    
    ! 宣言文; Declareration statements
    !    
    real(DP), intent(out) :: xyz_BVFreq(IA,JA,KA)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyz_H(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)

    ! 局所変数
    ! Local variables
    !
    real(DP) :: xy_alpha(IA,JA)
    real(DP) :: xy_beta(IA,JA)
    real(DP) :: r_m1(KA)
    real(DP) :: r_m2(KA)
    real(DP) :: xyr_PTemp(IA,JA,KA)
    real(DP) :: xyr_Salt(IA,JA,KA)
    integer :: j
    integer :: k

    ! 実行文; Executable statement
    !

    !$omp parallel private(xy_alpha, xy_beta, j, k)
    !$omp do
    do k=KS, KE-1
       r_m1(k) = z_CDK(k+1)/(z_CDK(k) + z_CDK(k+1))
       r_m2(k) = z_CDK(k  )/(z_CDK(k) + z_CDK(k+1))

       xyr_PTemp(:,:,k) = r_m1(k)*xyz_PTemp(:,:,k) + r_m2(k)*xyz_PTemp(:,:,k+1)
       xyr_Salt(:,:,k) = r_m1(k)*xyz_Salt(:,:,k) + r_m2(k)*xyz_Salt(:,:,k+1)       
    end do
    do j = JS, JE
       xyr_PTemp(:,:,KS-1) = xyz_PTemp(:,:,KS)
       xyr_PTemp(:,:,KE  ) = xyz_PTemp(:,:,KE)
       xyr_Salt (:,:,KS-1) = xyz_Salt(:,:,KS)
       xyr_Salt (:,:,KE  ) = xyz_Salt(:,:,KE)       
    end do

    !$omp do
    do k=KS, KE
       call EOSDriver_alpha_beta( alpha=xy_alpha, beta=xy_beta,   & ! (out)
            & theta=xyz_PTemp(:,:,k), S=xyz_Salt(:,:,k),          & ! (in)
            & p=-RefDens*Grav*xyz_Z(:,:,k) )                        ! (in)

       xyz_BVFreq(:,:,k) = - Grav/(xyz_H(:,:,k)*z_CDK(k)) * (          &
            &   xy_alpha(:,:)*(xyr_PTemp(:,:,k) - xyr_PTemp(:,:,k-1))  &
            & - xy_beta (:,:)*(xyr_Salt (:,:,k) - xyr_Salt (:,:,k-1))  &
            & )
    end do
    !$omp end parallel
    
  end subroutine DiagUtil_GetBVFreq
    
  subroutine DiagUtil_GetDensPot( xyz_DensPot,         &
       & xyz_PTemp, xyz_Salt, xyz_Z   )

    use EOSDriver_mod, only: &
       & EOSDriver_Eval

    ! 宣言文; Declareration statements
    !    
    real(DP), intent(out) :: xyz_DensPot(IA,JA,KA)
    real(DP), intent(in) :: xyz_PTemp(IA,JA,KA)
    real(DP), intent(in) :: xyz_Salt(IA,JA,KA)
    real(DP), intent(in) :: xyz_Z(IA,JA,KA)

    real(DP) :: xyz_Pres(IA,JA,KA)
    integer :: k
    
    ! 実行文; Executable statement
    !

    xyz_Pres(:,:,:) = 0d0
    call EOSDriver_Eval( rhoEdd = xyz_DensPot, & ! (out)
         & theta = xyz_PTemp, S = xyz_Salt,    & ! (in)
         & p     = xyz_Pres )                    ! (in)
    
  end subroutine DiagUtil_GetDensPot
  
end module OcnDiag_DiagUtil_mod
