!-------------------------------------------------------------
! Copyright (c) 2013-2016 Kawai Yuta. All rights reserved.
!-------------------------------------------------------------
!> @brief The module which provides some  operators in calculus using spectral method. 
!!
!! 
!! @author Yuta Kawai
!! @since 2016
!!
!!
module VSpmUtil_mod

  ! モジュール引用; Use statements
  !

  !* gtool
  !

  use dc_types, only: DP

  use dc_message, only: &
       & MessageNotify

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: VSpmUtil_Init, VSpmUtil_Final
  public :: VSpm_IntBtmToTop_xyz
  
  interface VSpm_Int_BtmToTop
     module procedure VSpm_IntBtmToTop_xyz
  end interface VSpm_Int_BtmToTop
  public :: VSpm_Int_BtmToTop
  
  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'SpmlUtil_mod' !< Module Name

  integer :: IA
  integer :: JA
  integer :: KS
  integer :: KE
  integer :: KA
  real(DP), allocatable :: CDK(:)

  logical :: isInitialzed = .false.
contains

  subroutine VSpmUtil_Init(KS_, KE_, KA_, CDK_, IA_, JA_)
    integer, intent(in) :: KS_
    integer, intent(in) :: KE_
    integer, intent(in) :: KA_
    real(DP), intent(in) :: CDK_(KA)
    integer, intent(in) :: IA_
    integer, intent(in) :: JA_
    
    KS = KS_
    KE = KE_
    KA = KA_

    IA = IA_
    JA = JA_
    
    allocate(CDK(KA))
    CDK(:) = CDK_

    isInitialzed = .true.
    
  end subroutine VSpmUtil_Init

  !------------------------------

  subroutine VSpmUtil_Final()

    if(isInitialzed) then
       deallocate(CDK)
    end if
    isInitialzed = .false.
    
  end subroutine VSpmUtil_Final

  !------------------------------

  function VSpm_IntBtmToTop_xyz( xyz, xyz_e3 ) result(xy_Int_ret)

    use SpmlUtil_mod, only: g_Sig_Weight
    
    real(DP), intent(in) :: xyz(:,:,:)
    real(DP), intent(in) :: xyz_e3(:,:,:)
    real(DP) :: xy_Int_ret(size(xyz,1), size(xyz,2))

    real(DP) :: xy_Int(size(xyz,1),size(xyz,2))
    integer :: k
    
    xy_Int(:,:) = 0d0
    !$omp parallel do reduction(+: xy_Int)
    do k=KS, KE
       xy_Int(:,:) = xy_Int + xyz(:,:,k)*xyz_e3(:,:,k)*g_Sig_Weight(k)
    end do

    xy_Int_ret(:,:) =  xy_Int
    
  end function VSpm_IntBtmToTop_xyz
  
end module VSpmUtil_mod
