!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DGCalcusUtil_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP
  use dc_message, only: &
       & MessageNotify

  use LagrangePolyn_mod, only: &
       & TriNk_interpolate, &
       & TriNk_sinteg

  use GridSet_mod, only: &
       & nDGElement, nDGNodePerElem, nDGNodePerFace, &
       & nDGSIntNodePerElem, &
       & get_DGElemJacobian

  use DGHelper_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DGCalcusUtil_Init, DGCalcusUtil_Final
  public :: integrate_over_globalRigion

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DGCalcusUtil_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DGCalcusUtil_Init()

    ! 実行文; Executable statements
    !

  end subroutine DGCalcusUtil_Init

  !>
  !!
  !!
  subroutine DGCalcusUtil_Final()

    ! 実行文; Executable statements
    !

  end subroutine DGCalcusUtil_Final

function integrate_over_globalRigion(wc) result(val)

  real(DP),intent(in) :: wc(:,:)
  real(DP) :: val

  integer :: nk, nc
  real(DP) :: Jacobi(nDGNodePerElem)

  val = 0d0

  !$omp parallel do private(Jacobi) reduction(+:val)
  do nc=1, nDGElement
     Jacobi(:) = get_DGElemJacobian(nc)
     val =  val + TriNk_sinteg(wc(:,nc)*Jacobi(:))
  end do

end function integrate_over_globalRigion

end module DGCalcusUtil_mod

