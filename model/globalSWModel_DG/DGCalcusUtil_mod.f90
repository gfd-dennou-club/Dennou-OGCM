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

  use VectorSpace_mod

  use LagrangePolyn_mod, only: &
       & TriNk_interpolate, &
       & TriNk_sinteg, TriNk_sinteg_dotProdWt

  use GridSet_mod, only: &
       & nDGElement, nDGNodePerElem, nDGNodePerFace, &
       & nDGSIntNodePerElem, &
       & DGElemInfo, &
       & get_DGElemSIntNodeJacobian

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
  real(DP), dimension(nDGSIntNodePerElem) :: &
       & s_Jacobi, s_y1, s_y2

  val = 0d0
  s_y1 = DGElemInfo%sIntNode_y1
  s_y2 = DGElemInfo%sIntNode_y2

  !$omp parallel do private(s_Jacobi) reduction(+:val)
  do nc=1, nDGElement
     s_Jacobi(:) = get_DGElemSIntNodeJacobian(nc)
     val =  val + TriNk_sinteg_dotProdWt( &
          & s_Jacobi(:)*TriNk_interpolate(s_y1, s_y2, wc(:,nc)) &
          & )
  end do

end function integrate_over_globalRigion

end module DGCalcusUtil_mod

