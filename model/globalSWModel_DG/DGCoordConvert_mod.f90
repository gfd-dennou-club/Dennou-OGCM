!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DGCoordConvert_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP
  use dc_message, only: &
       & MessageNotify

  use VectorSpace_mod
  use SphericalCoord_mod
  use PolyMesh_mod

  use SimParameters_mod, only: &
       & Radius

  use GridSet_mod, only: &
       & DGElemInfo, nDGElement, nDGNodePerElem, &
       & wc_DGNodePos, &
       & get_DGElemCovariantBasis

  use DGHelper_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: DGCoordConvert_Init, DGCoordConvert_Final
  public :: convert_VecBasis_Local2GeoCoord

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DGCoordConvert_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine DGCoordConvert_Init()

    ! 実行文; Executable statements
    !

  end subroutine DGCoordConvert_Init

  !>
  !!
  !!
  subroutine DGCoordConvert_Final()

    ! 実行文; Executable statements
    !

  end subroutine DGCoordConvert_Final

!!!!!!!!!!!!!!!!!!!
  
  !> @brief 
  !!
  !!
  subroutine convert_VecBasis_local2GeoCoord(wc_U, wc_V, &
       & wc_U1, wc_U2 )
    
    ! 宣言文; Declaration statement
    !
    real(DP), dimension(nDGNodePerElem,nDGElement), intent(out) :: wc_U, wc_V
    real(DP), dimension(nDGNodePerElem,nDGElement), intent(in) :: wc_U1, wc_U2
    
    ! 局所変数
    ! Local variables
    !
    type(vector3d) :: b_1, b_2, cartVel, geoVel
    integer :: nk, nc

    ! 実行文; Executable statement
    !

    !$omp parallel do private(nk,b_1,b_2,geoVel)
    do nc=1,nDGElement
       do nk=1,nDGNodePerElem
          b_1 = get_DGElemCovariantBasis(1,DGElemInfo%node(nk),nc)
          b_2 = get_DGElemCovariantBasis(2,DGElemInfo%node(nk),nc)
          geoVel = CartToSphVec(wc_U1(nk,nc)*b_1 + wc_U2(nk,nc)*b_2, wc_DGNodePos(nk,nc))
          wc_U(nk,nc) = geoVel%v_(1)
          wc_V(nk,nc) = geoVel%v_(2)
       end do
    end do

  end subroutine convert_VecBasis_local2GeoCoord
  

end module DGCoordConvert_mod

