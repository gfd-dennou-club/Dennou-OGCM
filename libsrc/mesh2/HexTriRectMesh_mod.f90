!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module HexTriRectMesh_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN, STRING

  use dc_message, only: &
       & MessageNotify

  use VectorSpace_mod
  use SphericalCoord_mod
  use geometry_mod
  use PolyMesh_mod
  use fvMeshInfo_mod
  use GeometricField_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !
  public :: HexTriRectMesh_Init, HexTriRectMesh_Final

  type, public :: HexTriRectMesh
     type(PolyMesh), pointer :: mesh => null()
  end type HexTriRecMesh

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'HexTriRectMesh_mod' !< Module Name

  
contains

  !>
  !!
  !!
  subroutine HexTriRectMesh_Init(self, mesh)

    type(HexTriIcMesh), intent(inout) :: self
    type(PolyMesh), intent(in) :: mesh
    
    ! 実行文; Executable statements
    !
    
    self%mesh => mesh

  end subroutine HexTriRectMesh_Init

  !>
  !!
  !!
  subroutine HexTriRectMesh_Final()

    ! 実行文; Executable statements
    !

  end subroutine HexTriRectMesh_Final

  !> @brief 
  !!
  !!
  subroutine HexTriRectMesh_generare(self, x0, x1, y0, y1, nx, ny)
    
    ! 宣言文; Declaration statement
    !
    type(HexTriRectMesh), intent(inout) :: self
    real(DP), intent(in) :: x0, x1, y0, y2
    integer, intent(in) :: NX, NY

    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
  end subroutine HexTriRectMesh_generare


end module HexTriRectMesh_mod

