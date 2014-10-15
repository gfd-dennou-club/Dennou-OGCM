!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module MeshUsageMask_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: &
       & DP, TOKEN

  use VectorSpace_mod

  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  type, public :: MeshUsageMask
     integer, pointer :: masks(:,:) => null()
     type(vector2d):: screen_x, screen_y
     real(DP) :: DX, DY
  end type MeshUsageMask
  integer, parameter, public :: DEFAULT_MASKNX = 1800
  integer, parameter, public :: DEFAULT_MASKNY = 900

  public :: MeshUsageMask_Init, MeshUsageMask_Final
  public :: MaskToOriCoord, OriToMaskCoord

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'MeshUsageMask_mod' !< Module Name

contains

  !>
  !!
  !!
  subroutine MeshUsageMask_Init(self, &
       & x1, x2, y1, y2, nx, ny )

    type(MeshUsageMask), intent(inout) :: self
    real(DP), intent(in) :: x1, x2, y1, y2
    integer, intent(in) :: nx, ny

    ! 実行文; Executable statements
    !
    
    call MeshUsageMask_Final(self)

    self%screen_x = (/ x1, x2 /)
    self%screen_y = (/ y1, y2 /)
    self%DX = (x2 - x1)/nx
    self%DY = (y2 - y1)/ny
    allocate( self%masks(nx, ny) )

  end subroutine MeshUsageMask_Init

  !>
  !!
  !!
  subroutine MeshUsageMask_Final(self)

    type(MeshUsageMask), intent(inout) :: self

    ! 実行文; Executable statements
    !

    if(associated(self%masks)) deallocate(self%masks)

  end subroutine MeshUsageMask_Final

  !> @brief 
  !!
  !!
  subroutine MaskToOriCoord(self, p1, p2, tX, tY)
    
    ! 宣言文; Declaration statement
    !
    type(MeshUsageMask), intent(in) :: self
    real(DP), intent(out) :: p1, p2
    integer, intent(in) :: tX, tY
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    p1 = self%screen_x%v_(1) + (tX - 0.5d0)*self%DX
    p2 = self%screen_y%v_(1) + (tY - 0.5d0)*self%DY

  end subroutine MaskToOriCoord

  !> @brief 
  !!
  !!
  subroutine OriToMaskCoord(self, tX, tY, p1, p2)
    
    ! 宣言文; Declaration statement
    !
    type(MeshUsageMask), intent(in) :: self
    integer, intent(out) :: tX, tY
    real(DP), intent(in) :: p1, p2

    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    tX = ceiling( (p1 - self%screen_x%v_(1))/self%DX )
    tY = ceiling( (p2 - self%screen_y%v_(1))/self%DY )

  end subroutine OriToMaskCoord


end module MeshUsageMask_mod

