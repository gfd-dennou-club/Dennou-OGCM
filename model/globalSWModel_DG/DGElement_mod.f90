!-------------------------------------------------------------
! Copyright (c) 2013-2014 Yuta Kawai. All rights reserved.
!-------------------------------------------------------------
!> @brief a template module
!! 
!! @author Yuta Kawai
!!
!!
module DGElement_mod 

  ! モジュール引用; Use statements
  !
  use dc_types, only: DP
  use dc_message, only: MessageNotify

  use VectorSpace_mod
  use LagrangePolyn_mod
  
  ! 宣言文; Declareration statements
  !
  implicit none
  private

  ! 公開手続き
  ! Public procedure
  !

  public :: DGTriElement_Init, DGTriElement_Final

  type, public :: Triangle
     type(vector3d) :: node(3)
  end type Triangle

  type, public :: DGTriElement
     type(vector2d), pointer :: node(:) => null()
     real(DP), pointer :: node_y1(:) => null()
     real(DP), pointer :: node_y2(:) => null()
     type(vector2d), pointer :: sIntNode(:) => null()
     real(DP), pointer :: sIntNode_y1(:) => null()
     real(DP), pointer :: sIntnode_y2(:) => null()
     integer :: nNode
     integer :: nSIntNode
  end type DGTriElement

  ! 非公開手続き
  ! Private procedure
  !

  ! 非公開変数
  ! Private variable
  !
  character(*), parameter:: module_name = 'DGElement_mod' !< Module Name

contains

  !> @brief 
  !!
  !!
  subroutine DGTriElement_Init(DGElem, PolyDegree)
    
    ! 宣言文; Declaration statement
    !
    type(DGTriElement), intent(inout) :: DGElem
    integer, intent(in) :: PolyDegree
    
    ! 局所変数
    ! Local variables
    !
    integer :: Nk, NSIntNode, i

    ! 実行文; Executable statement
    !
    
    Nk = (PolyDegree + 1)*(PolyDegree + 2)/2
    NSIntNode = 3*2**(PolyDegree-1)
    allocate( DGElem%node(Nk) )
    allocate( DGElem%node_y1(Nk), DGElem%node_y2(Nk) )
    allocate( DGElem%sIntNode(NSIntNode) )
    allocate( DGElem%sIntNode_y1(NSIntNode), DGElem%sIntNode_y2(NSIntNode) )

    select case (PolyDegree)
    case(1)
       DGElem%node_y1 = TriNk1_y1; DGElem%node_y2 = TriNk1_y2
       DGElem%sIntNode_y1 = TriNk1_SInt_y1; DGElem%sIntNode_y2 = TriNk1_SInt_y2 
    case(2)
       DGElem%node_y1 = TriNk2_y1; DGElem%node_y2 = TriNk2_y2
       DGElem%sIntNode_y1 = TriNk2_SInt_y1; DGElem%sIntNode_y2 = TriNk2_SInt_y2 
    case(3)
       DGElem%node_y1 = TriNk3_y1; DGElem%node_y2 = TriNk3_y2
       DGElem%sIntNode_y1 = TriNk3_SInt_y1; DGElem%sIntNode_y2 = TriNk3_SInt_y2 
    case default
       call MessageNotify('E', module_name, &
            & "Triangle element associated with the degree '%d' of a polynominal is not supported.", &
            & i=(/ PolyDegree /))
    end select

    !
    DGElem%nNode = Nk
    DGElem%nSIntNode = NSIntNode
    do i=1,Nk
       DGElem%node(i) = (/ DGElem%node_y1(i), DGElem%node_y2(i) /) 
    end do

    do i=1, NSIntNode
       DGElem%sIntNode(i) = (/ DGElem%sIntNode_y1(i), DGElem%sIntNode_y2(i) /)
    end do

  end subroutine DGTriElement_Init

  !> @brief 
  !!
  !!
  subroutine DGTriElement_Final(DGElem)
    
    ! 宣言文; Declaration statement
    !
    type(DGTriElement), intent(inout) :: DGElem
    
    ! 局所変数
    ! Local variables
    !
    
    
    ! 実行文; Executable statement
    !
    
    deallocate(DGElem%node, DGElem%sIntNode)

  end subroutine DGTriElement_Final

end module DGElement_mod

