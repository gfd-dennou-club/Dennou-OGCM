module List_i_mod
  !
  use dc_types
  use VectorSpace_mod
  !
  implicit none
  private
  !
  type, public :: List_i
     integer, pointer :: v_(:,:) => null()
     integer :: refCount = 0
     integer :: vHaloSize
  end type List_i
  interface List_Init
    module procedure List_i_Init
  end interface List_Init
  interface List_Final
    module procedure List_i_Final
  end interface List_Final
  interface getHListSize
    module procedure List_i_getHListSize
  end interface getHListSize
  interface getVListSize
     module procedure List_i_getVListSize
  end interface getVListSize
  interface decRef
     module procedure List_i_decRef
  end interface decRef
  interface incRef
    module procedure List_i_incRef
  end interface incRef
  ! Public procedures
  public :: List_Init, List_Final
  public :: incRef, decRef, getHListSize, getVListSize
contains
!
!
subroutine List_i_Init (list, list_vsize, list_hsize, vHaloSize)
  type(List_i), intent(inout) :: list
  integer, intent(in) :: list_vsize, list_hsize
  integer, intent(in) :: vHaloSize
  call List_Final(list)
  allocate(list%v_(1-vHaloSize:list_vsize+vHaloSize, list_hsize))
  list%vHaloSize = vHaloSize
  list%refCount = 1
end subroutine List_i_Init
!
!
subroutine List_i_Final (list)
  type(List_i), intent(inout) :: list
  !
  if(.not. associated(list%v_)) return
  ! Check whether the reference counting for this object is equal to 0.
  ! If so, the resource of list will be release.
  if(list%refCount == 0) then
! write(*,*) "Safe to release a resource of list data.."
    deallocate(list%v_)
  end if
end subroutine List_i_Final
pure function List_i_getHListSize (list) result(listSize)
  type(List_i), intent(in) :: list
  integer :: listSize
  listSize = size(list%v_, 2)
end function List_i_getHListSize
pure function List_i_getVListSize (list) result(listSize)
  type(List_i), intent(in) :: list
  integer :: listSize
  listSize = size(list%v_,1) - 2*list%vHaloSize
end function List_i_getVListSize
!
subroutine List_i_incRef (list)
  type(List_i), intent(inout) :: list
  list%refCount = list%refCount + 1
end subroutine List_i_incRef
subroutine List_i_decRef (list, ret)
  type(List_i), intent(inout) :: list
  logical, optional :: ret
  ! Decrease the reference counting for this object
  list%refCount = list%refCount - 1
  if ( present(ret) ) then
    ret = ( list%refCount == 0 )
  end if
! write(*,*) "current RefCount=", list%refCount
end subroutine List_i_decRef
end module List_i_mod
