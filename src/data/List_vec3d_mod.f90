module List_vec3d_mod
  !
  use dc_types
  use VectorSpace_mod
  !
  implicit none
  private
  !
  type, public :: List_vec3d
     type(vector3d), pointer :: v_(:) => null()
     integer :: refCount = 0
  end type List_vec3d
  interface List_Init
    module procedure List_vec3d_Init
  end interface List_Init
  interface List_Final
    module procedure List_vec3d_Final
  end interface List_Final
  interface getListSize
    module procedure List_vec3d_getListSize
  end interface getListSize
  interface decRef
    module procedure List_vec3d_decRef
  end interface decRef
  interface incRef
    module procedure List_vec3d_incRef
  end interface incRef
  ! Public procedures
  public :: List_Init, List_Final
  public :: incRef, decRef, getListSize
contains
!
!
subroutine List_vec3d_Init (list, list_size)
  type(List_vec3d), intent(inout) :: list
  integer, intent(in) :: list_size
  call List_Final(list)
  allocate(list%v_(list_size))
  list%refCount = 1
end subroutine List_vec3d_Init
!
!
subroutine List_vec3d_Final (list)
  type(List_vec3d), intent(inout) :: list
  !
  if(.not. associated(list%v_)) return
  ! Check whether the reference counting for this object is equal to 0.
  ! If so, the resource of list will be release.
  if(list%refCount == 0) then
    write(*,*) "Safe to release a resource.."
    deallocate(list%v_)
  end if
end subroutine List_vec3d_Final
function List_vec3d_getListSize (list) result(listSize)
  type(List_vec3d), intent(in) :: list
  integer :: listSize
  listSize = size(list%v_(:))
end function List_vec3d_getListSize
!
subroutine List_vec3d_incRef (list)
  type(List_vec3d), intent(inout) :: list
  list%refCount = list%refCount + 1
end subroutine List_vec3d_incRef
subroutine List_vec3d_decRef (list, ret)
  type(List_vec3d), intent(inout) :: list
  logical, optional :: ret
  ! Decrease the reference counting for this object
  list%refCount = list%refCount - 1
  if ( present(ret) ) then
    ret = ( list%refCount == 0 )
  end if
end subroutine List_vec3d_decRef
end module List_vec3d_mod
