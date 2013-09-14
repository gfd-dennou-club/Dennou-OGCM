#define _eval1(f, v) f(v)
#define _eval2(f, v1, v2) f(v1, v2)
#define _moduleName(s) s ## _mod
#define moduleName _eval1(_moduleName, ListTypeName)
#define _funcNameSuffix(pref, suff) pref ## _ ## suff
#define funcNameSuffix(suff) _eval2(_funcNameSuffix, ListTypeName, suff)

module moduleName
  
  !
  use dc_types
  use VectorSpace_mod

  !
  implicit none
  private

  !
  type, public :: ListTypeName
     ListElemType, pointer :: v_(:,:) => null()
     integer :: refCount = 0
  end type ListTypeName

  interface List_Init
    module procedure funcNameSuffix(Init)
  end interface List_Init

  interface List_Final
    module procedure funcNameSuffix(Final)
  end interface List_Final
  
  interface getHListSize
    module procedure funcNameSuffix(getHListSize)
  end interface getHListSize

  interface getVListSize
     module procedure funcNameSuffix(getVListSize)
  end interface getVListSize

  interface decRef
     module procedure funcNameSuffix(decRef)
  end interface decRef

  interface incRef
    module procedure funcNameSuffix(incRef)
  end interface incRef
  
  ! Public procedures
  public :: List_Init, List_Final
  public :: incRef, decRef, getHListSize, getVListSize
  
contains

!
!
subroutine funcNameSuffix(Init) (list, list_vsize, list_hsize)
  type(ListTypeName), intent(inout) :: list
  integer, intent(in) :: list_vsize, list_hsize
 
  call List_Final(list)

  allocate(list%v_(list_vsize, list_hsize))
  list%refCount = 1

end subroutine funcNameSuffix(Init)

!
!
subroutine funcNameSuffix(Final) (list)
  type(ListTypeName), intent(inout) :: list

  !
  if(.not. associated(list%v_)) return
    
  ! Check whether the reference counting for this object is equal to 0.
  ! If so, the resource of list will be release. 
  if(list%refCount == 0) then
!    write(*,*) "Safe to release a resource of list data.."
    deallocate(list%v_)
  end if

end subroutine funcNameSuffix(Final)

pure function funcNameSuffix(getHListSize) (list) result(listSize)
  type(ListTypeName), intent(in) :: list
  integer :: listSize

  listSize = size(list%v_, 2)

end function funcNameSuffix(getHListSize)

pure function funcNameSuffix(getVListSize) (list) result(listSize)
  type(ListTypeName), intent(in) :: list
  integer :: listSize

  listSize = size(list%v_, 1)

end function funcNameSuffix(getVListSize)

!
subroutine funcNameSuffix(incRef) (list)
  type(ListTypeName), intent(inout) :: list
 
  list%refCount = list%refCount + 1

end subroutine funcNameSuffix(incRef)

subroutine funcNameSuffix(decRef) (list, ret)
  type(ListTypeName), intent(inout) :: list
  logical, optional :: ret

  ! Decrease the reference counting for this object 
  list%refCount = list%refCount - 1
  
  if ( present(ret) ) then
    ret = ( list%refCount == 0 )
  end if

!  write(*,*) "current RefCount=", list%refCount
end subroutine funcNameSuffix(decRef)

end module moduleName
